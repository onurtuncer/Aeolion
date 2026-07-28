// Geometry/HandoffContract.h
//
// Parser for the aeolion_geometry.json handoff produced by
// vbat-uav-notebooks -- the full design description: planform stations, CST
// airfoil sections, control surfaces, requested lattice topology, and
// propeller blade geometry, tagged with a content hash of the design.
//
// This is a STRICT parser. Every field is type-checked and every ordering /
// range invariant the downstream solver relies on is asserted at load time,
// so a malformed handoff fails here with a located message rather than
// surfacing later as a silently wrong lattice. Unknown top-level keys are
// tolerated (forward compatibility within a schema major version); malformed
// known keys are not.
//
// The parser does no interpolation, no meshing, and no CST evaluation -- it
// yields the contract verbatim, coefficients included, because those
// coefficients are the intended optimization design variables.
//
// This supersedes the earlier GeometryContract.h, which parsed a draft
// "surfaces[]" schema with an integer schema_version. No released producer
// ever emitted that shape, so it was retired rather than kept alongside.

#pragma once

#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Geometry/AirfoilSection.h"
#include "Aeolion/Geometry/BodyGeometry.h"
#include "Aeolion/Geometry/ControlSurface.h"
#include "Aeolion/Geometry/DuctGeometry.h"
#include "Aeolion/Geometry/MeshTopology.h"
#include "Aeolion/Geometry/PlanformStation.h"
#include "Aeolion/Geometry/Propeller.h"
#include "Aeolion/Geometry/PropulsionSpec.h"
#include "Aeolion/Geometry/WingPlacement.h"
#include "Aeolion/Solver/WingParams.h"

#include <charconv>
#include <cmath>
#include <cstddef>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <nlohmann/json.hpp>

namespace Aeolion::Geometry {

// --- accepted contract shape ----------------------------------------------
inline constexpr int SupportedHandoffSchemaMajor = 1;
inline constexpr std::string_view ExpectedReferenceFrame = "aetherion_body_frd";
inline constexpr std::string_view DesignIdPrefix = "sha256:";
inline constexpr std::size_t DesignIdHexLength = 64;

// Schema bounds on the CST expansion order. Fewer than 4 cannot resolve
// leading-edge radius and boat-tail angle independently of the mid-chord
// shape; more than 8 is outside what the fitter is validated for.
inline constexpr std::size_t MinCstCoefficients = 4;
inline constexpr std::size_t MaxCstCoefficients = 8;

inline constexpr std::size_t MinPlanformStations = 2;
inline constexpr std::size_t MinBladeStations = 2;
inline constexpr std::size_t MinBodyStations = 2;

// The body's stated length must agree with the span its stations actually
// cover; a mismatch means one of the two is stale and there is no way to
// tell which, so the contract is rejected rather than silently trusting one.
inline constexpr double BodyLengthConsistencyTolerance = 1e-6;

inline constexpr double EtaMin = 0.0;
inline constexpr double EtaMax = 1.0;
inline constexpr double EtaEndpointTolerance = 1e-12; // eta[0]==0 / eta[n]==1 slack
inline constexpr double MinAxisNorm = 1e-9; // below this, a stated axis vector is treated as degenerate.

/** Fully parsed, validated content of an aeolion_geometry.json handoff. */
struct HandoffContract {
    std::string SchemaVersion;  ///< As written, e.g. "1.0.0".
    std::string DesignId;       ///< 64 hex chars, "sha256:" prefix stripped.
    std::string ReferenceFrame; ///< e.g. "aetherion_body_frd" (x-fwd, y-right, z-down).

    double Span = 0.0; ///< [m], full tip-to-tip.
    std::vector<PlanformStation> Stations;
    std::vector<AirfoilSection> AirfoilSections;
    std::vector<ControlSurface> ControlSurfaces;
    MeshTopology Mesh;
    PropulsionSpec Propulsion;
    BodyGeometry Body;      ///< Absent before schema 1.4.0; check Body.IsPresent().
    WingPlacement Placement; ///< Absent before schema 1.5.0; check Placement.IsStated.
    DuctGeometry Duct;       ///< Absent before schema 1.8.0; check Duct.IsStated. Disjoint from Body.
    Math::Vec3 MomentReferencePoint{0.0, 0.0, 0.0}; ///< Absent before schema 1.8.0; check MomentReferencePointStated.
    bool MomentReferencePointStated = false;
};

/** Thrown for every rejected contract; the message names the offending field. */
class ContractError : public std::invalid_argument {
public:
    explicit ContractError(const std::string& message) : std::invalid_argument(message) {}
};

namespace Detail {

[[nodiscard]] inline const nlohmann::json& Field(const nlohmann::json& node, const std::string& key,
                                                std::string_view where) {
    if (!node.is_object())
        throw ContractError(std::format("{}: expected an object", where));
    const auto it = node.find(key);
    if (it == node.end())
        throw ContractError(std::format("{}: missing required field '{}'", where, key));
    return *it;
}

[[nodiscard]] inline double Number(const nlohmann::json& node, const std::string& key, std::string_view where) {
    const auto& value = Field(node, key, where);
    if (!value.is_number())
        throw ContractError(std::format("{}.{}: expected a number", where, key));
    return value.get<double>();
}

[[nodiscard]] inline int Integer(const nlohmann::json& node, const std::string& key, std::string_view where) {
    const auto& value = Field(node, key, where);
    if (!value.is_number_integer())
        throw ContractError(std::format("{}.{}: expected an integer", where, key));
    return value.get<int>();
}

[[nodiscard]] inline std::string String(const nlohmann::json& node, const std::string& key, std::string_view where) {
    const auto& value = Field(node, key, where);
    if (!value.is_string())
        throw ContractError(std::format("{}.{}: expected a string", where, key));
    return value.get<std::string>();
}

[[nodiscard]] inline const nlohmann::json& Array(const nlohmann::json& node, const std::string& key,
                                                std::string_view where) {
    const auto& value = Field(node, key, where);
    if (!value.is_array())
        throw ContractError(std::format("{}.{}: expected an array", where, key));
    return value;
}

[[nodiscard]] inline const nlohmann::json& Object(const nlohmann::json& node, const std::string& key,
                                                 std::string_view where) {
    const auto& value = Field(node, key, where);
    if (!value.is_object())
        throw ContractError(std::format("{}.{}: expected an object", where, key));
    return value;
}

inline void RequirePositive(double value, std::string_view what) {
    if (!(value > 0.0))
        throw ContractError(std::format("{} must be positive, got {}", what, value));
}

inline void RequireInRange(double value, double low, double high, std::string_view what) {
    if (!(value >= low && value <= high))
        throw ContractError(std::format("{} must lie in [{}, {}], got {}", what, low, high, value));
}

// Coefficient vectors of length n produce an order-(n-1) Bernstein expansion;
// both surfaces of a section must share the same order for the section to be
// a well-formed design-variable block.
[[nodiscard]] inline std::vector<double> Coefficients(const nlohmann::json& node, const std::string& key,
                                                     std::string_view where) {
    const auto& array = Array(node, key, where);
    if (array.size() < MinCstCoefficients || array.size() > MaxCstCoefficients)
        throw ContractError(std::format("{}.{}: expected {} to {} CST coefficients, got {}", where, key,
                                        MinCstCoefficients, MaxCstCoefficients, array.size()));
    std::vector<double> coefficients;
    coefficients.reserve(array.size());
    for (std::size_t i = 0; i < array.size(); ++i) {
        if (!array[i].is_number())
            throw ContractError(std::format("{}.{}[{}]: expected a number", where, key, i));
        coefficients.push_back(array[i].get<double>());
    }
    return coefficients;
}

// Every eta-indexed list must cover the whole semi-span exactly once, in
// order: strictly increasing, rooted at 0, reaching 1. Interpolation between
// stations is meaningless otherwise, and a duplicated eta would make it
// ambiguous which entry wins.
inline void RequireSpanningEtaSequence(const std::vector<double>& etas, std::string_view where) {
    for (std::size_t i = 0; i < etas.size(); ++i) {
        RequireInRange(etas[i], EtaMin, EtaMax, std::format("{}[{}].eta", where, i));
        if (i > 0 && !(etas[i] > etas[i - 1]))
            throw ContractError(std::format("{}: eta must strictly increase, but [{}]={} follows [{}]={}", where, i,
                                            etas[i], i - 1, etas[i - 1]));
    }
    if (std::abs(etas.front() - EtaMin) > EtaEndpointTolerance)
        throw ContractError(std::format("{}: must start at eta={}, got {}", where, EtaMin, etas.front()));
    if (std::abs(etas.back() - EtaMax) > EtaEndpointTolerance)
        throw ContractError(std::format("{}: must reach eta={}, got {}", where, EtaMax, etas.back()));
}

[[nodiscard]] inline int SchemaMajor(const std::string& version) {
    // "MAJOR.MINOR.PATCH", each component a non-negative integer.
    std::size_t start = 0;
    int components[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
        const std::size_t dot = version.find('.', start);
        const bool last = (i == 2);
        if (last != (dot == std::string::npos))
            throw ContractError(std::format("schema_version: expected MAJOR.MINOR.PATCH, got '{}'", version));
        const std::size_t end = last ? version.size() : dot;
        const auto* first = version.data() + start;
        const auto* stop = version.data() + end;
        const auto parsed = std::from_chars(first, stop, components[i]);
        if (parsed.ec != std::errc{} || parsed.ptr != stop)
            throw ContractError(std::format("schema_version: expected MAJOR.MINOR.PATCH, got '{}'", version));
        start = end + 1;
    }
    return components[0];
}

[[nodiscard]] inline std::string DesignIdDigest(const std::string& raw) {
    if (!raw.starts_with(DesignIdPrefix))
        throw ContractError(std::format("design_id: expected a '{}' digest, got '{}'", DesignIdPrefix, raw));
    std::string digest = raw.substr(DesignIdPrefix.size());
    if (digest.size() != DesignIdHexLength)
        throw ContractError(std::format("design_id: expected {} hex characters, got {}", DesignIdHexLength,
                                        digest.size()));
    for (const char c : digest) {
        const bool hex = (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
        if (!hex) throw ContractError(std::format("design_id: '{}' is not a lowercase hex digest", digest));
    }
    return digest;
}

// Reads a 3-element numeric array field into a Vec3, unnormalized. Shared by
// every axis-shaped field the contract carries (hinge axes, the propeller's
// rotation axis) so the "3 components, all numeric" check is written once.
[[nodiscard]] inline Math::Vec3 Vec3Array(const nlohmann::json& node, const std::string& key,
                                          std::string_view where) {
    const auto& array = Array(node, key, where);
    if (array.size() != 3)
        throw ContractError(std::format("{}.{}: expected 3 components, got {}", where, key, array.size()));
    for (std::size_t i = 0; i < 3; ++i)
        if (!array[i].is_number())
            throw ContractError(std::format("{}.{}[{}]: expected a number", where, key, i));
    return Math::Vec3(array[0].get<double>(), array[1].get<double>(), array[2].get<double>());
}

[[nodiscard]] inline Math::Vec3 HingeAxis(const nlohmann::json& node, std::string_view where) {
    const Math::Vec3 axis = Vec3Array(node, "hinge_axis", where);
    if (axis.Norm() < MinAxisNorm)
        throw ContractError(std::format("{}.hinge_axis: must not be the zero vector", where));
    return axis.Normalized(); // direction is all that matters; magnitude is not part of the contract
}

// The solver trails its wake along +x and has no relaxed-wake iteration, so
// "relaxed" is refused outright. Accepting it and solving a frozen wake
// anyway would answer a question the caller did not ask, and nothing in the
// result would reveal the substitution.
inline void RequireSupportedWakeModel(const std::string& name, std::string_view where) {
    if (name == "frozen") return;
    if (name == "relaxed")
        throw ContractError(std::format(
            "{}.wake_model: 'relaxed' is not implemented; this solver trails a frozen wake only", where));
    throw ContractError(std::format("{}.wake_model: unknown model '{}'", where, name));
}

[[nodiscard]] inline ControlSurfaceBinding ParseBindingField(const std::string& surface, std::string_view where) {
    if (surface == "wing") return ControlSurfaceBinding::Wing;
    if (surface == "duct_jet") return ControlSurfaceBinding::DuctJet;
    throw ContractError(std::format("{}.surface: unknown binding '{}'", where, surface));
}

// Schema 1.0.0 has no 'surface' field, so the binding is carried implicitly by
// the control's name (see ControlSurfaceBinding). Codifying the mapping here
// rather than leaving it to each call site is what keeps duct-jet vanes out of
// the wing lattice; ADR-0016 states the identical rule on the producer side.
//
// An unrecognized name is rejected on purpose. There is no safe default: the
// two bindings put a surface on different bodies and give its eta band
// different meanings, so guessing would bind a new control to the wrong
// geometry silently. Emitting an explicit 'surface' field (schema >= 1.1.0)
// overrides this table and is the escape hatch for new names.
[[nodiscard]] inline ControlSurfaceBinding BindingFromName(const std::string& name, std::string_view where) {
    if (name == "aileron" || name == "rudder") return ControlSurfaceBinding::Wing;
    if (name == "vane") return ControlSurfaceBinding::DuctJet;
    throw ContractError(std::format(
        "{}.name: '{}' has no known surface binding; emit an explicit 'surface' field", where, name));
}

// Travel limits arrived in schema 1.4.0. Older documents omit them, leaving
// the limits zeroed -- a caller must therefore check for a degenerate range
// rather than assume every contract states one.
[[nodiscard]] inline DeflectionLimits ParseDeflectionLimits(const nlohmann::json& node, std::string_view where) {
    DeflectionLimits limits;
    if (!node.contains("deflection_limits_deg")) return limits;

    const auto& block = Object(node, "deflection_limits_deg", where);
    const auto blockWhere = std::format("{}.deflection_limits_deg", where);
    limits.MinDeg = Number(block, "min", blockWhere);
    limits.MaxDeg = Number(block, "max", blockWhere);
    if (!(limits.MaxDeg > limits.MinDeg))
        throw ContractError(std::format("{}: max ({}) must exceed min ({})", blockWhere, limits.MaxDeg,
                                        limits.MinDeg));

    if (node.contains("deflection_soft_limit_deg")) {
        limits.HasSoftLimit = true;
        limits.SoftLimitDeg = Number(node, "deflection_soft_limit_deg", where);
        RequirePositive(limits.SoftLimitDeg, std::format("{}.deflection_soft_limit_deg", where));
        // A soft limit outside the mechanical stops is not a conservative
        // bound at all -- it would authorize travel the hardware cannot make.
        if (limits.SoftLimitDeg > limits.MaxDeg || -limits.SoftLimitDeg < limits.MinDeg)
            throw ContractError(std::format(
                "{}.deflection_soft_limit_deg: +/-{} lies outside the hard limits [{}, {}]", where,
                limits.SoftLimitDeg, limits.MinDeg, limits.MaxDeg));
    }
    return limits;
}

// The body of revolution (schema 1.4.0+). Stations run nose to tail, which
// in the contract's x-forward frame means x strictly DECREASING.
[[nodiscard]] inline BodyGeometry ParseBody(const nlohmann::json& node) {
    BodyGeometry body;
    body.Length = Number(node, "length", "body");
    RequirePositive(body.Length, "body.length");

    const auto& stations = Array(node, "stations", "body");
    if (stations.size() < MinBodyStations)
        throw ContractError(std::format("body.stations: expected at least {} stations, got {}", MinBodyStations,
                                        stations.size()));

    for (std::size_t i = 0; i < stations.size(); ++i) {
        const auto where = std::format("body.stations[{}]", i);
        BodyStation station;
        station.x = Number(stations[i], "x", where);
        station.Radius = Number(stations[i], "radius", where);
        // Zero is legitimate at a sharp nose or a closed tail; negative is not.
        if (station.Radius < 0.0)
            throw ContractError(std::format("{}.radius must not be negative, got {}", where, station.Radius));
        if (i > 0 && !(station.x < body.Stations.back().x))
            throw ContractError(std::format("{}: x must strictly decrease nose to tail, but {} follows {}", where,
                                            station.x, body.Stations.back().x));
        body.Stations.push_back(station);
    }

    const double covered = body.Stations.front().x - body.Stations.back().x;
    if (std::fabs(covered - body.Length) > BodyLengthConsistencyTolerance * body.Length)
        throw ContractError(std::format(
            "body.length {} disagrees with the {} its stations span; one of the two is stale", body.Length,
            covered));
    return body;
}

// The duct (schema 1.8.0+): a single annular ring, disjoint from the main
// body -- see DuctGeometry.h for why this replaced the open-tail convention.
[[nodiscard]] inline DuctGeometry ParseDuct(const nlohmann::json& node) {
    DuctGeometry duct;
    duct.IsStated = true;
    duct.Chord = Number(node, "chord", "duct");
    RequirePositive(duct.Chord, "duct.chord");
    duct.InnerDiameter = Number(node, "inner_diameter", "duct");
    RequirePositive(duct.InnerDiameter, "duct.inner_diameter");
    duct.OuterDiameter = Number(node, "outer_diameter", "duct");
    RequirePositive(duct.OuterDiameter, "duct.outer_diameter");
    if (!(duct.OuterDiameter > duct.InnerDiameter))
        throw ContractError(std::format("duct.outer_diameter ({}) must exceed duct.inner_diameter ({})",
                                        duct.OuterDiameter, duct.InnerDiameter));

    const auto& placement = Object(node, "placement", "duct");
    const auto& center = Object(placement, "center", "duct.placement");
    duct.Center = Math::Vec3(Number(center, "x", "duct.placement.center"),
                             Number(center, "y", "duct.placement.center"),
                             Number(center, "z", "duct.placement.center"));
    return duct;
}

} // namespace Detail

/** Parse and strictly validate an already-loaded handoff document; throws ContractError on any violation. */
[[nodiscard]] inline HandoffContract ParseHandoff(const nlohmann::json& root) {
    using namespace Detail;

    HandoffContract contract;

    contract.SchemaVersion = String(root, "schema_version", "$");
    if (const int major = SchemaMajor(contract.SchemaVersion); major != SupportedHandoffSchemaMajor)
        throw ContractError(std::format("schema_version '{}': unsupported major version {}, expected {}",
                                        contract.SchemaVersion, major, SupportedHandoffSchemaMajor));

    // The contract carries no unit metadata beyond this block, so anything
    // other than metres/degrees would silently rescale the whole geometry.
    const auto& units = Object(root, "units", "$");
    if (const auto length = String(units, "length", "units"); length != "m")
        throw ContractError(std::format("units.length: expected 'm', got '{}'", length));
    if (const auto angle = String(units, "angle", "units"); angle != "deg")
        throw ContractError(std::format("units.angle: expected 'deg', got '{}'", angle));

    contract.ReferenceFrame = String(root, "reference_frame", "$");
    if (contract.ReferenceFrame != ExpectedReferenceFrame)
        throw ContractError(std::format("reference_frame: expected '{}', got '{}'", ExpectedReferenceFrame,
                                        contract.ReferenceFrame));

    contract.DesignId = DesignIdDigest(String(root, "design_id", "$"));

    // --- planform ---------------------------------------------------------
    {
        const auto& planform = Object(root, "planform", "$");
        contract.Span = Number(planform, "span", "planform");
        RequirePositive(contract.Span, "planform.span");

        const auto& stations = Array(planform, "stations", "planform");
        if (stations.size() < MinPlanformStations)
            throw ContractError(std::format("planform.stations: expected at least {} stations, got {}",
                                            MinPlanformStations, stations.size()));

        std::vector<double> etas;
        etas.reserve(stations.size());
        for (std::size_t i = 0; i < stations.size(); ++i) {
            const auto where = std::format("planform.stations[{}]", i);
            PlanformStation station;
            station.Eta = Number(stations[i], "eta", where);
            station.Chord = Number(stations[i], "chord", where);
            station.TwistDeg = Number(stations[i], "twist", where);
            station.SweepQuarterChordDeg = Number(stations[i], "sweep_qc", where);
            station.DihedralDeg = Number(stations[i], "dihedral", where);
            RequirePositive(station.Chord, std::format("{}.chord", where));
            etas.push_back(station.Eta);
            contract.Stations.push_back(station);
        }
        RequireSpanningEtaSequence(etas, "planform.stations");

        // Optional: absent before schema 1.5.0. Without it the wing and the
        // body have no stated relationship at all (see WingPlacement.h).
        if (planform.contains("placement")) {
            const auto& placement = Object(planform, "placement", "planform");
            const auto& anchor = Object(placement, "root_leading_edge", "planform.placement");
            contract.Placement.IsStated = true;
            contract.Placement.RootLeadingEdge =
                Math::Vec3(Number(anchor, "x", "planform.placement.root_leading_edge"),
                           Number(anchor, "y", "planform.placement.root_leading_edge"),
                           Number(anchor, "z", "planform.placement.root_leading_edge"));
        }
    }

    // --- airfoil sections -------------------------------------------------
    {
        const auto& sections = Array(root, "airfoil_sections", "$");
        if (sections.empty()) throw ContractError("airfoil_sections: must not be empty");

        std::vector<double> etas;
        etas.reserve(sections.size());
        for (std::size_t i = 0; i < sections.size(); ++i) {
            const auto where = std::format("airfoil_sections[{}]", i);
            const auto parameterization = String(sections[i], "parameterization", where);
            if (parameterization != "CST")
                throw ContractError(std::format("{}.parameterization: expected 'CST', got '{}'", where,
                                                parameterization));

            AirfoilSection section;
            section.Eta = Number(sections[i], "eta", where);
            section.CoefficientsUpper = Coefficients(sections[i], "coefficients_upper", where);
            section.CoefficientsLower = Coefficients(sections[i], "coefficients_lower", where);
            if (section.CoefficientsUpper.size() != section.CoefficientsLower.size())
                throw ContractError(std::format(
                    "{}: upper and lower surfaces must share a CST order, got {} and {} coefficients", where,
                    section.CoefficientsUpper.size(), section.CoefficientsLower.size()));
            etas.push_back(section.Eta);
            contract.AirfoilSections.push_back(std::move(section));
        }
        RequireSpanningEtaSequence(etas, "airfoil_sections");
    }

    // --- control surfaces -------------------------------------------------
    // Optional: a clean wing has none. Names are not unique (see
    // ControlSurface.h), so nothing here dedupes them.
    if (root.contains("control_surfaces")) {
        const auto& surfaces = Array(root, "control_surfaces", "$");
        for (std::size_t i = 0; i < surfaces.size(); ++i) {
            const auto where = std::format("control_surfaces[{}]", i);
            ControlSurface surface;
            surface.Name = String(surfaces[i], "name", where);
            if (surface.Name.empty()) throw ContractError(std::format("{}.name: must not be empty", where));
            surface.ChordFraction = Number(surfaces[i], "chord_fraction", where);
            surface.EtaStart = Number(surfaces[i], "eta_start", where);
            surface.EtaEnd = Number(surfaces[i], "eta_end", where);
            surface.HingeAxis = HingeAxis(surfaces[i], where);
            surface.Limits = ParseDeflectionLimits(surfaces[i], where);
            // Forward-compatible with schema 1.1.0's explicit 'surface' field:
            // honour it when present, fall back to the 1.0.0 name convention.
            surface.Binding = surfaces[i].contains("surface")
                                  ? ParseBindingField(String(surfaces[i], "surface", where), where)
                                  : BindingFromName(surface.Name, where);

            // chord_fraction == 1 is legitimate: an all-moving vane.
            RequireInRange(surface.ChordFraction, 0.0, 1.0, std::format("{}.chord_fraction", where));
            RequirePositive(surface.ChordFraction, std::format("{}.chord_fraction", where));
            RequireInRange(surface.EtaStart, EtaMin, EtaMax, std::format("{}.eta_start", where));
            RequireInRange(surface.EtaEnd, EtaMin, EtaMax, std::format("{}.eta_end", where));
            if (!(surface.EtaEnd > surface.EtaStart))
                throw ContractError(std::format("{}: eta_end ({}) must exceed eta_start ({})", where,
                                                surface.EtaEnd, surface.EtaStart));
            contract.ControlSurfaces.push_back(std::move(surface));
        }
    }

    // --- mesh topology ----------------------------------------------------
    {
        const auto& mesh = Object(root, "mesh_topology", "$");
        contract.Mesh.ChordwisePanels = Integer(mesh, "chordwise_panels", "mesh_topology");
        contract.Mesh.SpanwisePanelsPerSection = Integer(mesh, "spanwise_panels_per_section", "mesh_topology");
        RequireSupportedWakeModel(String(mesh, "wake_model", "mesh_topology"), "mesh_topology");
        if (contract.Mesh.ChordwisePanels < 1)
            throw ContractError(std::format("mesh_topology.chordwise_panels must be >= 1, got {}",
                                            contract.Mesh.ChordwisePanels));
        if (contract.Mesh.SpanwisePanelsPerSection < 1)
            throw ContractError(std::format("mesh_topology.spanwise_panels_per_section must be >= 1, got {}",
                                            contract.Mesh.SpanwisePanelsPerSection));
    }

    // --- body -------------------------------------------------------------
    // Optional: absent before schema 1.4.0, and a flying wing has none.
    if (root.contains("body")) contract.Body = ParseBody(Object(root, "body", "$"));

    // --- duct ---------------------------------------------------------
    // Optional: absent before schema 1.8.0, and an unducted propeller has
    // none. Disjoint from `body` -- see DuctGeometry.h.
    if (root.contains("duct")) contract.Duct = ParseDuct(Object(root, "duct", "$"));

    // --- moment reference point ---------------------------------------
    // Optional: absent before schema 1.8.0.
    if (root.contains("moment_reference_point")) {
        const auto& mrp = Object(root, "moment_reference_point", "$");
        contract.MomentReferencePointStated = true;
        contract.MomentReferencePoint = Math::Vec3(Number(mrp, "x", "moment_reference_point"),
                                                    Number(mrp, "y", "moment_reference_point"),
                                                    Number(mrp, "z", "moment_reference_point"));
    }

    // --- propulsion -------------------------------------------------------
    // Optional: an unpowered glider has no propeller block.
    if (root.contains("propulsion_bemt")) {
        const auto& propulsion = Object(root, "propulsion_bemt", "$");
        contract.Propulsion.DiskRadius = Number(propulsion, "disk_radius", "propulsion_bemt");
        contract.Propulsion.ReferenceRpm = Number(propulsion, "reference_rpm", "propulsion_bemt");
        RequirePositive(contract.Propulsion.DiskRadius, "propulsion_bemt.disk_radius");
        RequirePositive(contract.Propulsion.ReferenceRpm, "propulsion_bemt.reference_rpm");

        // Blade count arrived in 1.4.0; older documents leave it zero for the
        // consumer to supply, as they always had to.
        if (propulsion.contains("n_blades")) {
            contract.Propulsion.BladeCount = Integer(propulsion, "n_blades", "propulsion_bemt");
            if (contract.Propulsion.BladeCount < 2)
                throw ContractError(std::format("propulsion_bemt.n_blades must be >= 2, got {}",
                                                contract.Propulsion.BladeCount));
        }

        const auto& stations = Array(propulsion, "blade_stations", "propulsion_bemt");
        if (stations.size() < MinBladeStations)
            throw ContractError(std::format("propulsion_bemt.blade_stations: expected at least {} stations, got {}",
                                            MinBladeStations, stations.size()));

        for (std::size_t i = 0; i < stations.size(); ++i) {
            const auto where = std::format("propulsion_bemt.blade_stations[{}]", i);
            BladeStationSpec station;
            station.RadiusFraction = Number(stations[i], "r_over_R", where);
            station.Chord = Number(stations[i], "chord", where);
            station.TwistDeg = Number(stations[i], "twist", where);

            RequireInRange(station.RadiusFraction, 0.0, 1.0, std::format("{}.r_over_R", where));
            // Chord tapers to exactly zero at a rounded tip, so this is >= 0,
            // not > 0 -- unlike the planform chords above.
            if (station.Chord < 0.0)
                throw ContractError(std::format("{}.chord must not be negative, got {}", where, station.Chord));
            if (i > 0 && !(station.RadiusFraction > contract.Propulsion.BladeStations.back().RadiusFraction))
                throw ContractError(std::format("{}: r_over_R must strictly increase, but {} follows {}", where,
                                                station.RadiusFraction,
                                                contract.Propulsion.BladeStations.back().RadiusFraction));
            contract.Propulsion.BladeStations.push_back(station);
        }

        // RotationAxis arrived in schema 1.8.0; older documents leave it the
        // zero vector, same as the rest of this block's 1.8.0 additions.
        if (propulsion.contains("rotation_axis")) {
            const Math::Vec3 axis = Vec3Array(propulsion, "rotation_axis", "propulsion_bemt");
            if (axis.Norm() < MinAxisNorm)
                throw ContractError("propulsion_bemt.rotation_axis: must not be the zero vector");
            contract.Propulsion.RotationAxis = axis.Normalized();
        }

        // Blade airfoil shape arrived in schema 1.8.0; older documents state
        // only chord/twist per station, with no section shape at all. Sections
        // are keyed by r_over_R rather than eta -- radial position along the
        // blade, not semi-span position along the wing -- so they are not
        // parsed as AirfoilSection/RequireSpanningEtaSequence: unlike the wing,
        // a blade need not have a section at r/R = 0 (the hub cutout starts it
        // higher) or r/R = 1.
        if (propulsion.contains("airfoil_sections")) {
            const auto& sections = Array(propulsion, "airfoil_sections", "propulsion_bemt");
            for (std::size_t i = 0; i < sections.size(); ++i) {
                const auto where = std::format("propulsion_bemt.airfoil_sections[{}]", i);
                const auto parameterization = String(sections[i], "parameterization", where);
                if (parameterization != "CST")
                    throw ContractError(std::format("{}.parameterization: expected 'CST', got '{}'", where,
                                                    parameterization));

                BladeAirfoilSection section;
                section.RadiusFraction = Number(sections[i], "r_over_R", where);
                RequireInRange(section.RadiusFraction, 0.0, 1.0, std::format("{}.r_over_R", where));
                section.CoefficientsUpper = Coefficients(sections[i], "coefficients_upper", where);
                section.CoefficientsLower = Coefficients(sections[i], "coefficients_lower", where);
                if (section.CoefficientsUpper.size() != section.CoefficientsLower.size())
                    throw ContractError(std::format(
                        "{}: upper and lower surfaces must share a CST order, got {} and {} coefficients", where,
                        section.CoefficientsUpper.size(), section.CoefficientsLower.size()));
                if (i > 0 &&
                    !(section.RadiusFraction > contract.Propulsion.BladeAirfoilSections.back().RadiusFraction))
                    throw ContractError(std::format(
                        "{}: r_over_R must strictly increase, but {} follows {}", where, section.RadiusFraction,
                        contract.Propulsion.BladeAirfoilSections.back().RadiusFraction));
                contract.Propulsion.BladeAirfoilSections.push_back(std::move(section));
            }
        }
    }

    return contract;
}

// --- bridge to the parametric single-trapezoid solver wing -----------------
// BuildWing() models ONE linearly-tapered, linearly-twisted trapezoid, which
// is strictly less expressive than the contract's per-station planform. This
// reduction therefore refuses anything it cannot represent exactly rather than
// quietly fitting a trapezoid through a shape that is not one -- a silently
// approximated planform is the kind of error that only shows up as a few
// percent of misattributed CL much later.
//
// A station-resolved lattice builder would supersede this entirely; until then
// this is the honest narrow path.
inline constexpr double TrapezoidFitTolerance = 1e-9;

/**
 * Reduce a per-station planform to the parametric solver's single trapezoid.
 * Throws ContractError if the planform is not exactly one linearly-tapered,
 * linearly-twisted trapezoid (sweep/dihedral constant, chord/twist linear in
 * eta, zero root incidence).
 */
[[nodiscard]] inline Solver::WingParams ToWingParams(const HandoffContract& contract) {
    const auto& stations = contract.Stations;
    const PlanformStation& root = stations.front();
    const PlanformStation& tip = stations.back();

    for (std::size_t i = 0; i < stations.size(); ++i) {
        const auto& station = stations[i];
        if (std::fabs(station.SweepQuarterChordDeg - root.SweepQuarterChordDeg) > TrapezoidFitTolerance ||
            std::fabs(station.DihedralDeg - root.DihedralDeg) > TrapezoidFitTolerance)
            throw ContractError(std::format(
                "planform.stations[{}]: sweep/dihedral vary along the span; not a single trapezoid", i));

        // Chord and twist must be exactly the linear root->tip interpolants.
        const double chord = root.Chord + (tip.Chord - root.Chord) * station.Eta;
        const double twist = root.TwistDeg + (tip.TwistDeg - root.TwistDeg) * station.Eta;
        if (std::fabs(station.Chord - chord) > TrapezoidFitTolerance * root.Chord)
            throw ContractError(std::format(
                "planform.stations[{}]: chord {} breaks linear taper (expected {})", i, station.Chord, chord));
        if (std::fabs(station.TwistDeg - twist) > TrapezoidFitTolerance)
            throw ContractError(std::format(
                "planform.stations[{}]: twist {} is not linear in eta (expected {})", i, station.TwistDeg, twist));
    }
    // BuildWing() hard-codes zero twist at the root and measures TwistTipDeg
    // relative to it, so a rigged root incidence cannot be represented.
    if (std::fabs(root.TwistDeg) > TrapezoidFitTolerance)
        throw ContractError(std::format("planform.stations[0].twist: root incidence {} is not representable",
                                        root.TwistDeg));

    Solver::WingParams wing;
    wing.Span = contract.Span;
    wing.RootChord = root.Chord;
    wing.TipChord = tip.Chord;
    wing.SweepQuarterChordDeg = root.SweepQuarterChordDeg;
    wing.DihedralDeg = root.DihedralDeg;
    wing.TwistTipDeg = tip.TwistDeg;
    wing.NPanelsSemiSpan =
        contract.Mesh.SpanwisePanelsPerSection * static_cast<int>(stations.size() - 1);
    // Spacing is not part of the schema; the lattice keeps its uniform default.
    wing.CosineSpacing = false;
    return wing;
}

// --- contract-to-metric propeller conversion ---------------------------
// PropulsionSpec states radii as fractions of the disk radius, because
// that is the form the blade is designed in; Geometry::Propeller wants
// metres, because that is the form every consumer works in. This is that
// conversion, and it is the only place it should happen. (The BEMT solver
// that used to sit behind this bridge is now a separate project this repo
// does not depend on; the metric Propeller struct is what any in-repo
// consumer -- the viewer today, a calculation method tomorrow -- poses
// itself on.)
//
// The hub radius is the innermost blade station: the contract defines the
// blade only outboard of it (see PropulsionSpec.h), so that station IS the
// cutout rather than a point on a blade that continues inward.
//
// Blade count came into the schema at 1.4.0. Older contracts leave it zero,
// so the caller must supply one; that is refused loudly here rather than
// silently defaulting, because a wrong blade count changes everything a
// consumer derives from the propeller and would still look plausible.
inline constexpr int MinBladeCount = 2;

/**
 * Convert a contract's fraction-of-disk-radius blade stations to the metric
 * Geometry::Propeller every consumer works with. bladeCountOverride is
 * required for contracts older than schema 1.4.0, which carry no blade
 * count of their own.
 */
[[nodiscard]] inline Propeller ToPropeller(const PropulsionSpec& spec, int bladeCountOverride = 0) {
    if (spec.BladeStations.size() < MinBladeStations)
        throw ContractError("propulsion_bemt: no blade stations to build a propeller from");
    Detail::RequirePositive(spec.DiskRadius, "propulsion_bemt.disk_radius");

    const int blades = (bladeCountOverride > 0) ? bladeCountOverride : spec.BladeCount;
    if (blades < MinBladeCount)
        throw ContractError(std::format(
            "propulsion_bemt: blade count {} is unusable; this contract states none (schema < 1.4.0), "
            "so the caller must supply one",
            blades));

    Propeller prop;
    prop.BladeCount = blades;
    prop.Radius = spec.DiskRadius;
    prop.HubRadius = spec.BladeStations.front().RadiusFraction * spec.DiskRadius;
    prop.Stations.reserve(spec.BladeStations.size());
    for (const BladeStationSpec& station : spec.BladeStations)
        prop.Stations.push_back({station.RadiusFraction * spec.DiskRadius, station.Chord, station.TwistDeg});
    return prop;
}

/** Read and parse a handoff document from disk. */
[[nodiscard]] inline HandoffContract LoadHandoff(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open geometry handoff: " + path);
    nlohmann::json root;
    try {
        input >> root;
    } catch (const nlohmann::json::parse_error& error) {
        throw ContractError(std::format("{}: malformed JSON: {}", path, error.what()));
    }
    return ParseHandoff(root);
}

} // namespace Aeolion::Geometry
