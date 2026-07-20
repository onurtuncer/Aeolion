#include "Aeolion/Geometry/HandoffContract.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using Aeolion::Geometry::ContractError;
using Aeolion::Geometry::ControlSurfaceBinding;

namespace {

int g_Failures = 0;

void Check(bool condition, const std::string& what) {
    if (!condition) {
        std::cerr << "FAIL: " << what << "\n";
        ++g_Failures;
    }
}

void CheckClose(double actual, double expected, const std::string& what, double tolerance = 1e-12) {
    if (!(std::fabs(actual - expected) <= tolerance)) {
        std::cerr << "FAIL: " << what << " -- expected " << expected << ", got " << actual << "\n";
        ++g_Failures;
    }
}

std::string FixturePath(const std::string& version) {
    return std::string(AEOLION_TEST_DATA_DIR) + "/AeolionGeometryHandoff-" + version + ".json";
}

// The 1.0.0 document is the baseline for the invariant tests; 1.1.0 is the
// same geometry re-emitted with explicit 'surface' fields.
const std::string& FixturePath() {
    static const std::string path = FixturePath("1.0.0");
    return path;
}

nlohmann::json JsonAt(const std::string& path) {
    std::ifstream input(path);
    nlohmann::json root;
    input >> root;
    return root;
}

nlohmann::json FixtureJson() { return JsonAt(FixturePath()); }

// Applies `mutate` to a fresh copy of the fixture and asserts the parser
// rejects the result. A contract parser that accepts malformed input is worse
// than no parser at all, so every invariant gets a negative test.
template <typename Mutate>
void CheckRejects(const std::string& what, Mutate mutate) {
    nlohmann::json root = FixtureJson();
    mutate(root);
    try {
        (void)Aeolion::Geometry::ParseHandoff(root);
    } catch (const ContractError&) {
        return;
    }
    std::cerr << "FAIL: expected rejection -- " << what << "\n";
    ++g_Failures;
}

void TestParsesFixture() {
    const auto contract = Aeolion::Geometry::LoadHandoff(FixturePath());

    Check(contract.SchemaVersion == "1.0.0", "schema_version round-trips");
    Check(contract.ReferenceFrame == "aetherion_body_frd", "reference_frame round-trips");
    Check(contract.DesignId == "2732d99b6276483926f2a4047cf3dc0be31699a3aad09c93d3dddf516ec24a63",
          "design_id digest has its 'sha256:' prefix stripped");

    CheckClose(contract.Span, 1.0628798195068672, "planform.span");
    Check(contract.Stations.size() == 5, "5 planform stations");
    CheckClose(contract.Stations.front().Eta, 0.0, "root station eta");
    CheckClose(contract.Stations.back().Eta, 1.0, "tip station eta");
    CheckClose(contract.Stations[2].Chord, 0.17714663658447782, "mid station chord");
    CheckClose(contract.Stations[2].TwistDeg, 0.0, "mid station twist");

    Check(contract.AirfoilSections.size() == 5, "5 airfoil sections");
    for (const auto& section : contract.AirfoilSections) {
        Check(section.CoefficientsUpper.size() == 6, "6 upper CST coefficients");
        Check(section.CoefficientsLower.size() == 6, "6 lower CST coefficients");
    }
    // A_0 governs the leading-edge radius and A_n the trailing-edge boat-tail
    // angle, so the endpoints are the ones worth pinning exactly.
    CheckClose(contract.AirfoilSections[0].CoefficientsUpper.front(), 0.21809304448357872, "upper A_0");
    CheckClose(contract.AirfoilSections[0].CoefficientsUpper.back(), 0.269483355187719, "upper A_n");
    CheckClose(contract.AirfoilSections[0].CoefficientsLower.front(), -0.13411949093207398, "lower A_0");
    CheckClose(contract.AirfoilSections[0].CoefficientsLower.back(), -0.011931451366098689, "lower A_n");

    Check(contract.ControlSurfaces.size() == 5, "5 control surfaces");
    Check(contract.ControlSurfaces[0].Name == "aileron", "first control surface is the aileron");
    CheckClose(contract.ControlSurfaces[0].ChordFraction, 0.12, "aileron chord fraction");
    CheckClose(contract.ControlSurfaces[0].EtaStart, 0.88, "aileron eta_start");
    CheckClose(contract.ControlSurfaces[0].EtaEnd, 1.0, "aileron eta_end");
    // Four entries share the name "vane"; they must all survive as distinct
    // entries rather than collapsing into one.
    int vanes = 0;
    for (const auto& surface : contract.ControlSurfaces)
        if (surface.Name == "vane") ++vanes;
    Check(vanes == 4, "all four identically-named vanes are kept");
    CheckClose(contract.ControlSurfaces[1].ChordFraction, 1.0, "an all-moving vane may span the full chord");
    CheckClose(contract.ControlSurfaces[2].HingeAxis.y, -1.0, "hinge axis sign is preserved");

    Check(contract.Mesh.ChordwisePanels == 8, "chordwise panel count");
    Check(contract.Mesh.SpanwisePanelsPerSection == 6, "spanwise panel count per section");

    Check(contract.Propulsion.BladeStations.size() == 13, "13 blade stations");
    CheckClose(contract.Propulsion.DiskRadius, 0.1015, "disk radius");
    CheckClose(contract.Propulsion.ReferenceRpm, 12211.259999999998, "reference rpm");
    CheckClose(contract.Propulsion.BladeStations.front().RadiusFraction, 0.42, "hub cutout r/R");
    CheckClose(contract.Propulsion.BladeStations.back().RadiusFraction, 1.0, "tip r/R");
    CheckClose(contract.Propulsion.BladeStations.back().Chord, 0.0, "tip chord tapers to zero");
    CheckClose(contract.Propulsion.BladeStations.back().TwistDeg, 13.427041758485208, "tip twist");
}

// The 1.0.0 schema carries no surface association, so the parser has to
// recover it from the control's name. Getting this wrong would put all-moving
// duct-jet vanes into the wing lattice -- with eta bands that are radius
// fractions, not semi-span fractions.
void TestSurfaceBinding() {
    const auto contract = Aeolion::Geometry::LoadHandoff(FixturePath());

    Check(contract.ControlSurfaces[0].Binding == ControlSurfaceBinding::Wing, "the aileron binds to the wing");
    for (std::size_t i = 1; i < contract.ControlSurfaces.size(); ++i)
        Check(contract.ControlSurfaces[i].Binding == ControlSurfaceBinding::DuctJet,
              "every vane binds to the duct jet, not the wing");

    // An explicit 'surface' field (schema >= 1.1.0) overrides the name rule,
    // and is the escape hatch for names the table does not know.
    nlohmann::json root = FixtureJson();
    root["control_surfaces"][1]["surface"] = "wing";
    const auto overridden = Aeolion::Geometry::ParseHandoff(root);
    Check(overridden.ControlSurfaces[1].Binding == ControlSurfaceBinding::Wing,
          "an explicit 'surface' field overrides the name convention");

    root = FixtureJson();
    root["control_surfaces"][0]["name"] = "elevon";
    root["control_surfaces"][0]["surface"] = "wing";
    try {
        const auto named = Aeolion::Geometry::ParseHandoff(root);
        Check(named.ControlSurfaces[0].Binding == ControlSurfaceBinding::Wing,
              "an unknown name is accepted when 'surface' states the binding");
    } catch (const ContractError& error) {
        std::cerr << "FAIL: 'surface' should admit unknown names -- " << error.what() << "\n";
        ++g_Failures;
    }
}

// Schema 1.1.0 states each control's binding explicitly instead of leaving it
// implicit in the name. The two must never disagree: the name table lives in
// this repo and the explicit field is written by the exporter, so if they ever
// drift, one of the two repos is silently mis-binding controls. Parsing the
// 1.1.0 document with and without its 'surface' fields must therefore give
// identical bindings.
void TestSchema110AgreesWithNameConvention() {
    const auto explicitJson = JsonAt(FixturePath("1.1.0"));

    nlohmann::json impliedJson = explicitJson;
    bool sawSurfaceField = false;
    for (auto& surface : impliedJson["control_surfaces"]) {
        if (surface.contains("surface")) sawSurfaceField = true;
        surface.erase("surface");
    }
    Check(sawSurfaceField, "the 1.1.0 fixture actually carries 'surface' fields");

    const auto fromField = Aeolion::Geometry::ParseHandoff(explicitJson);
    const auto fromName = Aeolion::Geometry::ParseHandoff(impliedJson);

    Check(fromField.SchemaVersion == "1.1.0", "the 1.1.0 document parses");
    Check(fromField.ControlSurfaces.size() == fromName.ControlSurfaces.size(), "control surface counts agree");
    for (std::size_t i = 0; i < fromField.ControlSurfaces.size(); ++i)
        Check(fromField.ControlSurfaces[i].Binding == fromName.ControlSurfaces[i].Binding,
              "explicit 'surface' agrees with the name convention for control_surfaces[" + std::to_string(i) + "]");

    // 1.1.0 added metadata only. The geometry -- and so the aerodynamics --
    // must be untouched relative to 1.0.0, even though design_id moved.
    const auto baseline = Aeolion::Geometry::LoadHandoff(FixturePath("1.0.0"));
    CheckClose(fromField.Span, baseline.Span, "span is unchanged across the version bump");
    Check(fromField.DesignId != baseline.DesignId, "design_id moves with the schema bump");
    const auto wing = Aeolion::Geometry::ToWingParams(fromField);
    const auto baselineWing = Aeolion::Geometry::ToWingParams(baseline);
    CheckClose(wing.RootChord, baselineWing.RootChord, "root chord is unchanged");
    CheckClose(wing.TipChord, baselineWing.TipChord, "tip chord is unchanged");
    Check(wing.NPanelsSemiSpan == baselineWing.NPanelsSemiSpan, "lattice density is unchanged");
}

// Reducing the per-station planform to BuildWing()'s single trapezoid is lossy,
// so it must refuse everything it cannot represent exactly.
void TestToWingParams() {
    const auto contract = Aeolion::Geometry::LoadHandoff(FixturePath());
    const auto wing = Aeolion::Geometry::ToWingParams(contract);

    CheckClose(wing.Span, 1.0628798195068672, "wing span");
    CheckClose(wing.RootChord, 0.17714663658447782, "wing root chord");
    CheckClose(wing.TipChord, 0.17714663658447782, "wing tip chord");
    CheckClose(wing.SweepQuarterChordDeg, 0.0, "wing sweep");
    CheckClose(wing.DihedralDeg, 0.0, "wing dihedral");
    CheckClose(wing.TwistTipDeg, 0.0, "wing tip twist");
    // 4 station intervals x 6 panels per section.
    Check(wing.NPanelsSemiSpan == 24, "panels per semi-span come from the mesh topology");

    const auto rejects = [](const std::string& what, auto mutate) {
        nlohmann::json root = FixtureJson();
        mutate(root);
        try {
            (void)Aeolion::Geometry::ToWingParams(Aeolion::Geometry::ParseHandoff(root));
        } catch (const ContractError&) {
            return;
        }
        std::cerr << "FAIL: expected rejection -- " << what << "\n";
        ++g_Failures;
    };

    rejects("a planform with a kinked chord distribution",
            [](auto& root) { root["planform"]["stations"][2]["chord"] = 0.3; });
    rejects("a planform with a spanwise sweep break",
            [](auto& root) { root["planform"]["stations"][3]["sweep_qc"] = 10.0; });
    rejects("a planform with a spanwise dihedral break",
            [](auto& root) { root["planform"]["stations"][3]["dihedral"] = 5.0; });
    rejects("a non-linear twist distribution",
            [](auto& root) { root["planform"]["stations"][2]["twist"] = -3.0; });
    rejects("a rigged root incidence", [](auto& root) {
        // Uniform incidence is linear in eta, but BuildWing() pins the root to
        // zero, so it still is not representable.
        for (auto& station : root["planform"]["stations"]) station["twist"] = 2.0;
    });
}

void TestNormalizesHingeAxis() {
    nlohmann::json root = FixtureJson();
    root["control_surfaces"][0]["hinge_axis"] = {0.0, 2.0, 0.0};
    const auto contract = Aeolion::Geometry::ParseHandoff(root);
    CheckClose(contract.ControlSurfaces[0].HingeAxis.y, 1.0, "hinge axis is normalized to unit length");
}

void TestToleratesUnknownKeys() {
    nlohmann::json root = FixtureJson();
    root["future_block"] = {{"anything", 1}};
    try {
        (void)Aeolion::Geometry::ParseHandoff(root);
    } catch (const ContractError& error) {
        std::cerr << "FAIL: unknown top-level keys should be tolerated -- " << error.what() << "\n";
        ++g_Failures;
    }
}

void TestOptionalBlocks() {
    nlohmann::json root = FixtureJson();
    root.erase("control_surfaces");
    root.erase("propulsion_bemt");
    try {
        const auto contract = Aeolion::Geometry::ParseHandoff(root);
        Check(contract.ControlSurfaces.empty(), "a contract with no control surfaces parses");
        Check(contract.Propulsion.BladeStations.empty(), "a contract with no propeller parses");
    } catch (const ContractError& error) {
        std::cerr << "FAIL: control_surfaces and propulsion_bemt should be optional -- " << error.what() << "\n";
        ++g_Failures;
    }
}

// --- schema 1.4.0: body, blade count, deflection limits --------------------
// The 1.4.0 document adds three blocks that older ones lack. Because the
// parser tolerates unknown keys for forward compatibility, an unwired block
// parses SILENTLY -- so these checks are what distinguish "understood" from
// "ignored".
void TestSchema140Blocks() {
    const auto contract = Aeolion::Geometry::LoadHandoff(FixturePath("1.4.0"));

    Check(contract.Body.IsPresent(), "1.4.0 carries a body");
    CheckClose(contract.Body.Length, 0.49099, "body length");
    Check(contract.Body.Stations.size() == 25, "25 body stations");
    CheckClose(contract.Body.Stations.front().x, 0.0, "nose at x=0");
    CheckClose(contract.Body.Stations.front().Radius, 0.0, "sharp nose");
    CheckClose(contract.Body.Stations.back().x, -0.49099, "tail x");
    // The tail does NOT close: this body is truncated at the duct exit, and
    // anything that panels it has to cope with an open base.
    CheckClose(contract.Body.Stations.back().Radius, 0.0406, "open base radius");

    Check(contract.Propulsion.BladeCount == 3, "blade count comes from n_blades");

    const auto& aileron = contract.ControlSurfaces.front();
    Check(aileron.Binding == ControlSurfaceBinding::Wing, "aileron binds to the wing");
    CheckClose(aileron.Limits.MinDeg, -20.0, "aileron min deflection");
    CheckClose(aileron.Limits.MaxDeg, 20.0, "aileron max deflection");
    Check(!aileron.Limits.HasSoftLimit, "the aileron states no soft limit");

    const auto& vane = contract.ControlSurfaces[1];
    Check(vane.Binding == ControlSurfaceBinding::DuctJet, "vane binds to the duct jet");
    Check(vane.Limits.HasSoftLimit, "the vane states a soft limit");
    CheckClose(vane.Limits.SoftLimitDeg, 15.0, "vane soft limit");

    // Older documents must keep parsing, with the new fields simply absent.
    const auto older = Aeolion::Geometry::LoadHandoff(FixturePath("1.0.0"));
    Check(!older.Body.IsPresent(), "1.0.0 has no body");
    Check(older.Propulsion.BladeCount == 0, "1.0.0 states no blade count");
    Check(!older.ControlSurfaces.front().Limits.HasSoftLimit, "1.0.0 states no soft limit");
}

// Every invariant the 1.4.0 blocks rely on gets a negative test, same as the
// rest of the contract.
void TestRejectsMalformed140() {
    const auto rejects = [](const std::string& what, auto mutate) {
        nlohmann::json root = JsonAt(FixturePath("1.4.0"));
        mutate(root);
        try {
            (void)Aeolion::Geometry::ParseHandoff(root);
            std::cerr << "FAIL: expected rejection -- " << what << "\n";
            ++g_Failures;
        } catch (const ContractError&) {
        }
    };

    rejects("a body whose x does not decrease nose to tail",
            [](auto& root) { root["body"]["stations"][3]["x"] = 0.5; });
    rejects("a negative body radius", [](auto& root) { root["body"]["stations"][3]["radius"] = -0.01; });
    rejects("a body length disagreeing with its stations",
            [](auto& root) { root["body"]["length"] = 0.6; });
    rejects("a body with one station", [](auto& root) {
        root["body"]["stations"] = nlohmann::json::array({root["body"]["stations"][0]});
    });
    rejects("a single-bladed propeller", [](auto& root) { root["propulsion_bemt"]["n_blades"] = 1; });
    rejects("deflection limits whose max does not exceed min",
            [](auto& root) { root["control_surfaces"][0]["deflection_limits_deg"]["max"] = -30.0; });
    rejects("a soft limit outside the hard stops",
            [](auto& root) { root["control_surfaces"][1]["deflection_soft_limit_deg"] = 25.0; });
}

void TestRejectsMalformed() {
    CheckRejects("a future schema major version",
                 [](auto& root) { root["schema_version"] = "2.0.0"; });
    CheckRejects("a non-semver schema version",
                 [](auto& root) { root["schema_version"] = "1.0"; });
    CheckRejects("a numeric schema version",
                 [](auto& root) { root["schema_version"] = 1; });
    CheckRejects("units in feet",
                 [](auto& root) { root["units"]["length"] = "ft"; });
    CheckRejects("angles in radians",
                 [](auto& root) { root["units"]["angle"] = "rad"; });
    CheckRejects("an unexpected reference frame",
                 [](auto& root) { root["reference_frame"] = "ned"; });
    CheckRejects("a design_id without its algorithm prefix",
                 [](auto& root) { root["design_id"] = "2732d99b"; });
    CheckRejects("a truncated design_id digest",
                 [](auto& root) { root["design_id"] = "sha256:2732d99b"; });
    CheckRejects("a non-hex design_id digest", [](auto& root) {
        root["design_id"] = "sha256:zzzzd99b6276483926f2a4047cf3dc0be31699a3aad09c93d3dddf516ec24a63";
    });

    CheckRejects("a missing planform", [](auto& root) { root.erase("planform"); });
    CheckRejects("a non-positive span", [](auto& root) { root["planform"]["span"] = 0.0; });
    CheckRejects("a single planform station", [](auto& root) {
        root["planform"]["stations"] = nlohmann::json::array({root["planform"]["stations"][0]});
    });
    CheckRejects("a zero chord", [](auto& root) { root["planform"]["stations"][1]["chord"] = 0.0; });
    CheckRejects("a missing chord", [](auto& root) { root["planform"]["stations"][1].erase("chord"); });
    CheckRejects("a chord given as a string", [](auto& root) { root["planform"]["stations"][1]["chord"] = "0.2"; });
    CheckRejects("out-of-order station eta", [](auto& root) { root["planform"]["stations"][1]["eta"] = 0.9; });
    CheckRejects("duplicated station eta", [](auto& root) { root["planform"]["stations"][1]["eta"] = 0.0; });
    CheckRejects("stations that do not start at the root",
                 [](auto& root) { root["planform"]["stations"][0]["eta"] = 0.1; });
    CheckRejects("stations that do not reach the tip",
                 [](auto& root) { root["planform"]["stations"][4]["eta"] = 0.9; });

    CheckRejects("no airfoil sections",
                 [](auto& root) { root["airfoil_sections"] = nlohmann::json::array(); });
    CheckRejects("a non-CST parameterization",
                 [](auto& root) { root["airfoil_sections"][0]["parameterization"] = "bspline"; });
    // 'relaxed' is a legal schema value that this solver cannot honour, so
    // it must be refused rather than silently solved with a frozen wake.
    CheckRejects("a relaxed wake, which is not implemented",
                 [](auto& root) { root["mesh_topology"]["wake_model"] = "relaxed"; });
    CheckRejects("an unknown wake model",
                 [](auto& root) { root["mesh_topology"]["wake_model"] = "banana"; });
    CheckRejects("too few CST coefficients",
                 [](auto& root) { root["airfoil_sections"][0]["coefficients_upper"] = {0.1, 0.2, 0.3}; });
    CheckRejects("too many CST coefficients", [](auto& root) {
        root["airfoil_sections"][0]["coefficients_upper"] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    });
    CheckRejects("mismatched upper/lower CST order",
                 [](auto& root) { root["airfoil_sections"][0]["coefficients_lower"] = {0.1, 0.2, 0.3, 0.4}; });
    CheckRejects("a non-numeric CST coefficient",
                 [](auto& root) { root["airfoil_sections"][0]["coefficients_upper"][2] = nullptr; });
    CheckRejects("out-of-order section eta",
                 [](auto& root) { root["airfoil_sections"][1]["eta"] = 0.9; });

    CheckRejects("an unnamed control surface",
                 [](auto& root) { root["control_surfaces"][0]["name"] = ""; });
    CheckRejects("a zero chord fraction",
                 [](auto& root) { root["control_surfaces"][0]["chord_fraction"] = 0.0; });
    CheckRejects("a chord fraction beyond the full chord",
                 [](auto& root) { root["control_surfaces"][0]["chord_fraction"] = 1.5; });
    CheckRejects("an inverted control surface span",
                 [](auto& root) { root["control_surfaces"][0]["eta_start"] = 1.0; });
    CheckRejects("a zero-width control surface",
                 [](auto& root) { root["control_surfaces"][0]["eta_end"] = 0.88; });
    CheckRejects("a control surface whose name implies no binding",
                 [](auto& root) { root["control_surfaces"][0]["name"] = "spoiler"; });
    CheckRejects("an unknown explicit surface binding",
                 [](auto& root) { root["control_surfaces"][0]["surface"] = "fuselage"; });
    CheckRejects("a zero hinge axis",
                 [](auto& root) { root["control_surfaces"][0]["hinge_axis"] = {0.0, 0.0, 0.0}; });
    CheckRejects("a two-component hinge axis",
                 [](auto& root) { root["control_surfaces"][0]["hinge_axis"] = {0.0, 1.0}; });

    CheckRejects("a missing mesh topology", [](auto& root) { root.erase("mesh_topology"); });
    CheckRejects("zero chordwise panels", [](auto& root) { root["mesh_topology"]["chordwise_panels"] = 0; });
    CheckRejects("a fractional panel count",
                 [](auto& root) { root["mesh_topology"]["chordwise_panels"] = 8.5; });
    CheckRejects("an unknown wake model",
                 [](auto& root) { root["mesh_topology"]["wake_model"] = "freewake"; });

    CheckRejects("a non-positive disk radius",
                 [](auto& root) { root["propulsion_bemt"]["disk_radius"] = 0.0; });
    CheckRejects("a non-positive reference rpm",
                 [](auto& root) { root["propulsion_bemt"]["reference_rpm"] = 0.0; });
    CheckRejects("a negative blade chord",
                 [](auto& root) { root["propulsion_bemt"]["blade_stations"][0]["chord"] = -0.01; });
    CheckRejects("an r/R above the tip",
                 [](auto& root) { root["propulsion_bemt"]["blade_stations"][0]["r_over_R"] = 1.5; });
    CheckRejects("out-of-order blade stations",
                 [](auto& root) { root["propulsion_bemt"]["blade_stations"][1]["r_over_R"] = 0.41; });
}

void TestRejectsMalformedJson() {
    try {
        (void)Aeolion::Geometry::ParseHandoff(nlohmann::json::parse("[]"));
        std::cerr << "FAIL: expected rejection -- a top-level array\n";
        ++g_Failures;
    } catch (const ContractError&) {
    }
}

} // namespace

int main() {
    TestParsesFixture();
    TestSurfaceBinding();
    TestSchema110AgreesWithNameConvention();
    TestToWingParams();
    TestNormalizesHingeAxis();
    TestToleratesUnknownKeys();
    TestOptionalBlocks();
    TestSchema140Blocks();
    TestRejectsMalformed140();
    TestRejectsMalformed();
    TestRejectsMalformedJson();

    if (g_Failures > 0) {
        std::cerr << g_Failures << " check(s) failed\n";
        return 1;
    }
    return 0;
}
