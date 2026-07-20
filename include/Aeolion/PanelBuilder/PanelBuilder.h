// PanelBuilder/PanelBuilder.h
//
// Builds the solver's lattice directly from a parsed geometry handoff --
// station by station, on the CAMBER SURFACE the handoff's CST sections
// describe, discretized both spanwise AND chordwise, with wing-bound
// control surfaces deflected about their stated hinge lines. This is the
// "station-resolved lattice builder" that HandoffContract.h's own comment
// calls out as what supersedes Geometry::ToWingParams(), which can only
// reduce a design to Solver::WingParams's single flat trapezoid.
//
// --- planform: how the quarter-chord line is traced ------------------------
// Geometry::PlanformStation.h states the rule this relies on: Eta, Chord,
// TwistDeg, SweepQuarterChordDeg, and DihedralDeg are all linearly
// interpolated between adjacent stations. Chord and twist are therefore
// exact by direct linear interpolation. Sweep and dihedral are angles, and
// a linearly-varying angle has no closed-form quarter-chord curve (the
// curve integrates tan(angle) over span), so the quarter-chord position is
// accumulated outward from the root one panel boundary at a time using each
// sub-interval's midpoint angle. That accumulation is EXACT whenever
// sweep/dihedral are constant across a station segment (the tangent is then
// genuinely constant and the sum telescopes to the closed form, regardless
// of spanwise_panels_per_section) and converges with panel count otherwise.
//
// --- camber: what makes this not a flat lattice ----------------------------
// Every lattice point is placed on the mean line of the local CST section
// (Geometry/CstSurface.h), interpolated across span between the handoff's
// airfoil_sections. The panel's flow-tangency normal is taken from the
// camber SLOPE at its own three-quarter-chord point, not merely from twist,
// so a cambered section develops lift at zero geometric incidence -- the
// physical effect a flat lattice structurally cannot represent.
//
// Panel::Area is the area of that curved surface, from the true camber-line
// arc length and the dihedral-tilted spanwise extent. Panel::PlanformArea
// is its projection, which is what every force coefficient normalizes by
// (see Solver::Panel for why the two must not be confused).
//
// --- chordwise discretization ----------------------------------------------
// mesh_topology.chordwise_panels rows are laid along each section. Row m
// spans some [psi_a, psi_b]; its bound vortex sits at that row's own quarter
// chord and its control point at that row's three-quarter chord. This is the
// classic horseshoe-per-panel arrangement (Katz & Plotkin), which the
// existing solver core already handles unmodified: it operates on an
// arbitrary panel list and never assumed one row. Every panel in a chordwise
// stack shares a StripIndex, so the solver reports the stack as ONE spanwise
// station with the summed section circulation.
//
// --- control surfaces -------------------------------------------------------
// Where a wing-bound control surface covers a strip, the chordwise rows are
// redistributed so that a panel EDGE falls exactly on the hinge line rather
// than a panel straddling it -- a straddling panel would have to be either
// wholly deflected or wholly not, making the hinge position an artifact of
// the mesh. Rows are then split between the fore and aft portions in
// proportion to their chord fractions.
//
// Panels aft of the hinge are rigidly rotated about the hinge line: both the
// lattice geometry AND the flow-tangency normal, so the deflection is a real
// change in shape rather than a linearized normal-only approximation. Since
// the rotation is rigid it preserves camber-line arc length, so a deflected
// panel's Area is unchanged -- only its projection moves.
//
// Duct-jet-bound surfaces are skipped entirely. ControlSurface.h is
// emphatic that they are not lifting surfaces and must never enter the wing
// lattice; their eta band measures a duct-exit radius fraction, a different
// quantity that merely shares the name.
//
// Hinge axes arrive in the contract's own aetherion_body_frd frame and are
// converted here, at ingest, exactly as ADR-0016 and ControlSurface.h
// prescribe -- rotating the axis vector and keeping the right-hand rule
// about it, never re-deriving a sign convention from prose.
//
// With one chordwise row, a symmetric section and no control surfaces, all
// of this reproduces Solver::BuildWing() exactly, which TestPanelBuilder
// asserts. Two knowing deviations from that builder, both toward being more
// geometrically consistent:
//   - The section frame here is orthonormal: the chordwise axis is twisted
//     about global y, and the surface-normal axis comes from crossing it
//     with the dihedral-tilted spanwise axis. BuildWing() instead rotates
//     the normal by twist and dihedral in sequence while leaving the
//     chordwise offset untwisted by dihedral. The two agree exactly when
//     only one of twist/dihedral is nonzero, and differ by O(sin(twist) *
//     (1 - cos(dihedral))) when both are.
//   - Spanwise panel edges follow the interpolated planform, so a cranked
//     wing gets a genuine kink rather than a fitted straight taper.

#pragma once

#include "Aeolion/Geometry/CstSurface.h"
#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Solver/Panel.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace Aeolion::PanelBuilder {

// A commanded deflection of one control surface, in degrees, positive by
// the right-hand rule about the surface's own hinge axis (for the usual
// spanwise axis that is trailing-edge down).
//
// The two semi-spans are commanded separately because ONE contract entry
// describes the band on BOTH of them -- its eta is a semi-span fraction --
// and the interesting control modes differ in exactly this sign. Equal
// angles give a flap/elevator; opposite angles give an aileron. There is no
// safe default that covers both, so both are stated.
//
// Surfaces are addressed by INDEX into HandoffContract::ControlSurfaces,
// never by name: ControlSurface.h warns that names are not unique (a
// four-vane duct emits four entries all called "vane").
struct ControlDeflection {
    std::size_t SurfaceIndex = 0;
    double RightAngleDeg = 0.0; // +y semi-span
    double LeftAngleDeg = 0.0;  // -y semi-span
};

[[nodiscard]] inline ControlDeflection Symmetric(std::size_t surfaceIndex, double angleDeg) {
    return {surfaceIndex, angleDeg, angleDeg};
}

[[nodiscard]] inline ControlDeflection Antisymmetric(std::size_t surfaceIndex, double angleDeg) {
    return {surfaceIndex, angleDeg, -angleDeg};
}

namespace Detail {

// A chordwise row count below this cannot resolve a hinge line at all --
// there is no interior edge to place it on -- so a covered strip keeps its
// uniform distribution and the surface simply does not articulate.
inline constexpr int MinRowsToResolveHinge = 2;

// One accumulated station on the (implicit right) semi-span: the
// chord-line quarter-chord point plus everything needed to place lattice
// points on the section there. All of these are functions of |y| only, so a
// single semi-span pass is mirrored for both sides rather than walked twice.
struct SpanStation {
    double Eta = 0.0;
    double y = 0.0;
    double QuarterChordX = 0.0;
    double QuarterChordZ = 0.0;
    double Chord = 0.0;
    double TwistDeg = 0.0;
    double DihedralDeg = 0.0;
    double SweepQuarterChordDeg = 0.0; // only needed to seed the next accumulation step
};

[[nodiscard]] inline std::vector<SpanStation> BuildSemiSpanStations(const Geometry::HandoffContract& contract) {
    const auto& stations = contract.Stations;
    const double halfSpan = contract.Span * Math::Half;
    const int panelsPerSegment = contract.Mesh.SpanwisePanelsPerSection;

    std::vector<SpanStation> result;
    result.reserve((stations.size() - 1) * static_cast<std::size_t>(panelsPerSegment) + 1);

    SpanStation previous;
    previous.Eta = stations.front().Eta;
    previous.Chord = stations.front().Chord;
    previous.TwistDeg = stations.front().TwistDeg;
    previous.DihedralDeg = stations.front().DihedralDeg;
    previous.SweepQuarterChordDeg = stations.front().SweepQuarterChordDeg;
    result.push_back(previous);

    for (std::size_t i = 0; i + 1 < stations.size(); ++i) {
        const Geometry::PlanformStation& a = stations[i];
        const Geometry::PlanformStation& b = stations[i + 1];
        for (int j = 1; j <= panelsPerSegment; ++j) {
            const double t = static_cast<double>(j) / panelsPerSegment;

            SpanStation current;
            current.Eta = a.Eta + (b.Eta - a.Eta) * t;
            current.y = current.Eta * halfSpan;
            current.Chord = a.Chord + (b.Chord - a.Chord) * t;
            current.TwistDeg = a.TwistDeg + (b.TwistDeg - a.TwistDeg) * t;
            current.DihedralDeg = a.DihedralDeg + (b.DihedralDeg - a.DihedralDeg) * t;
            current.SweepQuarterChordDeg =
                a.SweepQuarterChordDeg + (b.SweepQuarterChordDeg - a.SweepQuarterChordDeg) * t;

            const double dy = current.y - previous.y;
            const double sweepMidRad =
                Math::DegToRad(Math::Half * (previous.SweepQuarterChordDeg + current.SweepQuarterChordDeg));
            const double dihedralMidRad =
                Math::DegToRad(Math::Half * (previous.DihedralDeg + current.DihedralDeg));
            current.QuarterChordX = previous.QuarterChordX + std::tan(sweepMidRad) * dy;
            current.QuarterChordZ = previous.QuarterChordZ + std::tan(dihedralMidRad) * dy;

            result.push_back(current);
            previous = current;
        }
    }
    return result;
}

[[nodiscard]] inline SpanStation MidStation(const SpanStation& a, const SpanStation& b) {
    SpanStation mid;
    mid.Eta = Math::Half * (a.Eta + b.Eta);
    mid.y = Math::Half * (a.y + b.y);
    mid.QuarterChordX = Math::Half * (a.QuarterChordX + b.QuarterChordX);
    mid.QuarterChordZ = Math::Half * (a.QuarterChordZ + b.QuarterChordZ);
    mid.Chord = Math::Half * (a.Chord + b.Chord);
    mid.TwistDeg = Math::Half * (a.TwistDeg + b.TwistDeg);
    mid.DihedralDeg = Math::Half * (a.DihedralDeg + b.DihedralDeg);
    mid.SweepQuarterChordDeg = Math::Half * (a.SweepQuarterChordDeg + b.SweepQuarterChordDeg);
    return mid;
}

// The section's local orthonormal frame at a station. ChordDir runs aft
// along the (twisted) chord line; UpDir is perpendicular to it in the
// section plane, tilted by dihedral -- this is the direction camber is
// measured along, and the panel normal when the camber slope is zero.
struct SectionFrame {
    Solver::Vec3 ChordDir;
    Solver::Vec3 UpDir;
};

[[nodiscard]] inline SectionFrame FrameAt(const SpanStation& station, double sign) {
    const double twistRad = Math::DegToRad(station.TwistDeg);
    const double dihedralRad = Math::DegToRad(station.DihedralDeg);

    SectionFrame frame;
    frame.ChordDir = Solver::RotateAboutY(Solver::Vec3(1, 0, 0), twistRad);
    // Spanwise reference tilted by dihedral, mirrored for the left side so
    // a positive dihedral raises BOTH tips.
    const Solver::Vec3 spanDir = Solver::RotateAboutX(Solver::Vec3(0, 1, 0), dihedralRad * sign);
    frame.UpDir = Solver::Cross(frame.ChordDir, spanDir).Normalized();
    return frame;
}

// The hinge a given strip articulates about, already resolved to this
// semi-span's axis and angle. Inactive when no wing-bound control surface
// covers the strip.
struct HingeSpec {
    bool Active = false;
    double PsiHinge = 1.0;
    Solver::Vec3 Axis{0, 1, 0};
    double RightAngleRad = 0.0;
    double LeftAngleRad = 0.0;
};

// The contract states its frame once and the consumer converts explicitly
// at ingest (ADR-0016). aetherion_body_frd is x-forward / y-right / z-down;
// the solver is x-aft / y-right / z-up, a 180-degree rotation about y. That
// is a PROPER rotation, so a rotation axis transforms the same way an
// ordinary vector does and the right-hand rule survives intact.
[[nodiscard]] inline Solver::Vec3 HingeAxisToSolverFrame(const Math::Vec3& contractAxis) {
    return Solver::Vec3(-contractAxis.x, contractAxis.y, -contractAxis.z);
}

// Mirror of a rotation axis across the xz-plane. A rotation axis is a
// pseudovector, so mirroring negates the components in the plane and keeps
// the one normal to it -- which is why commanding both semi-spans the SAME
// angle produces a symmetric (flap-like) deflection, and opposite angles
// produce an aileron.
[[nodiscard]] inline Solver::Vec3 MirrorHingeAxis(const Solver::Vec3& axis) {
    return Solver::Vec3(-axis.x, axis.y, -axis.z);
}

[[nodiscard]] inline HingeSpec HingeForStrip(const Geometry::HandoffContract& contract,
                                             const std::vector<ControlDeflection>& deflections, double eta) {
    HingeSpec hinge;
    for (std::size_t i = 0; i < contract.ControlSurfaces.size(); ++i) {
        const Geometry::ControlSurface& surface = contract.ControlSurfaces[i];
        // Duct-jet vanes are not lifting surfaces and their eta means a
        // radius fraction, not a semi-span fraction -- never lattice them.
        if (surface.Binding != Geometry::ControlSurfaceBinding::Wing) continue;
        if (eta < surface.EtaStart || eta > surface.EtaEnd) continue;

        hinge.Active = true;
        hinge.PsiHinge = 1.0 - surface.ChordFraction;
        hinge.Axis = HingeAxisToSolverFrame(surface.HingeAxis);
        for (const ControlDeflection& deflection : deflections) {
            if (deflection.SurfaceIndex != i) continue;
            hinge.RightAngleRad = Math::DegToRad(deflection.RightAngleDeg);
            hinge.LeftAngleRad = Math::DegToRad(deflection.LeftAngleDeg);
            break;
        }
        // First covering surface wins. Overlapping bands over the same
        // chord region are not physically meaningful, and the contract
        // does not rank them.
        break;
    }
    return hinge;
}

// Chordwise row boundaries in psi. Uniform, except that a hinge splits the
// chord into fore and aft groups so that a row EDGE lands exactly on it.
[[nodiscard]] inline std::vector<std::pair<double, double>> ChordwiseRowBounds(int rows, const HingeSpec& hinge) {
    std::vector<std::pair<double, double>> bounds;
    bounds.reserve(static_cast<std::size_t>(rows));

    const bool splittable = hinge.Active && rows >= MinRowsToResolveHinge && hinge.PsiHinge > 0.0 &&
                            hinge.PsiHinge < 1.0;
    if (!splittable) {
        for (int m = 0; m < rows; ++m)
            bounds.emplace_back(static_cast<double>(m) / rows, static_cast<double>(m + 1) / rows);
        return bounds;
    }

    int aft = static_cast<int>(std::lround(rows * (1.0 - hinge.PsiHinge)));
    aft = std::clamp(aft, 1, rows - 1);
    const int fore = rows - aft;

    for (int m = 0; m < fore; ++m)
        bounds.emplace_back(hinge.PsiHinge * m / fore, hinge.PsiHinge * (m + 1) / fore);
    for (int m = 0; m < aft; ++m)
        bounds.emplace_back(hinge.PsiHinge + (1.0 - hinge.PsiHinge) * m / aft,
                            hinge.PsiHinge + (1.0 - hinge.PsiHinge) * (m + 1) / aft);
    return bounds;
}

// A point on the camber surface at chordwise fraction psi. The planform's
// quarter-chord line is a CHORD-LINE reference (psi = 1/4), so the chordwise
// offset is measured from there; camber is then added along the section
// normal, scaled by the local chord.
[[nodiscard]] inline Solver::Vec3 SurfacePoint(const SpanStation& station, double sign, double psi,
                                               const std::vector<Geometry::AirfoilSection>& sections) {
    const SectionFrame frame = FrameAt(station, sign);
    const Solver::Vec3 quarterChordPoint(station.QuarterChordX, sign * station.y, station.QuarterChordZ);
    const double camber = Geometry::CamberAt(sections, station.Eta, psi);
    return quarterChordPoint + frame.ChordDir * ((psi - Math::QuarterChord) * station.Chord) +
           frame.UpDir * (camber * station.Chord);
}

} // namespace Detail

// Builds the full-span wing lattice from the handoff's per-station planform
// and CST camber surface, discretized spanwise and chordwise, with the given
// control surfaces deflected. An empty deflection list still articulates the
// mesh around every hinge line (at zero angle), so the lattice topology does
// not change when a surface is later commanded.
[[nodiscard]] inline std::vector<Solver::Panel> BuildPanels(const Geometry::HandoffContract& contract,
                                                            const std::vector<ControlDeflection>& deflections = {}) {
    const std::vector<Detail::SpanStation> half = Detail::BuildSemiSpanStations(contract);
    const auto& sections = contract.AirfoilSections;
    const int rows = contract.Mesh.ChordwisePanels;

    std::vector<Solver::Panel> panels;
    panels.reserve((half.size() - 1) * 2 * static_cast<std::size_t>(rows));

    int nextStripIndex = 0;

    // `inner`/`outer` are two adjacent semi-span stations (inner has the
    // smaller |y|); `sign` is +1 for the right semi-span, -1 for the
    // mirrored left one.
    auto emitStrip = [&](const Detail::SpanStation& inner, const Detail::SpanStation& outer, double sign) {
        const Detail::SpanStation mid = Detail::MidStation(inner, outer);
        const Detail::SectionFrame midFrame = Detail::FrameAt(mid, sign);
        const Detail::HingeSpec hinge = Detail::HingeForStrip(contract, deflections, mid.Eta);
        const auto rowBounds = Detail::ChordwiseRowBounds(rows, hinge);
        const int stripIndex = nextStripIndex++;

        const Solver::Vec3 hingeAxis =
            (sign > 0) ? hinge.Axis : Detail::MirrorHingeAxis(hinge.Axis);
        const double hingeAngle = (sign > 0) ? hinge.RightAngleRad : hinge.LeftAngleRad;
        const bool articulates = hinge.Active && hingeAngle != 0.0;

        // Rigid rotation of the aft portion about the hinge line, taken at
        // this point's own span station so the hinge follows the planform.
        auto deflect = [&](const Detail::SpanStation& station, const Solver::Vec3& point, double psi) {
            if (!articulates || psi <= hinge.PsiHinge) return point;
            const Solver::Vec3 hingePoint = Detail::SurfacePoint(station, sign, hinge.PsiHinge, sections);
            return hingePoint + Solver::RotateAboutAxis(point - hingePoint, hingeAxis, hingeAngle);
        };

        for (int m = 0; m < rows; ++m) {
            const auto [psiStart, psiEnd] = rowBounds[static_cast<std::size_t>(m)];
            const double psiBound = psiStart + Math::QuarterChord * (psiEnd - psiStart);
            const double psiControl = psiStart + Math::ThreeQuarterChord * (psiEnd - psiStart);

            const Solver::Vec3 innerPoint =
                deflect(inner, Detail::SurfacePoint(inner, sign, psiBound, sections), psiBound);
            const Solver::Vec3 outerPoint =
                deflect(outer, Detail::SurfacePoint(outer, sign, psiBound, sections), psiBound);

            Solver::Panel panel;
            // Panel.h requires A.y < B.y; the right semi-span walks
            // root->tip with y increasing, the mirrored left one the
            // opposite way.
            if (sign > 0) { panel.A = innerPoint; panel.B = outerPoint; }
            else          { panel.A = outerPoint; panel.B = innerPoint; }

            panel.ControlPoint =
                deflect(mid, Detail::SurfacePoint(mid, sign, psiControl, sections), psiControl);

            // Flow tangency is applied against the camber surface: tilt the
            // section normal by the local mean-line slope. d(camber)/d(psi)
            // is already normalized by chord, so it IS the slope dz/dx. The
            // normal then rides the hinge rotation with the geometry.
            const double slope = Geometry::CamberSlopeAt(sections, mid.Eta, psiControl);
            Solver::Vec3 normal = (midFrame.UpDir - midFrame.ChordDir * slope).Normalized();
            if (articulates && psiControl > hinge.PsiHinge)
                normal = Solver::RotateAboutAxis(normal, hingeAxis, hingeAngle);
            panel.Normal = normal;

            panel.TrailDirA = Solver::Vec3(1, 0, 0);
            panel.TrailDirB = Solver::Vec3(1, 0, 0);
            panel.SpanwiseWidth = outer.y - inner.y;

            // Projection for coefficients; true curved surface for geometry.
            // Deflection is a rigid rotation, so it leaves arc length -- and
            // therefore Area -- untouched.
            const double arcFraction =
                Geometry::CamberArcLengthFraction(sections, mid.Eta, psiStart, psiEnd);
            const double dihedralRad = Math::DegToRad(mid.DihedralDeg);
            panel.PlanformArea = (mid.Chord * (psiEnd - psiStart)) * panel.SpanwiseWidth;
            panel.Area = (mid.Chord * arcFraction) * (panel.SpanwiseWidth / std::cos(dihedralRad));

            panel.Surface = "wing";
            panel.StripIndex = stripIndex;
            panels.push_back(panel);
        }
    };

    for (std::size_t k = 0; k + 1 < half.size(); ++k) emitStrip(half[k], half[k + 1], -1.0); // left semi-span
    for (std::size_t k = 0; k + 1 < half.size(); ++k) emitStrip(half[k], half[k + 1], +1.0); // right semi-span

    std::ranges::sort(panels, {}, [](const Solver::Panel& p) { return p.A.y; });
    return panels;
}

} // namespace Aeolion::PanelBuilder
