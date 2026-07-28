// PanelBuilder/PanelBuilder.cpp
//
// Implementation of LatticeBuilder. The rationale for what this code does
// -- the quarter-chord march, the camber surface, the breakpoint rules, the
// hinge treatment, and the measurements behind the spacing default -- lives
// in the header alongside the declarations it explains.

#include "Aeolion/PanelBuilder/PanelBuilder.h"

#include "Aeolion/Geometry/CstSurface.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

namespace Aeolion::PanelBuilder {

// Momentum-theory hover inflow ratio, lambda = v_i / (Omega * R): the
// prescribed-wake axial-convection floor BuildPropellerLattice() gives its
// trailing legs at low advance speed (see the comment at its use).
// lambda = sqrt(C_T / 2) ~ 0.07 for ordinary propeller disk loadings.
inline constexpr double HoverInflowRatio = 0.07;

LatticeBuilder::LatticeBuilder(Geometry::HandoffContract contract, LatticeOptions options)
    : m_Contract(std::move(contract)), m_Options(options) {
    // Breakpoints and the spanwise march depend only on the contract and
    // the spacing choice, never on a commanded deflection, so they are
    // computed once here rather than per Build().
    // Order matters. The trim station is a spanwise breakpoint, so it has to
    // exist before the boundary etas are laid out; the placement offset is a
    // rigid translation applied afterwards and affects neither.
    m_TrimEta = ComputeTrimEta();
    m_BoundaryEtas = ComputeBoundaryEtas();
    m_SemiSpan = ComputeSemiSpan();
    m_PlacementOffset = ComputePlacementOffset();
}

// Where the wing enters the fuselage, as a semi-span fraction. Zero when
// there is nothing to trim against.
double LatticeBuilder::ComputeTrimEta() const {
    if (!m_Options.TrimWingAtBody) return 0.0;
    if (!m_Contract.Body.IsPresent() || !m_Contract.Placement.IsStated) return 0.0;

    const double halfSpan = m_Contract.Span * Math::Half;
    if (halfSpan <= 0.0) return 0.0;

    // The anchor is already in the contract's frame, which is the frame the
    // body's radius law is stated in -- so no conversion here.
    const double radius = Geometry::RadiusAt(m_Contract.Body, m_Contract.Placement.RootLeadingEdge.x);
    const double eta = radius / halfSpan;
    // A body wider than the wing would leave nothing to solve; refuse to
    // trim it all away rather than return an empty lattice.
    return (eta < MaxTrimEta) ? eta : 0.0;
}

// Rigid translation putting the root leading edge where the contract says.
// The lattice is otherwise built with the root quarter chord at the origin.
Solver::Vec3 LatticeBuilder::ComputePlacementOffset() const {
    if (!m_Contract.Placement.IsStated || m_SemiSpan.empty()) return Solver::Vec3(0, 0, 0);

    const Math::Vec3& anchor = m_Contract.Placement.RootLeadingEdge;
    // Contract frame is x-forward / z-down, the solver x-aft / z-up: the
    // same 180-degree rotation about y used for hinge axes and body stations.
    const Solver::Vec3 target(-anchor.x, anchor.y, -anchor.z);

    // Where the untranslated root leading edge currently sits. Camber
    // vanishes at psi = 0, so the leading edge is on the chord line.
    const SpanStation& root = m_SemiSpan.front();
    const SectionFrame frame = FrameAt(root, 1.0);
    const Solver::Vec3 current = Solver::Vec3(root.QuarterChordX, 0.0, root.QuarterChordZ) -
                                 frame.ChordDir * (Math::QuarterChord * root.Chord);
    return target - current;
}

double LatticeBuilder::GrossPlanformArea() const {
    const auto& stations = m_Contract.Stations;
    const double halfSpan = m_Contract.Span * Math::Half;
    double half = 0.0;
    for (std::size_t i = 0; i + 1 < stations.size(); ++i)
        half += Math::Half * (stations[i].Chord + stations[i + 1].Chord) *
                (stations[i + 1].Eta - stations[i].Eta) * halfSpan;
    return 2.0 * half; // both semi-spans
}

LatticeBuilder& LatticeBuilder::Deflect(const ControlDeflection& deflection) {
    for (ControlDeflection& existing : m_Deflections) {
        if (existing.SurfaceIndex == deflection.SurfaceIndex) {
            existing = deflection;
            return *this;
        }
    }
    m_Deflections.push_back(deflection);
    return *this;
}

LatticeBuilder& LatticeBuilder::ClearDeflections() {
    m_Deflections.clear();
    return *this;
}

// --- pure geometry helpers --------------------------------------------------

double LatticeBuilder::SpacingFraction(Spacing spacing, int index, int count) {
    const double t = static_cast<double>(index) / count;
    if (spacing == Spacing::Uniform) return t;
    return Math::Half * (1.0 - std::cos(std::numbers::pi * t));
}

// Split `total` among sub-intervals in proportion to their widths, giving
// each at least `minimum`, and summing to exactly `total` (largest remainder
// takes the leftovers). This honours mesh_topology's
// "spanwise_panels_per_section" literally: subdividing a section at a
// control surface edge redistributes that section's panels rather than
// adding more, so panel count -- and the O(N^3) solve cost -- stays
// predictable no matter how many surfaces the contract carries.
std::vector<int> LatticeBuilder::AllocateProportional(const std::vector<double>& widths, int total,
                                                             int minimum) {
    const std::size_t count = widths.size();
    std::vector<int> allocation(count, minimum);
    const int assigned = minimum * static_cast<int>(count);
    if (assigned >= total) return allocation; // the minimums already exhaust the budget

    double totalWidth = 0.0;
    for (double width : widths) totalWidth += width;
    const int remaining = total - assigned;

    std::vector<double> remainder(count, 0.0);
    int handed = 0;
    for (std::size_t i = 0; i < count; ++i) {
        const double ideal = (totalWidth > 0.0) ? remaining * widths[i] / totalWidth : 0.0;
        const int whole = static_cast<int>(ideal);
        allocation[i] += whole;
        handed += whole;
        remainder[i] = ideal - whole;
    }
    for (int leftover = remaining - handed; leftover > 0; --leftover) {
        std::size_t best = 0;
        for (std::size_t i = 1; i < count; ++i)
            if (remainder[i] > remainder[best]) best = i;
        ++allocation[best];
        remainder[best] = -1.0; // spent
    }
    return allocation;
}

// Planform quantities at an arbitrary eta, linearly interpolated between the
// bracketing contract stations (Geometry::PlanformStation.h states that this
// linear rule is the contract's own).
Geometry::PlanformStation LatticeBuilder::InterpolatePlanform(
    const std::vector<Geometry::PlanformStation>& stations, double eta) {
    std::size_t high = 1;
    while (high + 1 < stations.size() && stations[high].Eta < eta) ++high;
    const Geometry::PlanformStation& a = stations[high - 1];
    const Geometry::PlanformStation& b = stations[high];
    const double width = b.Eta - a.Eta;
    const double t = (width > 0.0) ? std::clamp((eta - a.Eta) / width, 0.0, 1.0) : 0.0;

    Geometry::PlanformStation result;
    result.Eta = eta;
    result.Chord = a.Chord + (b.Chord - a.Chord) * t;
    result.TwistDeg = a.TwistDeg + (b.TwistDeg - a.TwistDeg) * t;
    result.SweepQuarterChordDeg =
        a.SweepQuarterChordDeg + (b.SweepQuarterChordDeg - a.SweepQuarterChordDeg) * t;
    result.DihedralDeg = a.DihedralDeg + (b.DihedralDeg - a.DihedralDeg) * t;
    return result;
}

// The contract states its frame once and the consumer converts explicitly
// at ingest (ADR-0016). aetherion_body_frd is x-forward / y-right / z-down;
// the solver is x-aft / y-right / z-up, a 180-degree rotation about y. That
// is a PROPER rotation, so a rotation axis transforms the same way an
// ordinary vector does and the right-hand rule survives intact.
Solver::Vec3 LatticeBuilder::HingeAxisToSolverFrame(const Math::Vec3& contractAxis) {
    return Solver::Vec3(-contractAxis.x, contractAxis.y, -contractAxis.z);
}

// Mirror of a rotation axis across the xz-plane. A rotation axis is a
// pseudovector, so mirroring negates the components in the plane and keeps
// the one normal to it -- which is why commanding both semi-spans the SAME
// angle produces a symmetric (flap-like) deflection, and opposite angles
// produce an aileron.
Solver::Vec3 LatticeBuilder::MirrorHingeAxis(const Solver::Vec3& axis) {
    return Solver::Vec3(-axis.x, axis.y, -axis.z);
}

LatticeBuilder::SectionFrame LatticeBuilder::FrameAt(const SpanStation& station, double sign) {
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

LatticeBuilder::SpanStation LatticeBuilder::MidStation(const SpanStation& a, const SpanStation& b) {
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

// --- construction-time work -------------------------------------------------

// Every eta at which the mesh MUST place a panel boundary: the planform
// stations, plus the edges of every wing-bound control surface.
std::vector<double> LatticeBuilder::ComputeBreakpoints() const {
    std::vector<double> points;
    points.reserve(m_Contract.Stations.size() + 2 * m_Contract.ControlSurfaces.size());
    // The parser guarantees these are sorted, unique, and span [0,1].
    for (const Geometry::PlanformStation& station : m_Contract.Stations) points.push_back(station.Eta);

    const auto insertIfDistinct = [&](double eta) {
        if (eta <= EtaInteriorMin || eta >= EtaInteriorMax) return; // 0 and 1 are already stations
        for (double existing : points)
            if (std::fabs(existing - eta) < BreakpointMergeTolerance) return;
        points.push_back(eta);
    };
    for (const Geometry::ControlSurface& surface : m_Contract.ControlSurfaces) {
        if (surface.Binding != Geometry::ControlSurfaceBinding::Wing) continue;
        insertIfDistinct(surface.EtaStart);
        insertIfDistinct(surface.EtaEnd);
    }
    // The fuselage cut is a breakpoint like any other: a panel edge must
    // land exactly on it, or the trim would fall wherever the mesh happened
    // to put a boundary -- the same artifact control surface edges had.
    insertIfDistinct(m_TrimEta);
    std::ranges::sort(points);
    return points;
}

// The ordered panel boundary etas across one semi-span: each planform
// section is subdivided at any control surface edge inside it, that
// section's panel budget is split proportionally between the pieces, and
// each piece is then filled according to the requested spacing.
std::vector<double> LatticeBuilder::ComputeBoundaryEtas() const {
    const auto& stations = m_Contract.Stations;
    const std::vector<double> breakpoints = ComputeBreakpoints();
    const int panelsPerSection = m_Contract.Mesh.SpanwisePanelsPerSection;

    std::vector<double> etas;
    etas.push_back(stations.front().Eta);

    for (std::size_t i = 0; i + 1 < stations.size(); ++i) {
        const double sectionStart = stations[i].Eta;
        const double sectionEnd = stations[i + 1].Eta;

        // Sub-interval boundaries: this section's ends plus any breakpoint
        // strictly inside it.
        std::vector<double> edges{sectionStart};
        for (double point : breakpoints)
            if (point > sectionStart + BreakpointMergeTolerance &&
                point < sectionEnd - BreakpointMergeTolerance)
                edges.push_back(point);
        edges.push_back(sectionEnd);

        std::vector<double> widths;
        widths.reserve(edges.size() - 1);
        for (std::size_t e = 0; e + 1 < edges.size(); ++e) widths.push_back(edges[e + 1] - edges[e]);

        const std::vector<int> allocation =
            AllocateProportional(widths, panelsPerSection, MinPanelsPerSubInterval);
        for (std::size_t e = 0; e + 1 < edges.size(); ++e) {
            const int count = allocation[e];
            for (int k = 1; k <= count; ++k)
                etas.push_back(edges[e] + widths[e] * SpacingFraction(m_Options.Spanwise, k, count));
        }
    }
    return etas;
}

std::vector<LatticeBuilder::SpanStation> LatticeBuilder::ComputeSemiSpan() const {
    const auto& stations = m_Contract.Stations;
    const double halfSpan = m_Contract.Span * Math::Half;

    std::vector<SpanStation> result;
    result.reserve(m_BoundaryEtas.size());

    const auto stationAt = [&](double eta) {
        const Geometry::PlanformStation planform = InterpolatePlanform(stations, eta);
        SpanStation station;
        station.Eta = eta;
        station.y = eta * halfSpan;
        station.Chord = planform.Chord;
        station.TwistDeg = planform.TwistDeg;
        station.DihedralDeg = planform.DihedralDeg;
        station.SweepQuarterChordDeg = planform.SweepQuarterChordDeg;
        return station;
    };

    SpanStation previous = stationAt(m_BoundaryEtas.front());
    result.push_back(previous);

    for (std::size_t i = 1; i < m_BoundaryEtas.size(); ++i) {
        SpanStation current = stationAt(m_BoundaryEtas[i]);

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
    return result;
}

// --- per-Build work ---------------------------------------------------------

LatticeBuilder::HingeSpec LatticeBuilder::HingeForStrip(double eta) const {
    HingeSpec hinge;
    for (std::size_t i = 0; i < m_Contract.ControlSurfaces.size(); ++i) {
        const Geometry::ControlSurface& surface = m_Contract.ControlSurfaces[i];
        // Duct-jet vanes are not lifting surfaces and their eta means a
        // radius fraction, not a semi-span fraction -- never lattice them.
        if (surface.Binding != Geometry::ControlSurfaceBinding::Wing) continue;
        if (eta < surface.EtaStart || eta > surface.EtaEnd) continue;

        hinge.Active = true;
        hinge.PsiHinge = 1.0 - surface.ChordFraction;
        hinge.Axis = HingeAxisToSolverFrame(surface.HingeAxis);
        for (const ControlDeflection& deflection : m_Deflections) {
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

// Chordwise row boundaries in psi. A hinge splits the chord into fore and
// aft groups so a row EDGE lands exactly on it; each group is then filled
// according to the requested spacing.
std::vector<std::pair<double, double>> LatticeBuilder::ChordwiseRowBounds(const HingeSpec& hinge) const {
    const int rows = m_Contract.Mesh.ChordwisePanels;
    std::vector<std::pair<double, double>> bounds;
    bounds.reserve(static_cast<std::size_t>(rows));

    const auto fill = [&](double start, double end, int count) {
        for (int m = 0; m < count; ++m)
            bounds.emplace_back(start + (end - start) * SpacingFraction(m_Options.Chordwise, m, count),
                                start + (end - start) * SpacingFraction(m_Options.Chordwise, m + 1, count));
    };

    const bool splittable = hinge.Active && rows >= MinRowsToResolveHinge && hinge.PsiHinge > 0.0 &&
                            hinge.PsiHinge < 1.0;
    if (!splittable) {
        fill(0.0, 1.0, rows);
        return bounds;
    }

    // Same proportional split as the spanwise budget: the requested row
    // count is divided between fore and aft, not added to.
    int aft = static_cast<int>(std::lround(rows * (1.0 - hinge.PsiHinge)));
    aft = std::clamp(aft, 1, rows - 1);
    fill(0.0, hinge.PsiHinge, rows - aft);
    fill(hinge.PsiHinge, 1.0, aft);
    return bounds;
}

// A point on the camber surface at chordwise fraction psi. The planform's
// quarter-chord line is a CHORD-LINE reference (psi = 1/4), so the chordwise
// offset is measured from there; camber is then added along the section
// normal, scaled by the local chord.
Solver::Vec3 LatticeBuilder::SurfacePoint(const SpanStation& station, double sign, double psi) const {
    const SectionFrame frame = FrameAt(station, sign);
    const Solver::Vec3 quarterChordPoint(station.QuarterChordX, sign * station.y, station.QuarterChordZ);
    const double camber = Geometry::CamberAt(m_Contract.AirfoilSections, station.Eta, psi);
    return quarterChordPoint + frame.ChordDir * ((psi - Math::QuarterChord) * station.Chord) +
           frame.UpDir * (camber * station.Chord) + m_PlacementOffset;
}

// `inner`/`outer` are two adjacent semi-span stations (inner has the smaller
// |y|); `sign` is +1 for the right semi-span, -1 for the mirrored left one.
void LatticeBuilder::EmitStrip(const SpanStation& inner, const SpanStation& outer, double sign,
                               int stripIndex, bool carryThrough,
                               std::vector<Solver::Panel>& panels) const {
    const SpanStation mid = MidStation(inner, outer);
    const SectionFrame midFrame = FrameAt(mid, sign);
    const HingeSpec hinge = HingeForStrip(mid.Eta);
    const auto rowBounds = ChordwiseRowBounds(hinge);

    const Solver::Vec3 hingeAxis = (sign > 0) ? hinge.Axis : MirrorHingeAxis(hinge.Axis);
    const double hingeAngle = (sign > 0) ? hinge.RightAngleRad : hinge.LeftAngleRad;
    const bool articulates = hinge.Active && hingeAngle != 0.0;

    // Rigid rotation of the aft portion about the hinge line, taken at this
    // point's own span station so the hinge follows the planform.
    const auto deflect = [&](const SpanStation& station, const Solver::Vec3& point, double psi) {
        if (!articulates || psi <= hinge.PsiHinge) return point;
        const Solver::Vec3 hingePoint = SurfacePoint(station, sign, hinge.PsiHinge);
        return hingePoint + Solver::RotateAboutAxis(point - hingePoint, hingeAxis, hingeAngle);
    };

    for (std::size_t m = 0; m < rowBounds.size(); ++m) {
        const auto [psiStart, psiEnd] = rowBounds[m];
        const double psiBound = psiStart + Math::QuarterChord * (psiEnd - psiStart);
        const double psiControl = psiStart + Math::ThreeQuarterChord * (psiEnd - psiStart);

        // On a carry-through strip the bound segment starts at the
        // CENTRELINE, not at the body surface: that extension is the
        // vorticity the wing carries through the fuselage. It follows the
        // planform's own quarter-chord line inboard, so a swept wing's
        // carry-through is swept too.
        const SpanStation& inboard = carryThrough ? m_SemiSpan.front() : inner;
        const Solver::Vec3 innerPoint = deflect(inboard, SurfacePoint(inboard, sign, psiBound), psiBound);
        const Solver::Vec3 outerPoint = deflect(outer, SurfacePoint(outer, sign, psiBound), psiBound);

        Solver::Panel panel;
        // Panel.h requires A.y < B.y; the right semi-span walks root->tip
        // with y increasing, the mirrored left one the opposite way.
        if (sign > 0) { panel.A = innerPoint; panel.B = outerPoint; }
        else          { panel.A = outerPoint; panel.B = innerPoint; }

        panel.ControlPoint = deflect(mid, SurfacePoint(mid, sign, psiControl), psiControl);

        // Flow tangency is applied against the camber surface: tilt the
        // section normal by the local mean-line slope. d(camber)/d(psi) is
        // already normalized by chord, so it IS the slope dz/dx. The normal
        // then rides the hinge rotation with the geometry.
        const double slope = Geometry::CamberSlopeAt(m_Contract.AirfoilSections, mid.Eta, psiControl);
        Solver::Vec3 normal = (midFrame.UpDir - midFrame.ChordDir * slope).Normalized();
        if (articulates && psiControl > hinge.PsiHinge)
            normal = Solver::RotateAboutAxis(normal, hingeAxis, hingeAngle);
        panel.Normal = normal;

        panel.TrailDirA = Solver::Vec3(1, 0, 0);
        panel.TrailDirB = Solver::Vec3(1, 0, 0);
        // Width matches the bound segment's extent, so lift-per-span and
        // sectional cl stay consistent with the load actually on it.
        panel.SpanwiseWidth = outer.y - inboard.y;

        // Projection for coefficients; true curved surface for geometry.
        // Deflection is a rigid rotation, so it leaves arc length -- and
        // therefore Area -- untouched.
        const double arcFraction =
            Geometry::CamberArcLengthFraction(m_Contract.AirfoilSections, mid.Eta, psiStart, psiEnd);
        const double dihedralRad = Math::DegToRad(mid.DihedralDeg);
        // Mean chord over the strip's full extent, which on a carry-through
        // strip reaches the centreline rather than the body surface.
        const double stripChord = carryThrough ? Math::Half * (inboard.Chord + outer.Chord) : mid.Chord;
        panel.PlanformArea = (stripChord * (psiEnd - psiStart)) * panel.SpanwiseWidth;
        panel.Area = (stripChord * arcFraction) * (panel.SpanwiseWidth / std::cos(dihedralRad));

        panel.Surface = "wing";
        panel.StripIndex = stripIndex;
        panels.push_back(panel);
    }
}

std::vector<Solver::Panel> LatticeBuilder::Build() const {
    const int rows = m_Contract.Mesh.ChordwisePanels;
    std::vector<Solver::Panel> panels;
    panels.reserve((m_SemiSpan.size() - 1) * 2 * static_cast<std::size_t>(rows));

    int stripIndex = 0;
    // Strips inboard of the fuselage cut are dropped WHOLE. Trimming part
    // of a chordwise stack would leave the solver reporting a sectional cl
    // over a section chord that no longer exists (see Panel::StripIndex).
    const auto buried = [this](const SpanStation& inner, const SpanStation& outer) {
        return m_TrimEta > 0.0 && Math::Half * (inner.Eta + outer.Eta) < m_TrimEta;
    };

    // The first strip to survive the trim on each side is the one whose
    // bound vortex carries through the body.
    std::size_t innermost = m_SemiSpan.size();
    for (std::size_t k = 0; k + 1 < m_SemiSpan.size(); ++k)
        if (!buried(m_SemiSpan[k], m_SemiSpan[k + 1])) { innermost = k; break; }
    const bool carryThrough = m_Options.CarryThroughLift && m_TrimEta > 0.0;

    for (std::size_t k = 0; k + 1 < m_SemiSpan.size(); ++k) { // left semi-span
        if (buried(m_SemiSpan[k], m_SemiSpan[k + 1])) continue;
        EmitStrip(m_SemiSpan[k], m_SemiSpan[k + 1], -1.0, stripIndex++, carryThrough && k == innermost,
                  panels);
    }
    for (std::size_t k = 0; k + 1 < m_SemiSpan.size(); ++k) { // right semi-span
        if (buried(m_SemiSpan[k], m_SemiSpan[k + 1])) continue;
        EmitStrip(m_SemiSpan[k], m_SemiSpan[k + 1], +1.0, stripIndex++, carryThrough && k == innermost,
                  panels);
    }

    std::ranges::sort(panels, {}, [](const Solver::Panel& p) { return p.A.y; });
    return panels;
}


// --- fuselage ---------------------------------------------------------------
// The body is a surface of revolution: sweep the contract's (x, radius)
// profile around the body axis, one quad per (axial segment x azimuthal
// sector). Sources rather than horseshoes -- see Lattice::SourcePanel for
// why a non-lifting closed volume takes one and not the other.
//
// FRAME. The contract states x FORWARD (aetherion_body_frd) so its stations
// run from a nose at x = 0 to a tail at x = -Length. The solver is x AFT.
// The conversion is the 180-degree rotation about y that ControlSurface.h
// describes: x_solver = -x_frd, y unchanged, z_solver = -z_frd. For an
// axisymmetric profile the z flip is invisible (a circle maps to itself),
// so in practice only the axial coordinate reverses -- but it is written out
// rather than assumed, because "only x changes" is true of THIS body and not
// of the frame.
//
// WINDING. The kernel's solid angle is winding-sensitive, so each quad's
// corner order is fixed against its own outward normal rather than trusted
// to come out right from the sweep direction. Getting this backwards would
// invert the body's whole pressure field, which is exactly the sort of error
// that still looks plausible in a plot.
//
// THE BASE. This airframe's profile does not close: its final radius is a
// duct exit. A source-panel body with an open end does not conserve mass, so
// the base is capped here and the cap is tagged separately, leaving the
// panels in place for the efflux boundary condition to claim later. Until
// then the cap is an ordinary solid wall, which is a closed-body answer --
// honest, and wrong in a known direction rather than unknown.
std::vector<Lattice::SourcePanel> LatticeBuilder::BuildBody() const {
    std::vector<Lattice::SourcePanel> panels;
    if (!m_Options.IncludeBody) return panels;

    const Geometry::BodyGeometry& body = m_Contract.Body;
    if (!body.IsPresent()) return panels;

    const int sectors = m_Options.BodyCircumferentialPanels;
    if (sectors < MinBodySectors) return panels;

    const auto& stations = body.Stations;

    // Contract frame -> solver frame.
    const auto axialPosition = [](double contractX) { return -contractX; };

    // A point on the surface of revolution at axial station value `x`
    // (already in solver frame) and radius `r`, at azimuth `angle`.
    const auto surfacePoint = [](double x, double r, double angle) {
        return Solver::Vec3(x, r * std::cos(angle), r * std::sin(angle));
    };

    // Fixes corner order so the quad's own normal points away from the body
    // axis, then fills in the derived quantities.
    const auto finishPanel = [&](Lattice::SourcePanel& panel, const Solver::Vec3& outwardHint) {
        Solver::Vec3 normal = Solver::Cross(panel.Corners[2] - panel.Corners[0],
                                            panel.Corners[3] - panel.Corners[1]);
        if (Solver::Dot(normal, outwardHint) < 0.0) {
            std::swap(panel.Corners[1], panel.Corners[3]);
            normal = normal * -1.0;
        }
        // Half the diagonal cross product is the quad's area, and stays
        // correct when a ring degenerates to triangles at a closed end.
        panel.Area = Math::Half * Solver::Cross(panel.Corners[2] - panel.Corners[0],
                                                panel.Corners[3] - panel.Corners[1])
                                      .Norm();
        panel.Normal = normal.Normalized();
    };

    // --- lateral surface ----------------------------------------------------
    for (std::size_t i = 0; i + 1 < stations.size(); ++i) {
        const double xNear = axialPosition(stations[i].x);
        const double xFar = axialPosition(stations[i + 1].x);
        const double rNear = stations[i].Radius;
        const double rFar = stations[i + 1].Radius;

        // A zero-length segment carries no surface.
        if (std::fabs(xFar - xNear) < BodyDegenerateLength && std::fabs(rFar - rNear) < BodyDegenerateLength)
            continue;

        for (int j = 0; j < sectors; ++j) {
            const double angle0 = Math::Two * std::numbers::pi * j / sectors;
            const double angle1 = Math::Two * std::numbers::pi * (j + 1) / sectors;

            Lattice::SourcePanel panel;
            panel.Corners = {surfacePoint(xNear, rNear, angle0), surfacePoint(xFar, rFar, angle0),
                             surfacePoint(xFar, rFar, angle1), surfacePoint(xNear, rNear, angle1)};

            Solver::Vec3 centroid(0, 0, 0);
            for (const Solver::Vec3& corner : panel.Corners) centroid = centroid + corner * 0.25;

            // Radially outward at the panel's own azimuth; degenerate on the
            // axis, where the ring has collapsed to a point.
            const double midAngle = Math::Half * (angle0 + angle1);
            const Solver::Vec3 outward(0.0, std::cos(midAngle), std::sin(midAngle));
            finishPanel(panel, outward);

            // Put the control point on the true frustum rather than inside
            // the chord the flat quad cuts across the circle.
            const double midRadius = Math::Half * (rNear + rFar);
            panel.ControlPoint = (midRadius > BodyDegenerateLength)
                                     ? Solver::Vec3(centroid.x, midRadius * std::cos(midAngle),
                                                    midRadius * std::sin(midAngle))
                                     : centroid;
            panel.Surface = BodySurfaceName;
            panel.StationIndex = static_cast<int>(i);
            if (panel.Area > BodyDegenerateArea) panels.push_back(panel);
        }
    }

    // --- base cap -----------------------------------------------------------
    // Only when the profile actually ends open.
    const double baseRadius = stations.back().Radius;
    if (baseRadius > BodyDegenerateLength) {
        const double baseX = axialPosition(stations.back().x);
        const Solver::Vec3 centre(baseX, 0.0, 0.0);
        // The base faces aft, which in solver axes is +x.
        const Solver::Vec3 outward(1.0, 0.0, 0.0);

        // Subdivided radially, not a single fan of full-radius triangles.
        // Potential flow has a velocity singularity at the sharp base rim,
        // and one control point per full-radius wedge resolves it so badly
        // that a CLOSED body stops satisfying d'Alembert -- which is how
        // this was caught.
        const int rings = BodyBaseRadialRings;
        for (int ring = 0; ring < rings; ++ring) {
            const double rInner = baseRadius * ring / rings;
            const double rOuter = baseRadius * (ring + 1) / rings;
            for (int j = 0; j < sectors; ++j) {
                const double angle0 = Math::Two * std::numbers::pi * j / sectors;
                const double angle1 = Math::Two * std::numbers::pi * (j + 1) / sectors;

                Lattice::SourcePanel panel;
                if (ring == 0) {
                    // Innermost ring is a triangle, written as a quad with a
                    // repeated apex; the kernel skips the zero-length edge.
                    panel.Corners = {centre, surfacePoint(baseX, rOuter, angle0),
                                     surfacePoint(baseX, rOuter, angle1), centre};
                } else {
                    panel.Corners = {surfacePoint(baseX, rInner, angle0), surfacePoint(baseX, rOuter, angle0),
                                     surfacePoint(baseX, rOuter, angle1), surfacePoint(baseX, rInner, angle1)};
                }
                finishPanel(panel, outward);

                const double midAngle = Math::Half * (angle0 + angle1);
                const double centroidRadius =
                    (ring == 0) ? Two3rds * rOuter : Math::Half * (rInner + rOuter);
                panel.ControlPoint = Solver::Vec3(baseX, centroidRadius * std::cos(midAngle),
                                                  centroidRadius * std::sin(midAngle));
                panel.Surface = BodyBaseSurfaceName;
                panel.StationIndex = static_cast<int>(stations.size()) - 1;
                if (panel.Area > BodyDegenerateArea) panels.push_back(panel);
            }
        }
    }

    return panels;
}


// --- duct ---------------------------------------------------------------
// The duct is disjoint from the fuselage (see Geometry/DuctGeometry.h) and,
// unlike it, is hollow all the way through: the bore is where the propeller
// and its slipstream sit, so nothing caps it. What DOES need panels is the
// material of the ring itself -- modeled here as two concentric cylinders
// (outer wall, inner bore wall) joined by flat annuli at the leading and
// trailing edges, since the schema states only a chord and an inner/outer
// diameter, not a shaped duct-lip profile.
//
// FRAME. Same 180-degree-about-y conversion as BuildBody(): contract x
// forward, solver x aft. The duct's leading edge is the LARGER contract x
// (further forward), which is the SMALLER solver x, so axial marching below
// runs leading -> trailing exactly as BuildBody() runs nose -> tail.
//
// WINDING. Each wall's outward normal faces away from the material: radially
// outward for the outer wall, radially INWARD (toward the axis) for the
// bore wall, forward for the leading cap, aft for the trailing cap --
// mirroring BuildBody()'s aft-facing base cap.
std::vector<Lattice::SourcePanel> LatticeBuilder::BuildDuct() const {
    std::vector<Lattice::SourcePanel> panels;
    if (!m_Options.IncludeDuct) return panels;

    const Geometry::DuctGeometry& duct = m_Contract.Duct;
    if (!duct.IsStated) return panels;

    const int sectors = m_Options.DuctCircumferentialPanels;
    if (sectors < MinDuctSectors) return panels;
    const int axialPanels = std::max(1, m_Options.DuctAxialPanels);

    const double rInner = duct.InnerDiameter * Math::Half;
    const double rOuter = duct.OuterDiameter * Math::Half;

    const double xLeadSolver = -(duct.Center.x + duct.Chord * Math::Half);
    const double xTrailSolver = -(duct.Center.x - duct.Chord * Math::Half);
    const double yOffset = duct.Center.y;
    const double zOffset = -duct.Center.z;

    const auto surfacePoint = [&](double x, double r, double angle) {
        return Solver::Vec3(x, yOffset + r * std::cos(angle), zOffset + r * std::sin(angle));
    };

    // Identical corner-order fixup to BuildBody()'s finishPanel.
    const auto finishPanel = [&](Lattice::SourcePanel& panel, const Solver::Vec3& outwardHint) {
        Solver::Vec3 normal = Solver::Cross(panel.Corners[2] - panel.Corners[0],
                                            panel.Corners[3] - panel.Corners[1]);
        if (Solver::Dot(normal, outwardHint) < 0.0) {
            std::swap(panel.Corners[1], panel.Corners[3]);
            normal = normal * -1.0;
        }
        panel.Area = Math::Half * Solver::Cross(panel.Corners[2] - panel.Corners[0],
                                                panel.Corners[3] - panel.Corners[1])
                                      .Norm();
        panel.Normal = normal.Normalized();
    };

    const auto centroidOf = [](const std::array<Solver::Vec3, 4>& corners) {
        Solver::Vec3 centroid(0, 0, 0);
        for (const Solver::Vec3& corner : corners) centroid = centroid + corner * 0.25;
        return centroid;
    };

    // --- cylindrical walls ---------------------------------------------------
    const auto emitCylinder = [&](double radius, bool facingOutward, const char* surfaceName) {
        for (int a = 0; a < axialPanels; ++a) {
            const double xNear = xLeadSolver + (xTrailSolver - xLeadSolver) * a / axialPanels;
            const double xFar = xLeadSolver + (xTrailSolver - xLeadSolver) * (a + 1) / axialPanels;

            for (int j = 0; j < sectors; ++j) {
                const double angle0 = Math::Two * std::numbers::pi * j / sectors;
                const double angle1 = Math::Two * std::numbers::pi * (j + 1) / sectors;

                Lattice::SourcePanel panel;
                panel.Corners = {surfacePoint(xNear, radius, angle0), surfacePoint(xFar, radius, angle0),
                                 surfacePoint(xFar, radius, angle1), surfacePoint(xNear, radius, angle1)};

                const double midAngle = Math::Half * (angle0 + angle1);
                const Solver::Vec3 radial(0.0, std::cos(midAngle), std::sin(midAngle));
                finishPanel(panel, facingOutward ? radial : radial * -1.0);

                panel.ControlPoint = centroidOf(panel.Corners);
                panel.Surface = surfaceName;
                panel.StationIndex = a;
                if (panel.Area > BodyDegenerateArea) panels.push_back(panel);
            }
        }
    };
    emitCylinder(rOuter, true, DuctOuterSurfaceName);
    emitCylinder(rInner, false, DuctInnerSurfaceName);

    // --- end caps: annuli from rInner to rOuter -----------------------------
    const auto emitCap = [&](double x, const Solver::Vec3& outward, const char* surfaceName,
                             int stationIndex) {
        for (int j = 0; j < sectors; ++j) {
            const double angle0 = Math::Two * std::numbers::pi * j / sectors;
            const double angle1 = Math::Two * std::numbers::pi * (j + 1) / sectors;

            Lattice::SourcePanel panel;
            panel.Corners = {surfacePoint(x, rInner, angle0), surfacePoint(x, rOuter, angle0),
                             surfacePoint(x, rOuter, angle1), surfacePoint(x, rInner, angle1)};
            finishPanel(panel, outward);

            panel.ControlPoint = centroidOf(panel.Corners);
            panel.Surface = surfaceName;
            panel.StationIndex = stationIndex;
            if (panel.Area > BodyDegenerateArea) panels.push_back(panel);
        }
    };
    emitCap(xLeadSolver, Solver::Vec3(-1.0, 0.0, 0.0), DuctLeadingCapSurfaceName, 0);
    emitCap(xTrailSolver, Solver::Vec3(1.0, 0.0, 0.0), DuctTrailingCapSurfaceName, axialPanels);

    return panels;
}


// --- base efflux ------------------------------------------------------------
int ApplyBaseEfflux(std::vector<Lattice::SourcePanel>& body,
                    const std::function<Math::Vec3(const Math::Vec3&)>& slipstream,
                    const Math::Vec3& referenceFreestream) {
    int affected = 0;
    for (Lattice::SourcePanel& panel : body) {
        // Only the cap transpires. The rest of the body is a wall and stays
        // one, so this cannot quietly perforate the whole fuselage.
        if (panel.Surface != BodyBaseSurfaceName) continue;

        const Math::Vec3 induced = slipstream ? slipstream(panel.ControlPoint) : Math::Vec3(0, 0, 0);
        panel.PrescribedNormalVelocity = Math::Dot(referenceFreestream + induced, panel.Normal);
        // The cap stops being a wall: it still imposes the efflux, but no
        // pressure force is integrated over it (Lattice::SourcePanel).
        panel.Permeable = true;
        ++affected;
    }
    return affected;
}

std::vector<Lattice::Panel> BuildPropellerLattice(const Geometry::Propeller& prop,
                                                  double axialSpeed, double omega) {
    std::vector<Lattice::Panel> panels;
    const std::size_t ns = prop.Stations.size();
    if (ns < 2 || prop.BladeCount < 1) return panels;

    panels.reserve(static_cast<std::size_t>(prop.BladeCount) * (ns - 1));

    const Math::Vec3 axial(1.0, 0.0, 0.0); // rotation axis, +x

    // Trailing-leg direction at a wake-root point: the local kinematic
    // velocity (axial inflow plus the -Omega x r sweep the rotation term
    // hands the blade), i.e. the linearized helix -- see the header note.
    //
    // The axial part is floored at the momentum-theory hover inflow,
    // lambda ~ 0.07 (v_i = lambda * Omega * R, the classic value for
    // ordinary propeller disk loadings, since C_T ~ 0.01 gives
    // lambda = sqrt(C_T/2) ~ 0.07). At exact hover the kinematic axial
    // inflow is zero and a straight leg would otherwise lie IN the rotor
    // plane forever, slicing past every strip outboard of its root (a
    // straight tangent line leaves its circle and crosses all larger
    // radii) -- a real wake convects out of the disk plane, and without
    // this floor the solve is pathologically sensitive to the station
    // layout.
    const double hoverInflow = HoverInflowRatio * std::fabs(omega) * prop.Radius;
    auto trailDirection = [&](const Math::Vec3& point) {
        const Math::Vec3 sweep = Math::Cross(axial * omega, point);
        const Math::Vec3 local = axial * std::max(axialSpeed, hoverInflow) - sweep;
        const Math::Vec3 unit = local.Normalized();
        return (unit.Norm() > 0.5) ? unit : axial; // no flow at all -> fall back to axial
    };

    for (int blade = 0; blade < prop.BladeCount; ++blade) {
        const double bladeAzimuth = Math::Two * std::numbers::pi * blade / prop.BladeCount;

        // Local chord frame at chord fraction f of a station, WRAPPED around
        // the cylinder of that station's radius. A blade element's chord
        // does not live in a flat tangent plane -- it lies on the annulus
        // the element sweeps, so moving along the chord changes azimuth by
        // s*cos(beta)/r. Near a small hub this is not a nicety: chord there
        // is comparable to radius, a straight tangent-plane chord subtends
        // tens of degrees of azimuth and leaves the swept surface entirely,
        // and the resulting phantom incidence overwhelms real camber/twist
        // loading (found empirically: the hub strips reversed the sign of
        // the whole blade's camber response).
        //
        // The blade moves toward +that, so the relative wind arrives from
        // +that: the leading edge faces +that (wrapped to positive azimuth
        // offset), and twist tilts the trailing edge aft (+x), which is
        // what makes the resultant blade force point upstream (-x) --
        // thrust -- at positive local incidence.
        //
        // CamberDir is the suction side the CST convention bows positive
        // camber toward: for a propeller that is the THRUST side, so a
        // positively-cambered section thrusts at zero twist the way a
        // wing's lifts at zero incidence (TestPropellerLattice pins the
        // sign).
        struct ChordFrame {
            Math::Vec3 Point;     // on the chord line, wrapped on the cylinder
            Math::Vec3 ChordDir;  // local LE -> TE tangent, unit
            Math::Vec3 CamberDir; // suction (thrust) side, unit
        };
        auto frameAt = [&](const Geometry::BladeStation& station, double f) {
            const double beta = Math::DegToRad(station.TwistDeg);
            const double s = (f - Math::Half) * station.Chord; // arc length from mid-chord, +aft
            const double radius = std::max(station.r, 1e-12);
            const double azimuth = bladeAzimuth - s * std::cos(beta) / radius;
            const Math::Vec3 rhat(0.0, std::cos(azimuth), std::sin(azimuth));
            const Math::Vec3 that(0.0, -std::sin(azimuth), std::cos(azimuth)); // axial x rhat
            ChordFrame frame;
            frame.ChordDir = -that * std::cos(beta) + axial * std::sin(beta);
            frame.CamberDir = -Math::Cross(frame.ChordDir, rhat); // unit: ChordDir and rhat are orthonormal
            frame.Point = rhat * station.r + axial * (s * std::sin(beta));
            return frame;
        };

        // A point on the CAMBER SURFACE at chord fraction f: the wrapped
        // chord-line point offset by the interpolated CST mean-line
        // ordinate toward the suction side, exactly the way the wing's
        // lattice sits on its camber surface. With no section data
        // CamberAt() is zero and this is the chord line itself.
        auto surfacePoint = [&](const Geometry::BladeStation& station, double f) {
            const ChordFrame frame = frameAt(station, f);
            const double eta = (prop.Radius > 0.0) ? station.r / prop.Radius : 0.0;
            const double camber = Geometry::CamberAt(prop.Sections, eta, f);
            return frame.Point + frame.CamberDir * (camber * station.Chord);
        };

        for (std::size_t i = 0; i + 1 < ns; ++i) {
            const Geometry::BladeStation& inboard = prop.Stations[i];
            const Geometry::BladeStation& outboard = prop.Stations[i + 1];
            const double width = outboard.r - inboard.r;
            if (width <= 0.0) continue; // degenerate/unordered station pair carries no strip
            const double etaMid = (prop.Radius > 0.0)
                                      ? Math::Half * (inboard.r + outboard.r) / prop.Radius
                                      : 0.0;

            {
                const double f0 = 0.0;
                const double f1 = 1.0;
                const double fBound = Math::QuarterChord;
                const double fControl = Math::ThreeQuarterChord;

                Lattice::Panel panel;
                panel.A = surfacePoint(inboard, fBound);
                panel.B = surfacePoint(outboard, fBound);
                panel.ControlPoint = (surfacePoint(inboard, fControl) + surfacePoint(outboard, fControl)) * Math::Half;
                panel.TrailDirA = trailDirection(panel.A);
                panel.TrailDirB = trailDirection(panel.B);

                // Control-point normal from the surface tangents there: the
                // radial tangent between the stations, and the chordwise
                // tangent tilted by the ANALYTIC camber slope (the wing's
                // convention -- the flow-tangency condition wants the mean
                // line's own slope at the control point, not a facet
                // average). Flat sections reduce this to the slice normal.
                const ChordFrame inFrame = frameAt(inboard, fControl);
                const ChordFrame outFrame = frameAt(outboard, fControl);
                const double slope = Geometry::CamberSlopeAt(prop.Sections, etaMid, fControl);
                const Math::Vec3 radialTangent = surfacePoint(outboard, fControl) - surfacePoint(inboard, fControl);
                const Math::Vec3 chordTangent =
                    ((inFrame.ChordDir + outFrame.ChordDir) + (inFrame.CamberDir + outFrame.CamberDir) * slope);
                panel.Normal = Math::Cross(radialTangent, chordTangent).Normalized();

                // Footprint from the CHORD-LINE quad (the camber-free
                // slice): that is the planform-style area coefficients
                // normalize by and the viewer reconstructs the drawn quad
                // from. The true curved surface is longer by the camber arc
                // ratio, exactly the wing's Area vs PlanformArea
                // distinction.
                const Math::Vec3 le0 = frameAt(inboard, f0).Point;
                const Math::Vec3 te0 = frameAt(inboard, f1).Point;
                const Math::Vec3 le1 = frameAt(outboard, f0).Point;
                const Math::Vec3 te1 = frameAt(outboard, f1).Point;
                // Planar quad: half the cross product of its diagonals.
                panel.PlanformArea = Math::Cross(te1 - le0, le1 - te0).Norm() * Math::Half;
                const double arcFraction = Geometry::CamberArcLengthFraction(prop.Sections, etaMid, f0, f1);
                panel.Area = panel.PlanformArea * (arcFraction / (f1 - f0));
                panel.SpanwiseWidth = width;
                panel.Surface = "blade" + std::to_string(blade);
                panel.StripIndex = static_cast<int>(blade * (ns - 1) + i);
                panels.push_back(panel);
            }
        }
    }
    return panels;
}

std::vector<Solver::StripSection> BuildPropellerStrips(const Geometry::Propeller& prop) {
    std::vector<Solver::StripSection> strips;
    const std::size_t ns = prop.Stations.size();
    if (ns < 2 || prop.BladeCount < 1) return strips;
    strips.reserve(static_cast<std::size_t>(prop.BladeCount) * (ns - 1));

    const Math::Vec3 axial(1.0, 0.0, 0.0);

    for (int blade = 0; blade < prop.BladeCount; ++blade) {
        const double bladeAzimuth = Math::Two * std::numbers::pi * blade / prop.BladeCount;

        // Mid-chord section frame of a station, wrapped on its radius
        // cylinder exactly the way BuildPropellerLattice wraps the panels
        // (mid-chord: azimuth offset zero, so the frame is the station's).
        const auto frameAt = [&](const Geometry::BladeStation& station) {
            const double beta = Math::DegToRad(station.TwistDeg);
            const Math::Vec3 rhat(0.0, std::cos(bladeAzimuth), std::sin(bladeAzimuth));
            const Math::Vec3 that(0.0, -std::sin(bladeAzimuth), std::cos(bladeAzimuth));
            const Math::Vec3 chordDir = -that * std::cos(beta) + axial * std::sin(beta);
            return std::pair{chordDir, -Math::Cross(chordDir, rhat)};
        };

        for (std::size_t i = 0; i + 1 < ns; ++i) {
            const Geometry::BladeStation& inboard = prop.Stations[i];
            const Geometry::BladeStation& outboard = prop.Stations[i + 1];
            const double width = outboard.r - inboard.r;
            if (width <= 0.0) continue; // must skip exactly what the lattice skips, to stay aligned

            const auto [chordIn, liftIn] = frameAt(inboard);
            const auto [chordOut, liftOut] = frameAt(outboard);
            const double etaMid = (prop.Radius > 0.0)
                                      ? Math::Half * (inboard.r + outboard.r) / prop.Radius
                                      : 0.0;

            Solver::StripSection strip;
            strip.ChordDir = (chordIn + chordOut).Normalized();
            strip.LiftDir = (liftIn + liftOut).Normalized();
            strip.Chord = Math::Half * (inboard.Chord + outboard.Chord);
            strip.Width = width;
            strip.Alpha0Deg = Geometry::SectionZeroLiftAngleDeg(prop.Sections, etaMid);
            strips.push_back(strip);
        }
    }
    return strips;
}

} // namespace Aeolion::PanelBuilder
