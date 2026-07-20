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
#include <utility>
#include <vector>

namespace Aeolion::PanelBuilder {

LatticeBuilder::LatticeBuilder(Geometry::HandoffContract contract, LatticeOptions options)
    : m_Contract(std::move(contract)), m_Options(options) {
    // Breakpoints and the spanwise march depend only on the contract and
    // the spacing choice, never on a commanded deflection, so they are
    // computed once here rather than per Build().
    m_BoundaryEtas = ComputeBoundaryEtas();
    m_SemiSpan = ComputeSemiSpan();
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
           frame.UpDir * (camber * station.Chord);
}

// `inner`/`outer` are two adjacent semi-span stations (inner has the smaller
// |y|); `sign` is +1 for the right semi-span, -1 for the mirrored left one.
void LatticeBuilder::EmitStrip(const SpanStation& inner, const SpanStation& outer, double sign,
                                      int stripIndex, std::vector<Solver::Panel>& panels) const {
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

        const Solver::Vec3 innerPoint = deflect(inner, SurfacePoint(inner, sign, psiBound), psiBound);
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
        panel.SpanwiseWidth = outer.y - inner.y;

        // Projection for coefficients; true curved surface for geometry.
        // Deflection is a rigid rotation, so it leaves arc length -- and
        // therefore Area -- untouched.
        const double arcFraction =
            Geometry::CamberArcLengthFraction(m_Contract.AirfoilSections, mid.Eta, psiStart, psiEnd);
        const double dihedralRad = Math::DegToRad(mid.DihedralDeg);
        panel.PlanformArea = (mid.Chord * (psiEnd - psiStart)) * panel.SpanwiseWidth;
        panel.Area = (mid.Chord * arcFraction) * (panel.SpanwiseWidth / std::cos(dihedralRad));

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
    for (std::size_t k = 0; k + 1 < m_SemiSpan.size(); ++k) // left semi-span
        EmitStrip(m_SemiSpan[k], m_SemiSpan[k + 1], -1.0, stripIndex++, panels);
    for (std::size_t k = 0; k + 1 < m_SemiSpan.size(); ++k) // right semi-span
        EmitStrip(m_SemiSpan[k], m_SemiSpan[k + 1], +1.0, stripIndex++, panels);

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

} // namespace Aeolion::PanelBuilder
