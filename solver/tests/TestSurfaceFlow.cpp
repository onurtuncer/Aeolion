// TestSurfaceFlow.cpp -- validates the body skin-flow analysis against
// answers that are known in closed form, not against itself.
//
//   - Sphere, exact stagnation points. Potential flow past a sphere has
//     surface velocity V_t = (3/2)(U - (U.n)n) EXACTLY, so it vanishes
//     precisely where the outward normal is parallel to the freestream:
//     an attachment point at n = -U_hat and a separation point at n = +U_hat.
//     That holds at ANY angle of attack and sideslip, which makes the sphere
//     a complete test of the alpha/beta sweep rather than of one condition.
//
//   - Sphere, exact strain rate. Near the front stagnation point
//     |V_t| = (3/2) U s/a to first order in arc length s, so both principal
//     eigenvalues are (3/2)U/a and the Jacobian trace is 3U/a. This pins
//     the classification machinery quantitatively, and it is the same
//     dU_e/ds a boundary layer would start its march from.
//
//   - Sphere, exact streamline spreading. Surface streamlines are meridians
//     through the stagnation point, so the metric factor h(s) that the
//     3-D momentum-integral equation needs must grow like sin(theta): the
//     ratio h(s2)/h(s1) must equal sin(theta2)/sin(theta1) regardless of
//     where the trace was seeded.
//
//   - Zero incidence, the coordinate-singular case. At alpha = beta = 0 the
//     stagnation point of a body of revolution sits on the nose apex, which
//     is upstream of every control point and is NOT an interior zero. The
//     analysis must say so rather than invent one nearby.
//
//   - Prolate spheroid at combined incidence. The stagnation point's exact
//     position is not elementary, but its AZIMUTH is: any body of
//     revolution about x is symmetric about the plane containing its axis
//     and the crossflow, so the attachment point must lie on the windward
//     meridian, phi = atan2(-sin(alpha)cos(beta), -sin(beta)). This is the
//     assertion that would catch an alpha/beta sign error.
#include "Aeolion/Solver/SurfaceFlow.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using Aeolion::Lattice::SourcePanel;
using Aeolion::Math::Cross;
using Aeolion::Math::Dot;
using Aeolion::Math::Vec3;
namespace S = Aeolion::Solver;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

constexpr double Pi = std::numbers::pi;

// A body of revolution about x, panelled as a structured (station, sector)
// grid. Station 0 is the UPSTREAM end -- the solver's freestream travels
// toward +x, so the flow arrives at the -x end, and this matches the
// fuselage builder's nose-first station ordering.
//
// `radiusAt(t)` is the body radius at t in [0,1] running nose to tail;
// `axialAt(t)` its x position.
template <typename RadiusFn, typename AxialFn>
std::vector<SourcePanel> BuildBodyOfRevolution(RadiusFn radiusAt, AxialFn axialAt, int stations,
                                               int sectors) {
    const auto position = [&](double t, double phi) {
        const double r = radiusAt(t);
        return Vec3(axialAt(t), r * std::cos(phi), r * std::sin(phi));
    };

    std::vector<SourcePanel> panels;
    panels.reserve(static_cast<std::size_t>(stations) * static_cast<std::size_t>(sectors));

    for (int i = 0; i < stations; ++i) {
        const double t0 = static_cast<double>(i) / stations;
        const double t1 = static_cast<double>(i + 1) / stations;
        for (int j = 0; j < sectors; ++j) {
            const double phi0 = 2.0 * Pi * j / sectors;
            const double phi1 = 2.0 * Pi * (j + 1) / sectors;

            SourcePanel panel;
            panel.Corners = {position(t0, phi0), position(t1, phi0), position(t1, phi1),
                             position(t0, phi1)};

            Vec3 centroid(0, 0, 0);
            for (const Vec3& corner : panel.Corners) centroid = centroid + corner * 0.25;
            panel.ControlPoint = centroid;

            Vec3 normal = Cross(panel.Corners[2] - panel.Corners[0], panel.Corners[3] - panel.Corners[1]);
            // Outward means away from the axis at this panel's own azimuth.
            const double midPhi = 0.5 * (phi0 + phi1);
            const Vec3 outward(0.0, std::cos(midPhi), std::sin(midPhi));
            if (Dot(normal, outward) < 0.0) {
                std::swap(panel.Corners[1], panel.Corners[3]);
                normal = normal * -1.0;
            }
            panel.Area = 0.5 * Cross(panel.Corners[2] - panel.Corners[0],
                                     panel.Corners[3] - panel.Corners[1])
                                   .Norm();
            panel.Normal = normal.Normalized();
            panel.Surface = "body";
            panel.StationIndex = i;
            panel.SectorIndex = j;
            if (panel.Area > 1e-18) panels.push_back(panel);
        }
    }
    return panels;
}

std::vector<SourcePanel> BuildSphere(double radius, int stations, int sectors) {
    // theta measured from the UPSTREAM pole, so t = 0 is the nose.
    return BuildBodyOfRevolution([radius](double t) { return radius * std::sin(Pi * t); },
                                 [radius](double t) { return -radius * std::cos(Pi * t); }, stations,
                                 sectors);
}

std::vector<SourcePanel> BuildSpheroid(double semiAxial, double semiRadial, int stations, int sectors) {
    return BuildBodyOfRevolution([semiRadial](double t) { return semiRadial * std::sin(Pi * t); },
                                 [semiAxial](double t) { return -semiAxial * std::cos(Pi * t); }, stations,
                                 sectors);
}

// Solve a source-panelled body and hand back its flow field. The
// PreparedSystem must outlive the field, so the caller owns it.
S::FlowField SolveBody(const S::PreparedSystem& prepared, const S::FreestreamConditions& fc) {
    const S::SolveResult result = S::SolveWithSystem(prepared, fc, S::ReferenceGeometry{});
    return S::MakeFlowField(prepared, fc, result.gamma, result.sigma);
}

// ---------------------------------------------------------------- tests ----

void TestSphereAtZeroIncidence() {
    const double radius = 1.0;
    const S::PanelSystem system{{}, BuildSphere(radius, 24, 24)};
    const S::PreparedSystem prepared = S::Prepare(system, 0.0);

    S::FreestreamConditions fc;
    fc.Vinf = 10.0;
    fc.alphaDeg = 0.0;
    fc.betaDeg = 0.0;

    const S::FlowField field = SolveBody(prepared, fc);
    const S::SurfaceGrid grid = S::BuildSurfaceGrid(field, "body");
    CHECK(grid.Valid(), "a fully indexed body must produce a valid surface grid");

    const S::SurfaceFlowTopology topology = S::AnalyzeSurfaceFlow(grid);
    CHECK(topology.AttachesUpstream,
          "at zero incidence the stagnation point is the nose apex, upstream of every control point");
    CHECK(S::PrimaryAttachment(topology) == nullptr,
          "no INTERIOR attachment point should be invented at zero incidence");
    CHECK((topology.NoseAttachment - Vec3(-radius, 0, 0)).Norm() < 0.05 * radius,
          "the reported nose attachment must be the upstream pole");

    // The rear stagnation point IS interior to nothing either -- it sits on
    // the aft apex -- so the topology should be clean of spurious zeros.
    for (const S::CriticalPoint& point : topology.CriticalPoints)
        CHECK(false, "zero incidence should produce no interior critical point, got one at station " +
                         std::to_string(point.Station));
}

// The core test: the sphere's stagnation points are exact at any attitude.
void TestSphereStagnationAcrossIncidence() {
    const double radius = 1.0;
    const S::PanelSystem system{{}, BuildSphere(radius, 40, 40)};
    const S::PreparedSystem prepared = S::Prepare(system, 0.0);

    struct Case { double alphaDeg, betaDeg; };
    const Case cases[] = {{10.0, 0.0}, {-12.0, 0.0}, {0.0, 9.0}, {8.0, -6.0}, {15.0, 12.0}};

    for (const Case& c : cases) {
        S::FreestreamConditions fc;
        fc.Vinf = 10.0;
        fc.alphaDeg = c.alphaDeg;
        fc.betaDeg = c.betaDeg;

        const S::FlowField field = SolveBody(prepared, fc);
        const S::SurfaceGrid grid = S::BuildSurfaceGrid(field, "body");
        const S::SurfaceFlowTopology topology = S::AnalyzeSurfaceFlow(grid);

        const std::string label =
            "alpha=" + std::to_string(c.alphaDeg) + " beta=" + std::to_string(c.betaDeg);

        const Vec3 stream = S::FreestreamVelocity(fc).Normalized();
        const S::CriticalPoint* attachment = S::PrimaryAttachment(topology);
        const S::CriticalPoint* separation = S::PrimarySeparation(topology);

        CHECK(attachment != nullptr, "an attachment point must be found at " + label);
        CHECK(separation != nullptr, "a separation point must be found at " + label);
        if (!attachment || !separation) continue;

        // Compare DIRECTIONS: the panel control points sit fractionally
        // inside the sphere (they are quad centroids), so the radius carries
        // a faceting offset the direction does not.
        const Vec3 attachDir = attachment->Point.Normalized();
        const Vec3 separateDir = separation->Point.Normalized();
        const double attachError = (attachDir - stream * -1.0).Norm();
        const double separateError = (separateDir - stream).Norm();

        // One panel spans 180/32 deg; half a panel of slack.
        const double tolerance = 0.06;
        CHECK(attachError < tolerance,
              "attachment must sit at n = -Vinf_hat at " + label + ", direction error " +
                  std::to_string(attachError));
        CHECK(separateError < tolerance,
              "separation must sit at n = +Vinf_hat at " + label + ", direction error " +
                  std::to_string(separateError));

        CHECK(attachment->Type == S::CriticalPointType::AttachmentNode,
              "the forward stagnation point of a sphere is a NODE, not a focus or saddle, at " + label);
        CHECK(separation->Type == S::CriticalPointType::SeparationNode,
              "the aft stagnation point of a sphere is a node at " + label);

        // Exact strain rate: both eigenvalues are (3/2)U/a, so tr(J) = 3U/a.
        // This is the quantity a boundary layer would start its march from
        // (theta_0^2 = 0.075 nu / (dU_e/ds)), so its scale matters, not just
        // its sign. Panelling error here is about 4% at N = 48 and does not
        // fall much further, so 15% is honest slack at this resolution.
        const double expectedTrace = 3.0 * fc.Vinf / radius;
        const double traceError = std::fabs(attachment->Trace - expectedTrace) / expectedTrace;
        CHECK(traceError < 0.15, "attachment strain rate must be 3U/a at " + label + ", got " +
                                     std::to_string(attachment->Trace) + " vs " +
                                     std::to_string(expectedTrace));
        CHECK(attachment->Trace > 0.0 && separation->Trace < 0.0,
              "attachment must have positive trace and separation negative at " + label);

        // A sphere's stagnation point is a STAR node -- two equal
        // eigenvalues -- so the measured split is a weak, noisy quantity and
        // is asserted loosely on purpose. It is the ratio of two nearly
        // equal numbers, it scatters with where the zero happens to fall
        // inside a cell, and it converges only slowly (0.28, 0.67, 0.85,
        // 0.94 at alpha = 10 deg for N = 24, 32, 48, 64). The sharp
        // assertions in this test are the POSITION and the TRACE; this one
        // exists to catch a gross anisotropy, not to pin a number.
        if (attachment->Eigenvalues.size() == 2 && attachment->Eigenvalues[0] != 0.0) {
            const double ratio = attachment->Eigenvalues[1] / attachment->Eigenvalues[0];
            CHECK(ratio > 0.4, "a sphere's stagnation node is near-isotropic: eigenvalue ratio " +
                                   std::to_string(ratio) + " at " + label);
        }
    }
}

// Streamlines from the stagnation point are great circles; the spreading
// metric the 3-D boundary layer needs must follow sin(theta).
void TestSphereStreamlineSpreading() {
    const double radius = 1.0;
    const S::PanelSystem system{{}, BuildSphere(radius, 40, 40)};
    const S::PreparedSystem prepared = S::Prepare(system, 0.0);

    S::FreestreamConditions fc;
    fc.Vinf = 10.0;
    fc.alphaDeg = 12.0;
    fc.betaDeg = 0.0;

    const S::FlowField field = SolveBody(prepared, fc);
    const S::SurfaceGrid grid = S::BuildSurfaceGrid(field, "body");
    const S::SurfaceFlowTopology topology = S::AnalyzeSurfaceFlow(grid);
    const S::CriticalPoint* attachment = S::PrimaryAttachment(topology);
    CHECK(attachment != nullptr, "spreading test needs an attachment point");
    if (!attachment) return;

    const std::vector<S::SurfaceStreamline> lines = S::TraceFromCriticalPoint(grid, *attachment);
    CHECK(!lines.empty(), "an attachment node must emit distinguished streamlines");
    if (lines.empty()) return;

    const Vec3 stagnationDir = attachment->Point.Normalized();

    // Take the longest line -- the one that ran furthest before leaving the
    // patch -- and check its spreading against sin(theta).
    const S::SurfaceStreamline* longest = &lines.front();
    for (const S::SurfaceStreamline& line : lines)
        if (line.Length > longest->Length) longest = &line;

    CHECK(longest->Length > 0.5 * Pi * radius,
          "a streamline from the forward stagnation point should run most of a meridian, got " +
              std::to_string(longest->Length));

    // Polar angle of a point measured from the stagnation direction. On a
    // sphere the streamline IS a meridian through both stagnation points, so
    // theta parameterizes it exactly.
    const auto polarAngle = [&](const Vec3& p) {
        const double c = std::clamp(Dot(p.Normalized(), stagnationDir), -1.0, 1.0);
        return std::acos(c);
    };

    // Sample by POLAR ANGLE rather than by point index. Output points are
    // not uniformly distributed along the line -- the integrator slows into
    // the aft stagnation node -- so an index-based sample would compare two
    // points that are both effectively at theta = pi.
    const auto pointNearAngle = [&](double target) -> const S::StreamlinePoint* {
        const S::StreamlinePoint* best = nullptr;
        double bestError = 1e30;
        for (const S::StreamlinePoint& p : longest->Points) {
            const double error = std::fabs(polarAngle(p.Point) - target);
            if (error < bestError) { bestError = error; best = &p; }
        }
        return (bestError < Aeolion::Math::DegToRad(5.0)) ? best : nullptr;
    };

    const S::StreamlinePoint* a = pointNearAngle(Pi / 4.0);  // 45 deg
    const S::StreamlinePoint* b = pointNearAngle(Pi / 2.0);  // 90 deg, the equator
    CHECK(a != nullptr && b != nullptr,
          "the traced meridian must pass through 45 and 90 degrees of polar angle");
    if (!a || !b) return;

    // h(s) must grow like sin(theta). The RATIO is what is checked, because
    // it cancels the arbitrary normalization h = 1 at the seed.
    const double thetaA = polarAngle(a->Point), thetaB = polarAngle(b->Point);
    if (a->Divergence > 1e-12) {
        const double measured = b->Divergence / a->Divergence;
        const double expected = std::sin(thetaB) / std::sin(thetaA);
        const double error = std::fabs(measured - expected) / expected;
        CHECK(error < 0.15, "streamline spreading h(s) must follow sin(theta): measured ratio " +
                                std::to_string(measured) + " vs expected " + std::to_string(expected));
    }

    // And the edge speed along it must follow (3/2) U sin(theta).
    const double expectedSpeed = 1.5 * fc.Vinf * std::sin(thetaB);
    const double speedError = std::fabs(b->Speed - expectedSpeed) / expectedSpeed;
    CHECK(speedError < 0.10, "edge speed along the streamline must be (3/2)U sin(theta): got " +
                                 std::to_string(b->Speed) + " vs " + std::to_string(expectedSpeed));
}

// A fuselage-like body: the stagnation point's position is not elementary,
// but its azimuth is fixed exactly by symmetry.
void TestSpheroidWindwardMeridian() {
    // Fineness and incidence are chosen so the stagnation point is actually
    // RESOLVED. On a body of revolution it sits roughly one nose radius of
    // curvature times the crossflow angle back from the apex, and the first
    // ring of control points sits half a station aft of the apex -- so a
    // slender body at small incidence hides its stagnation point inside the
    // first panel, where the honest answer is AttachesUpstream and not a
    // located point. Here R_nose = b^2/a = 0.32 m and the crossflow angles
    // are 10 deg and up, which puts the stagnation point several control
    // points aft of the nose.
    const S::PanelSystem system{{}, BuildSpheroid(2.0, 0.8, 40, 32)};
    const S::PreparedSystem prepared = S::Prepare(system, 0.0);

    struct Case { double alphaDeg, betaDeg; };
    const Case cases[] = {{12.0, 0.0}, {0.0, 11.0}, {9.0, 9.0}, {-10.0, 7.0}};

    for (const Case& c : cases) {
        S::FreestreamConditions fc;
        fc.Vinf = 30.0;
        fc.alphaDeg = c.alphaDeg;
        fc.betaDeg = c.betaDeg;

        const S::FlowField field = SolveBody(prepared, fc);
        const S::SurfaceGrid grid = S::BuildSurfaceGrid(field, "body");
        const S::SurfaceFlowTopology topology = S::AnalyzeSurfaceFlow(grid);
        const S::CriticalPoint* attachment = S::PrimaryAttachment(topology);

        const std::string label =
            "alpha=" + std::to_string(c.alphaDeg) + " beta=" + std::to_string(c.betaDeg);
        CHECK(attachment != nullptr, "a spheroid at incidence must have an attachment point, " + label);
        if (!attachment) continue;

        // Windward meridian: opposite the crossflow, which is
        // (V_y, V_z) = V (sin beta, sin alpha cos beta).
        const double alpha = Aeolion::Math::DegToRad(c.alphaDeg);
        const double beta = Aeolion::Math::DegToRad(c.betaDeg);
        const double expectedPhi = std::atan2(-std::sin(alpha) * std::cos(beta), -std::sin(beta));
        const double actualPhi = std::atan2(attachment->Point.z, attachment->Point.y);

        double difference = actualPhi - expectedPhi;
        while (difference > Pi) difference -= 2.0 * Pi;
        while (difference < -Pi) difference += 2.0 * Pi;

        // One sector spans 360/32 = 11.25 deg; allow a sector of slack.
        CHECK(std::fabs(difference) < Aeolion::Math::DegToRad(12.0),
              "the attachment point must lie on the windward meridian at " + label + ": got phi=" +
                  std::to_string(Aeolion::Math::RadToDeg(actualPhi)) + " deg, expected " +
                  std::to_string(Aeolion::Math::RadToDeg(expectedPhi)) + " deg");

        // And it must sit on the forward part of the body.
        CHECK(attachment->Point.x < 0.0,
              "the attachment point belongs on the upstream half of the body at " + label);
    }
}

// An unindexed surface must be declined, not guessed at.
void TestUnstructuredSurfaceIsDeclined() {
    std::vector<SourcePanel> panels = BuildSphere(1.0, 12, 12);
    for (SourcePanel& panel : panels) panel.SectorIndex = -1; // topology withheld

    const S::PanelSystem system{{}, panels};
    const S::PreparedSystem prepared = S::Prepare(system, 0.0);
    S::FreestreamConditions fc;
    fc.alphaDeg = 5.0;

    const S::FlowField field = SolveBody(prepared, fc);
    const S::SurfaceGrid grid = S::BuildSurfaceGrid(field, "body");
    CHECK(!grid.Valid(), "a surface without a stated topology must be declined, not reconstructed");
    CHECK(S::AnalyzeSurfaceFlow(grid).CriticalPoints.empty(),
          "an invalid grid must yield no critical points");
}

} // namespace

int main() {
    TestSphereAtZeroIncidence();
    TestSphereStagnationAcrossIncidence();
    TestSphereStreamlineSpreading();
    TestSpheroidWindwardMeridian();
    TestUnstructuredSurfaceIsDeclined();

    if (failures == 0) std::cout << "PASS: TestSurfaceFlow\n";
    return failures == 0 ? 0 : 1;
}
