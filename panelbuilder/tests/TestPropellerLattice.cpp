// TestPropellerLattice.cpp -- the panel-method propeller: blades meshed by
// BuildPropellerLattice() and spun through Solver::Solve's body-rate term.
//
// The checks are physical, not structural: a lifting blade schedule spun at
// +Omega must thrust upstream (-x), absorb torque against its own rotation
// (-Mx), cancel its side/vertical forces by blade symmetry, and lose thrust
// as axial inflow rises at fixed rpm (a fixed-pitch prop unloading with
// advance ratio). A sign error anywhere in the blade orientation, the
// rotation sense, or the force bookkeeping fails these loudly.
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/Solver.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>

using namespace Aeolion;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

constexpr double Rpm = 6000.0;
constexpr double Rho = 1.225;
constexpr double HoverFloorSpeed = 0.01; // wind-axis force directions need Vinf > 0
constexpr double Omega = Rpm * 2.0 * std::numbers::pi / 60.0;

Geometry::Propeller MakeTestProp() {
    Geometry::Propeller prop;
    prop.BladeCount = 2;
    prop.Radius = 0.15;
    prop.HubRadius = 0.02;
    const int n = 12;
    for (int i = 0; i < n; ++i) {
        const double r = prop.HubRadius + (prop.Radius - prop.HubRadius) * (i + 0.5) / (n + 1);
        const double eta = r / prop.Radius;
        prop.Stations.push_back({r, 0.025 * (1.0 - 0.3 * eta), 30.0 - 20.0 * eta});
    }
    return prop;
}

Solver::SolveResult SolveProp(const Geometry::Propeller& prop, const std::vector<Solver::Panel>& panels,
                              double axialSpeed) {
    Solver::FreestreamConditions fc;
    fc.Vinf = std::max(axialSpeed, HoverFloorSpeed);
    fc.alphaDeg = 0.0; // the default is a wing-study 5 degrees -- axial inflow here
    fc.betaDeg = 0.0;
    fc.rho = Rho;
    fc.p = Rpm * 2.0 * std::numbers::pi / Math::SecondsPerMinute;
    fc.RefPoint = Solver::Vec3(0.0, 0.0, 0.0);

    Solver::ReferenceGeometry ref;
    for (const auto& panel : panels) ref.Area += panel.PlanformArea;
    ref.Span = 2.0 * prop.Radius;
    ref.Chord = ref.Area / ref.Span;

    return Solver::Solve(panels, fc, ref, Solver::DefaultTrailSpanFactor * 2.0 * prop.Radius);
}

void TestLatticeStructure() {
    const auto prop = MakeTestProp();
    const int rows = 4;
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, rows, 0.0, Omega);

    const std::size_t expected =
        static_cast<std::size_t>(prop.BladeCount) * (prop.Stations.size() - 1) * rows;
    CHECK(panels.size() == expected, "panel count should be blades * strips * chordwise rows");

    for (const auto& panel : panels) {
        CHECK(std::fabs(panel.Normal.Norm() - 1.0) < 1e-9, "every panel normal must be unit length");
        CHECK(panel.Area > 0.0, "every panel must have positive area");
        CHECK(panel.SpanwiseWidth > 0.0, "every panel must have positive radial width");
        const double rA = std::hypot(panel.A.y, panel.A.z);
        const double rB = std::hypot(panel.B.y, panel.B.z);
        CHECK(rA < rB, "A must be the inboard bound-vortex endpoint");
        CHECK(rB <= prop.Radius + 1e-12, "no panel may reach outboard of the tip");
        CHECK(std::fabs(panel.TrailDirA.Norm() - 1.0) < 1e-9, "trailing-leg directions must be unit");
        // At hover the local wind is almost purely tangential: the leg must
        // leave along it (opposing the blade's own +Omega sweep), not axially.
        const Solver::Vec3 sweepA = Solver::Cross(Solver::Vec3(Omega, 0.0, 0.0), panel.A);
        const Solver::Vec3 expectedA = (-sweepA).Normalized();
        CHECK(Solver::Dot(panel.TrailDirA, expectedA) > 0.999,
              "each trailing leg must follow the local relative wind");
    }
}

void TestHoverForcesArePhysical() {
    const auto prop = MakeTestProp();
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, 4, 0.0, Omega);
    const auto result = SolveProp(prop, panels, 0.0);

    const double thrust = -result.Di;
    const double torque = -result.Mx;
    std::cout << "hover: T=" << thrust << " N  Q=" << torque << " N*m  L=" << result.L
              << " Y=" << result.Y << " My=" << result.My << " Mz=" << result.Mz << "\n";
    CHECK(std::isfinite(thrust) && std::isfinite(torque), "hover solve must be finite");
    CHECK(thrust > 0.0, "a positively-twisted blade schedule at +Omega must thrust upstream (-x)");
    CHECK(torque > 0.0, "the prop must absorb torque against its own rotation (Mx < 0)");

    const double omega = Rpm * 2.0 * std::numbers::pi / Math::SecondsPerMinute;
    CHECK(torque * omega > 0.0, "absorbed power must be positive");

    // Two blades 180 degrees apart: the in-plane force components must
    // cancel by symmetry, leaving pure thrust + torque.
    CHECK(std::fabs(result.L) < 0.02 * thrust, "blade symmetry must cancel the vertical force");
    CHECK(std::fabs(result.Y) < 0.02 * thrust, "blade symmetry must cancel the side force");
}

void TestThrustUnloadsWithAdvance() {
    const auto prop = MakeTestProp();
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, 4, 0.0, Omega);

    double previous = 1e30;
    for (double speed : {0.0, 5.0, 10.0}) {
        const auto result = SolveProp(prop, panels, speed);
        const double thrust = -result.Di;
        std::cout << "V=" << speed << "  T=" << thrust << " N\n";
        CHECK(thrust < previous,
              "fixed-pitch thrust must fall as axial inflow rises, V=" + std::to_string(speed));
        previous = thrust;
    }
}

} // namespace

int main() {
    TestLatticeStructure();
    TestHoverForcesArePhysical();
    TestThrustUnloadsWithAdvance();

    if (failures == 0) { std::cout << "PASS: TestPropellerLattice\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestPropellerLattice\n";
    return 1;
}
