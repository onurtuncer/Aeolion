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
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega);

    const std::size_t expected =
        static_cast<std::size_t>(prop.BladeCount) * (prop.Stations.size() - 1);
    CHECK(panels.size() == expected, "panel count should be blades * strips (one Weissinger row each)");

    for (const auto& panel : panels) {
        CHECK(std::fabs(panel.Normal.Norm() - 1.0) < 1e-9, "every panel normal must be unit length");
        CHECK(panel.Area > 0.0, "every panel must have positive area");
        CHECK(panel.SpanwiseWidth > 0.0, "every panel must have positive radial width");
        const double rA = std::hypot(panel.A.y, panel.A.z);
        const double rB = std::hypot(panel.B.y, panel.B.z);
        CHECK(rA < rB, "A must be the inboard bound-vortex endpoint");
        CHECK(rB <= prop.Radius + 1e-12, "no panel may reach outboard of the tip");
        CHECK(std::fabs(panel.TrailDirA.Norm() - 1.0) < 1e-9, "trailing-leg directions must be unit");
        // At hover the local wind is dominantly tangential: the leg must
        // leave along it (opposing the blade's own +Omega sweep), with the
        // prescribed momentum-theory inflow as its axial part -- never
        // purely axial, and never lying in the disk plane.
        const Solver::Vec3 sweepA = Solver::Cross(Solver::Vec3(Omega, 0.0, 0.0), panel.A);
        const Solver::Vec3 expectedA =
            (Solver::Vec3(0.07 * Omega * prop.Radius, 0.0, 0.0) - sweepA).Normalized();
        CHECK(Solver::Dot(panel.TrailDirA, expectedA) > 0.999,
              "each trailing leg must follow the local relative wind plus prescribed inflow");
        CHECK(panel.TrailDirA.x > 0.0, "the wake must convect out of the disk plane (+x)");
    }
}

void TestHoverForcesArePhysical() {
    const auto prop = MakeTestProp();
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega);
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

void TestCamberThrustsAtZeroTwist() {
    // The blade analogue of the wing's "camber lifts at zero incidence"
    // check, which also pins the camber SIGN: positive CST camber bows
    // toward the suction side, and for a propeller the suction side is the
    // thrust side. An untwisted flat blade at hover sees ~zero incidence
    // and produces ~nothing; the same blade with positively-cambered
    // sections must thrust. A flipped camber direction would make it
    // thrust BACKWARD and fail loudly.
    // Posed OUTBOARD of the test prop's tiny hub deliberately: the sign pin
    // wants a well-conditioned blade-element region (chord well under the
    // local radius). The hub region has its own honest weirdness -- the
    // chord there wraps a large azimuth arc -- and that belongs to the
    // hover-magnitude caveats, not to a sign test.
    auto prop = MakeTestProp();
    prop.HubRadius = 0.07;
    std::erase_if(prop.Stations,
                  [&](const Geometry::BladeStation& st) { return st.r < prop.HubRadius; });
    for (auto& st : prop.Stations) st.TwistDeg = 0.0;

    const auto flatPanels = PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega);
    const double flatThrust = -SolveProp(prop, flatPanels, 0.0).Di;

    // Pure camber, no thickness: equal upper/lower coefficient sets give a
    // mean line that IS the surface. Modest ~2% camber.
    prop.Sections.push_back({0.0, {0.05, 0.06, 0.05}, {0.05, 0.06, 0.05}});
    prop.Sections.push_back({1.0, {0.05, 0.06, 0.05}, {0.05, 0.06, 0.05}});
    const auto camberedPanels = PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega);
    const double camberedThrust = -SolveProp(prop, camberedPanels, 0.0).Di;

    std::cout << "zero twist: flat T=" << flatThrust << " N  cambered T=" << camberedThrust << " N\n";
    {
        const auto r = SolveProp(prop, camberedPanels, 0.0);
        std::cout << "  cambered gamma:";
        for (double g : r.gamma) std::cout << " " << g;
        std::cout << "\n";
    }
    CHECK(std::fabs(flatThrust) < 1.0, "an untwisted flat blade at hover should produce ~no thrust");
    CHECK(camberedThrust > flatThrust + 1.0,
          "positively-cambered sections must thrust at zero twist -- this pins the camber sign");

    // Camber also lengthens the true wetted surface relative to the flat
    // footprint, the wing's Area vs PlanformArea distinction.
    CHECK(camberedPanels.front().Area > camberedPanels.front().PlanformArea,
          "a cambered panel's surface area must exceed its flat footprint");
}

void TestThrustUnloadsWithAdvance() {
    const auto prop = MakeTestProp();
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega);

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
    TestCamberThrustsAtZeroTwist();
    TestThrustUnloadsWithAdvance();

    if (failures == 0) { std::cout << "PASS: TestPropellerLattice\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestPropellerLattice\n";
    return 1;
}
