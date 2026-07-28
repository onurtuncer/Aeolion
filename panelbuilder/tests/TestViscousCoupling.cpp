// TestViscousCoupling.cpp -- the Level-2 sectional lift-feedback coupling
// (Solver::SolveViscousCoupled) on the propeller lattice.
//
// The checks pin what the coupling must DO, not what it happens to
// compute: the fixed point must converge with the residual driven under
// tolerance; stall saturation must cut hover thrust below the inviscid
// lattice's (whose root strips sit far past any real section's stall);
// the section drag must add profile torque the inviscid solve cannot
// have; zeroing the drag polar must zero exactly that contribution; the
// CST camber must surface in the strips as a negative zero-lift angle;
// and the physical advance trend must survive the coupling.
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

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
constexpr double HoverFloorSpeed = 0.01;
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

Solver::FreestreamConditions Conditions(double axialSpeed) {
    Solver::FreestreamConditions fc;
    fc.Vinf = std::max(axialSpeed, HoverFloorSpeed);
    fc.alphaDeg = 0.0;
    fc.betaDeg = 0.0;
    fc.rho = Rho;
    fc.p = Omega;
    fc.RefPoint = Solver::Vec3(0.0, 0.0, 0.0);
    return fc;
}

Solver::ViscousCoupledResult SolveCoupled(const Geometry::Propeller& prop, double axialSpeed,
                                          const Solver::SectionModel& model,
                                          const Solver::ViscousCouplingOptions& options = {}) {
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, axialSpeed, Omega);
    const auto strips = PanelBuilder::BuildPropellerStrips(prop);

    Solver::ReferenceGeometry ref;
    for (const auto& panel : panels) ref.Area += panel.PlanformArea;
    ref.Span = 2.0 * prop.Radius;
    ref.Chord = ref.Area / ref.Span;

    return Solver::SolveViscousCoupled(panels, strips, Conditions(axialSpeed), ref,
                                       Solver::DefaultTrailSpanFactor * 2.0 * prop.Radius, model,
                                       options);
}

void TestConvergesAndBoundsThrust() {
    const auto prop = MakeTestProp();
    const auto coupled = SolveCoupled(prop, 0.0, Solver::AnalyticSectionModel{});

    std::cout << "coupled hover: T=" << -coupled.Base.Di << " N  Q=" << -coupled.Base.Mx
              << " N*m  (induced " << -coupled.InducedMoment.x << " + profile "
              << -coupled.ProfileMoment.x << ")  iters=" << coupled.Iterations
              << "  residual=" << coupled.MaxResidual << "\n";

    CHECK(coupled.Converged, "the relaxed fixed point must converge at hover");
    CHECK(coupled.MaxResidual < Solver::DefaultCouplingTolerance,
          "the cl residual must be driven under tolerance");

    const double coupledThrust = -coupled.Base.Di;
    CHECK(coupledThrust > 0.0, "coupled hover thrust must stay positive");

    // The inviscid lattice runs the root strips at 20-30 deg incidence with
    // an unbounded linear lift slope; a real section stalls. Saturation
    // must therefore cut hover thrust, substantially.
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega);
    Solver::ReferenceGeometry ref;
    for (const auto& panel : panels) ref.Area += panel.PlanformArea;
    ref.Span = 2.0 * prop.Radius;
    ref.Chord = ref.Area / ref.Span;
    const double inviscidThrust =
        -Solver::Solve(panels, Conditions(0.0), ref, Solver::DefaultTrailSpanFactor * 2.0 * prop.Radius).Di;
    CHECK(coupledThrust < inviscidThrust,
          "stall saturation must cut hover thrust below the inviscid lattice's");

    // Every strip's converged cl must respect the model's own bound.
    for (const auto& strip : coupled.Strips)
        CHECK(std::fabs(strip.cl) <= Solver::AnalyticSectionModel{}.ClMax + 1e-9,
              "no strip may exceed the section model's ClMax");
}

void TestProfileTorqueIsRealAndSwitchable() {
    const auto prop = MakeTestProp();

    const auto viscous = SolveCoupled(prop, 0.0, Solver::AnalyticSectionModel{});
    const double profileTorque = -viscous.ProfileMoment.x;
    CHECK(profileTorque > 0.0, "section drag must add profile torque opposing the rotation");

    // Zero the drag polar: the profile contribution -- and only it -- must
    // vanish. (cl unchanged, so the induced part stays comparable.)
    Solver::AnalyticSectionModel dragless;
    dragless.Cd0 = 0.0;
    dragless.KCd = 0.0;
    const auto inviscidDrag = SolveCoupled(prop, 0.0, dragless);
    CHECK(std::fabs(-inviscidDrag.ProfileMoment.x) < 1e-12,
          "a zero drag polar must produce exactly zero profile torque");
    CHECK(-inviscidDrag.InducedMoment.x > 0.0, "the induced torque must remain");
}

void TestCamberSurfacesInStrips() {
    // Positive CST camber -> negative thin-airfoil zero-lift angle, carried
    // per strip so the section model represents the camber the coupling
    // supersedes in the lattice's boundary condition.
    auto prop = MakeTestProp();
    auto flat = PanelBuilder::BuildPropellerStrips(prop);
    for (const auto& strip : flat)
        CHECK(strip.Alpha0Deg == 0.0, "an uncambered blade's strips carry zero alpha0");

    prop.Sections.push_back({0.0, {0.05, 0.06, 0.05}, {0.05, 0.06, 0.05}});
    prop.Sections.push_back({1.0, {0.05, 0.06, 0.05}, {0.05, 0.06, 0.05}});
    auto cambered = PanelBuilder::BuildPropellerStrips(prop);
    CHECK(cambered.size() == flat.size(), "strips must align with the lattice regardless of camber");
    for (const auto& strip : cambered)
        CHECK(strip.Alpha0Deg < -0.5,
              "positive camber must give a clearly negative zero-lift angle");
}

void TestBoundaryLayerModelCouples() {
    // The Level-3 section model (transpiration-coupled boundary layer)
    // must drive the same outer fixed point to convergence and keep the
    // load structure physical. Its thrust will differ from the analytic
    // polar's -- different section physics -- but the signs and the
    // profile torque must not.
    auto prop = MakeTestProp();
    prop.Sections.push_back({0.0, {0.05, 0.06, 0.05}, {0.05, 0.06, 0.05}});
    prop.Sections.push_back({1.0, {0.05, 0.06, 0.05}, {0.05, 0.06, 0.05}});

    // The BL model's interior is piecewise (the transition station is
    // discrete), so its cl carries jumps a fixed point cannot iterate
    // below; the tolerance acknowledges that (see theory.rst).
    Solver::ViscousCouplingOptions options;
    options.Tolerance = 2e-2;
    options.MaxIterations = 400;
    const auto coupled = SolveCoupled(prop, 0.0, PanelBuilder::MakePropellerSectionModel(prop), options);
    std::cout << "BL-coupled hover: T=" << -coupled.Base.Di << " N  Q=" << -coupled.Base.Mx
              << " N*m  (profile " << -coupled.ProfileMoment.x << ")  iters=" << coupled.Iterations
              << "  residual=" << coupled.MaxResidual << "\n";
    for (std::size_t i = 0; i < coupled.Strips.size() / 2; ++i)
        std::cout << "  strip " << i << ": alphaEff=" << coupled.Strips[i].alphaEffDeg
                  << "  cl=" << coupled.Strips[i].cl << "  R=" << coupled.Strips[i].Residual
                  << "  Re=" << coupled.Strips[i].Re << "\n";

    CHECK(coupled.Converged, "the outer fixed point must converge with the BL section model");
    CHECK(-coupled.Base.Di > 0.0, "BL-coupled hover thrust must stay positive");
    CHECK(-coupled.Base.Mx > 0.0, "BL-coupled torque must oppose the rotation");
    CHECK(-coupled.ProfileMoment.x > 0.0, "the BL's Squire-Young drag must add profile torque");
}

void TestThrustStillUnloadsWithAdvance() {
    const auto prop = MakeTestProp();
    double previous = 1e30;
    for (double speed : {0.0, 5.0, 10.0}) {
        const auto coupled = SolveCoupled(prop, speed, Solver::AnalyticSectionModel{});
        const double thrust = -coupled.Base.Di;
        std::cout << "V=" << speed << "  T=" << thrust << " N  (iters=" << coupled.Iterations << ")\n";
        CHECK(coupled.Converged, "the coupling must converge at V=" + std::to_string(speed));
        CHECK(thrust < previous,
              "coupled fixed-pitch thrust must fall as inflow rises, V=" + std::to_string(speed));
        previous = thrust;
    }
}

} // namespace

int main() {
    TestConvergesAndBoundsThrust();
    TestProfileTorqueIsRealAndSwitchable();
    TestCamberSurfacesInStrips();
    TestBoundaryLayerModelCouples();
    TestThrustStillUnloadsWithAdvance();

    if (failures == 0) { std::cout << "PASS: TestViscousCoupling\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestViscousCoupling\n";
    return 1;
}
