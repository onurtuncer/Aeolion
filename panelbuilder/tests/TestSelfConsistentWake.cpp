// TestSelfConsistentWake.cpp -- the Level-A wake: the trailing-leg helix
// pitch iterated to consistency with the SOLVED loading, radially varying,
// through the outer driver's RotorBuilder hook.
//
// What must hold: the outer iteration converges; the converged wake pitch
// is a fixed point (rebuilding the lattice from the final solved inflow
// reproduces the final mesh's trailing legs); the pitch VARIES radially
// with the loading rather than following one global constant; and a more
// heavily loaded blade earns a steeper wake -- the load-dependence the
// hardcoded lambda could not have.
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/RotorVaneCoupling.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using namespace Aeolion;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

constexpr double Rpm = 6000.0;
constexpr double Rho = 1.225;
constexpr double HoverFloorSpeed = 0.01;
constexpr double Omega = Rpm * 2.0 * std::numbers::pi / 60.0;

Geometry::Propeller MakeTestProp(double chordScale) {
    Geometry::Propeller prop;
    prop.BladeCount = 2;
    prop.Radius = 0.15;
    prop.HubRadius = 0.02;
    const int n = 12;
    for (int i = 0; i < n; ++i) {
        const double r = prop.HubRadius + (prop.Radius - prop.HubRadius) * (i + 0.5) / (n + 1);
        const double eta = r / prop.Radius;
        prop.Stations.push_back({r, chordScale * 0.025 * (1.0 - 0.3 * eta), 30.0 - 20.0 * eta});
    }
    return prop;
}

Solver::RotorVaneResult SolveWake(const Geometry::Propeller& prop,
                                  int wakeRevolutions = PanelBuilder::DefaultPropellerWakeRevolutions) {
    const auto strips = PanelBuilder::BuildPropellerStrips(prop);
    const auto seed = PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega, {}, wakeRevolutions);

    Solver::FreestreamConditions fc;
    fc.Vinf = HoverFloorSpeed;
    fc.alphaDeg = 0.0;
    fc.betaDeg = 0.0;
    fc.rho = Rho;
    fc.p = Omega;
    fc.RefPoint = Solver::Vec3(0.0, 0.0, 0.0);

    Solver::ReferenceGeometry ref;
    for (const auto& panel : seed) ref.Area += panel.PlanformArea;
    ref.Span = 2.0 * prop.Radius;
    ref.Chord = ref.Area / ref.Span;
    const double trail = Solver::DefaultTrailSpanFactor * 2.0 * prop.Radius;

    const Solver::RotorBuilder rotorBuilder = [&](const std::function<double(double)>& inflow) {
        return PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega, inflow, wakeRevolutions);
    };

    return Solver::SolveRotorVaneCoupled(seed, strips, fc, ref, trail,
                                         Solver::AnalyticSectionModel{}, {}, nullptr, fc, ref,
                                         trail, Solver::AnalyticSectionModel{}, 0.0, Rho, {},
                                         rotorBuilder);
}

// The wake helix pitch of a panel's A-leg at DEPARTURE: the axial
// (out-of-disk-plane) component of the leg's first direction -- the first
// helix segment when the curved near wake is present.
double AxialPitch(const Solver::Panel& panel) {
    const Solver::Vec3 direction = panel.TrailPathA.empty()
                                       ? panel.TrailDirA
                                       : (panel.TrailPathA.front() - panel.A).Normalized();
    return direction.x;
}

void TestConvergesToAFixedPoint() {
    const auto prop = MakeTestProp(1.0);
    const auto wake = SolveWake(prop);
    std::cout << "self-consistent wake: outer=" << wake.OuterIterations
              << " residual=" << wake.Residual << "  T=" << -wake.Rotor.Base.Di
              << " N  Q=" << -(wake.Rotor.InducedMoment.x + wake.Rotor.ProfileMoment.x) << " N*m\n";
    CHECK(wake.Converged, "the wake iteration must converge");
    CHECK(wake.OuterIterations >= 2, "self-consistency needs at least a second pass");

    // Fixed point: rebuilding the lattice from the FINAL solved inflow must
    // reproduce the final mesh's trailing legs.
    const auto inflow = Solver::AxialInflowFromBands(
        Solver::ComputeSlipstreamBands(wake.Rotor, 0.0, Rho));
    const auto rebuilt = PanelBuilder::BuildPropellerLattice(
        prop, 0.0, Omega, inflow, PanelBuilder::DefaultPropellerWakeRevolutions);
    CHECK(rebuilt.size() == wake.BladePanels.size(), "the rebuild must match the mesh");
    double worst = 0.0;
    for (std::size_t i = 0; i < rebuilt.size(); ++i)
        worst = std::max(worst,
                         (rebuilt[i].TrailDirA - wake.BladePanels[i].TrailDirA).Norm());
    std::cout << "  fixed-point trailing-leg mismatch: " << worst << "\n";
    CHECK(worst < 0.05, "the converged wake pitch must be a fixed point of the loading");

    // The pitch must VARY radially with the loading -- the point of vi(r)
    // over one global lambda.
    double pitchMin = 1.0, pitchMax = 0.0;
    for (const auto& panel : wake.BladePanels) {
        pitchMin = std::min(pitchMin, AxialPitch(panel));
        pitchMax = std::max(pitchMax, AxialPitch(panel));
    }
    std::cout << "  axial pitch range across the blade: " << pitchMin << " .. " << pitchMax << "\n";
    CHECK(pitchMax / std::max(pitchMin, 1e-9) > 1.2,
          "the wake pitch must vary radially with the loading");
}

void TestHeavierLoadingSteepensTheWake() {
    const auto light = SolveWake(MakeTestProp(0.6));
    const auto heavy = SolveWake(MakeTestProp(1.4));
    CHECK(light.Converged && heavy.Converged, "both loadings must converge");

    const auto meanPitch = [](const Solver::RotorVaneResult& wake) {
        double sum = 0.0;
        for (const auto& panel : wake.BladePanels) sum += AxialPitch(panel);
        return sum / static_cast<double>(wake.BladePanels.size());
    };
    const double lightPitch = meanPitch(light);
    const double heavyPitch = meanPitch(heavy);
    std::cout << "mean axial pitch: light=" << lightPitch << "  heavy=" << heavyPitch
              << "  (T light=" << -light.Rotor.Base.Di << " N, heavy=" << -heavy.Rotor.Base.Di
              << " N)\n";
    CHECK(-heavy.Rotor.Base.Di > -light.Rotor.Base.Di, "the heavier blade must thrust more");
    CHECK(heavyPitch > lightPitch * 1.05,
          "a heavier disk loading must steepen the wake helix -- the load dependence the fixed "
          "lambda could not express");
}

void TestPolylineLegMatchesStraight() {
    // Kernel consistency: a curved-path leg whose points lie exactly on
    // the straight direction must reproduce the single-segment leg's
    // induction (zero core, so the regularization cannot differ).
    Solver::Panel straight;
    straight.A = Solver::Vec3(0.0, 0.0, 0.0);
    straight.B = Solver::Vec3(0.0, 1.0, 0.0);
    straight.TrailDirA = Solver::Vec3(1.0, 0.0, 0.0);
    straight.TrailDirB = Solver::Vec3(1.0, 0.0, 0.0);

    Solver::Panel polyline = straight;
    polyline.TrailPathA = {Solver::Vec3(0.3, 0.0, 0.0), Solver::Vec3(0.7, 0.0, 0.0)};
    polyline.TrailPathB = {Solver::Vec3(0.3, 1.0, 0.0), Solver::Vec3(0.7, 1.0, 0.0)};

    // Identical total extent: the polyline's far tail starts 0.7 further
    // along, so the straight comparison leg gets that much more trail.
    const Solver::Vec3 probe(0.4, 0.5, 0.8);
    const Solver::Vec3 a = Solver::HorseshoeVelocity(probe, straight, 1.0, 50.7);
    const Solver::Vec3 b = Solver::HorseshoeVelocity(probe, polyline, 1.0, 50.0);
    CHECK((a - b).Norm() < 1e-9 * a.Norm() + 1e-15,
          "a collinear polyline leg must reproduce the straight leg's induction");
}

void TestCurvedWakeAddsInducedPower() {
    // The point of the curved wake: straight legs cannot build the
    // contracted downstream cylinder that generates real induced power, so
    // the straight-leg hover figure of merit reads optimistic. The helical
    // wake must pull it DOWN -- more induced torque for comparable thrust.
    const auto prop = MakeTestProp(1.0);
    const auto straight = SolveWake(prop, 0);
    const auto curved = SolveWake(prop, PanelBuilder::DefaultPropellerWakeRevolutions);
    CHECK(straight.Converged && curved.Converged, "both wake models must converge");

    const auto figureOfMerit = [&](const Solver::RotorVaneResult& wake) {
        const double thrust = -wake.Rotor.Base.Di;
        const double torque = -(wake.Rotor.InducedMoment.x + wake.Rotor.ProfileMoment.x);
        const double power = torque * Omega;
        const double diskArea = std::numbers::pi * prop.Radius * prop.Radius;
        return thrust * std::sqrt(thrust / (2.0 * Rho * diskArea)) / std::max(power, 1e-12);
    };
    const double fomStraight = figureOfMerit(straight);
    const double fomCurved = figureOfMerit(curved);
    std::cout << "figure of merit: straight legs=" << fomStraight << "  helical wake=" << fomCurved
              << "  (T " << -straight.Rotor.Base.Di << " -> " << -curved.Rotor.Base.Di << " N)\n";
    CHECK(fomCurved < fomStraight,
          "the helical contracted wake must add induced power the straight legs cannot carry");
    CHECK(-curved.Rotor.Base.Di > 0.0, "and the rotor must still thrust");
}

} // namespace

int main() {
    TestConvergesToAFixedPoint();
    TestHeavierLoadingSteepensTheWake();
    TestPolylineLegMatchesStraight();
    TestCurvedWakeAddsInducedPower();

    if (failures == 0) { std::cout << "PASS: TestSelfConsistentWake\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestSelfConsistentWake\n";
    return 1;
}
