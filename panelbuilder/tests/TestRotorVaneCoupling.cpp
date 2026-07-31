// TestRotorVaneCoupling.cpp -- the two-way rotor-vane coupling
// (Solver::SolveRotorVaneCoupled): the block Gauss-Seidel alternation of
// the rotating-frame rotor(+duct) solve and the static-frame vane solve,
// with the swirl momentum budget.
//
// What must hold, and did not under the one-way scheme: the vanes cannot
// recover more torque than the jet's angular-momentum flux -- the shaft
// torque -- delivers (the one-way model was measured recovering three
// times that). And what makes it two-way at all: the rotor's loading must
// shift measurably when the vanes are present, because their blockage
// and upwash now enter its boundary condition.
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/RotorVaneCoupling.h"
#include "Aeolion/Solver/SectionBoundaryLayer.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

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

std::vector<Geometry::ControlSurface> MakeVaneSet() {
    std::vector<Geometry::ControlSurface> surfaces;
    const Math::Vec3 axes[] = {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, -1.0}};
    for (const Math::Vec3& axis : axes) {
        Geometry::ControlSurface vane;
        vane.Name = "vane";
        vane.Binding = Geometry::ControlSurfaceBinding::DuctJet;
        vane.ChordFraction = 1.0;
        vane.EtaStart = 0.4;
        vane.EtaEnd = 1.0;
        vane.HingeAxis = axis;
        surfaces.push_back(vane);
    }
    return surfaces;
}

Solver::RotorVaneResult SolveJoint(const std::vector<double>& deflectionsDeg) {
    const auto prop = MakeTestProp();
    const double shroudInner = 1.04 * prop.Radius;
    const double shroudChord = 0.5 * prop.Radius;
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega);
    const auto strips = PanelBuilder::BuildPropellerStrips(prop);
    const auto duct = PanelBuilder::BuildPropellerDuct(shroudInner, 1.2 * prop.Radius, shroudChord);
    const auto surfaces = MakeVaneSet();

    Solver::FreestreamConditions fc;
    fc.Vinf = HoverFloorSpeed;
    fc.alphaDeg = 0.0;
    fc.betaDeg = 0.0;
    fc.rho = Rho;
    fc.p = Omega;
    fc.RefPoint = Solver::Vec3(0.0, 0.0, 0.0);

    Solver::ReferenceGeometry ref;
    for (const auto& panel : panels) ref.Area += panel.PlanformArea;
    ref.Span = 2.0 * prop.Radius;
    ref.Chord = ref.Area / ref.Span;
    const double trail = Solver::DefaultTrailSpanFactor * 2.0 * prop.Radius;

    const auto vaneBuilder = [&](const std::function<Solver::Vec3(const Solver::Vec3&)>& flow)
        -> std::pair<std::vector<Solver::Panel>, std::vector<Solver::StripSection>> {
        return {PanelBuilder::BuildDuctVanes(surfaces, shroudInner, 0.5 * shroudChord, shroudChord,
                                             deflectionsDeg, flow),
                PanelBuilder::BuildDuctVaneStrips(surfaces, shroudInner, shroudChord, deflectionsDeg)};
    };

    Solver::FreestreamConditions vaneFc = fc;
    vaneFc.p = 0.0;
    Solver::ReferenceGeometry vaneRef;
    vaneRef.Area = 1.0;
    vaneRef.Span = 2.0 * shroudInner;
    vaneRef.Chord = shroudChord;

    // Both sides on the analytic polar: the helical-wake jet's swirl
    // angles sit beyond the BL section model's convergent envelope (it
    // hands such strips to this same analytic polar anyway), and the
    // budget iteration needs the smooth model to close. BL-vane coverage
    // lives in TestPropellerVanes' mild-jet case, inside the envelope.
    return Solver::SolveRotorVaneCoupled(panels, strips, fc, ref, trail,
                                         Solver::AnalyticSectionModel{}, duct, vaneBuilder, vaneFc,
                                         vaneRef, Solver::DefaultTrailSpanFactor * 2.0 * shroudInner,
                                         Solver::AnalyticSectionModel{}, 0.0, Rho);
}

void TestConvergesAndBudgets() {
    const auto joint = SolveJoint({});
    const double jetTorque = -(joint.Rotor.InducedMoment.x + joint.Rotor.ProfileMoment.x);

    std::cout << "two-way hover: outer=" << joint.OuterIterations << " residual=" << joint.Residual
              << " swirlFactor=" << joint.SwirlFactor << "\n  rotor T=" << -joint.Rotor.Base.Di
              << " N  jet Q=" << jetTorque << " N*m  vane Mx=" << joint.Vanes.Base.Mx
              << " N*m  vane drag=" << joint.Vanes.Base.Di << " N\n";

    CHECK(joint.Converged, "the outer rotor-vane iteration must converge");
    CHECK(joint.OuterIterations >= 2, "convergence must be judged across at least two passes");

    // The headline: the momentum budget. The vanes' recovered torque may
    // not exceed the angular-momentum flux the jet delivers.
    CHECK(joint.Vanes.Base.Mx <= jetTorque * 1.05 + 1e-9,
          "the vanes cannot recover more torque than the jet carries");
    CHECK(joint.Vanes.Base.Mx > 0.0, "they must still recover SOME of the swirl");
    CHECK(joint.SwirlFactor < 1.0,
          "on this high-solidity cruciform the budget must actually bind");

    // Cruciform symmetry must survive the coupling.
    const double thrust = -(joint.Rotor.Base.Di + joint.Vanes.Base.Di);
    CHECK(std::fabs(joint.Vanes.Base.Y) < 0.02 * thrust + 1e-6, "symmetry must cancel vane Fy");
    CHECK(std::fabs(joint.Vanes.Base.L) < 0.02 * thrust + 1e-6, "symmetry must cancel vane Fz");
}

void TestRotorFeelsTheVanes() {
    // The point of two-way: the rotor loading with vanes present must
    // differ from the vane-free rotor. Solve the same rotor+duct alone
    // and compare.
    const auto joint = SolveJoint({});

    const auto prop = MakeTestProp();
    const double shroudInner = 1.04 * prop.Radius;
    const double shroudChord = 0.5 * prop.Radius;
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega);
    const auto strips = PanelBuilder::BuildPropellerStrips(prop);
    const auto duct = PanelBuilder::BuildPropellerDuct(shroudInner, 1.2 * prop.Radius, shroudChord);

    Solver::FreestreamConditions fc;
    fc.Vinf = HoverFloorSpeed;
    fc.alphaDeg = 0.0;
    fc.betaDeg = 0.0;
    fc.rho = Rho;
    fc.p = Omega;
    fc.RefPoint = Solver::Vec3(0.0, 0.0, 0.0);
    Solver::ReferenceGeometry ref;
    for (const auto& panel : panels) ref.Area += panel.PlanformArea;
    ref.Span = 2.0 * prop.Radius;
    ref.Chord = ref.Area / ref.Span;
    const auto alone =
        Solver::SolveViscousCoupled(panels, strips, fc, ref,
                                    Solver::DefaultTrailSpanFactor * 2.0 * prop.Radius,
                                    Solver::AnalyticSectionModel{}, {}, duct);

    const double shift = std::fabs((-joint.Rotor.Base.Di) - (-alone.Base.Di));
    std::cout << "rotor thrust alone=" << -alone.Base.Di << " N  with vanes=" << -joint.Rotor.Base.Di
              << " N  (shift " << shift << " N)\n";
    CHECK(shift > 1e-3, "the rotor must feel the vanes' blockage/upwash -- that IS two-way");
    CHECK(shift < 0.3 * std::fabs(alone.Base.Di),
          "but the feedback must remain a perturbation, not an upheaval");
}

void TestControlResponseSurvives() {
    const auto neutral = SolveJoint({});
    const auto plus = SolveJoint({8.0, 0.0, 0.0, 0.0});
    const auto minus = SolveJoint({-8.0, 0.0, 0.0, 0.0});
    std::cout << "deflected outer: +8 residual=" << plus.Residual << " outer=" << plus.OuterIterations
              << " s=" << plus.SwirlFactor << " vaneIters=" << plus.Vanes.Iterations
              << " vaneRes=" << plus.Vanes.MaxResidual << "\n                 -8 residual="
              << minus.Residual << " outer=" << minus.OuterIterations << " s=" << minus.SwirlFactor
              << " vaneIters=" << minus.Vanes.Iterations << " vaneRes=" << minus.Vanes.MaxResidual
              << "\n";
    CHECK(plus.Converged && minus.Converged, "deflected joint solves must converge");

    const double thrust = -(neutral.Rotor.Base.Di + neutral.Vanes.Base.Di);
    const double dFzPlus = plus.Vanes.Base.L - neutral.Vanes.Base.L;
    const double dFzMinus = minus.Vanes.Base.L - neutral.Vanes.Base.L;
    std::cout << "coupled control response: dFz(+8)=" << dFzPlus << " N  dFz(-8)=" << dFzMinus
              << " N  (T=" << thrust << " N)\n";
    CHECK(dFzPlus * dFzMinus < 0.0, "opposite commands must still pull opposite ways");
    CHECK(std::fabs(dFzPlus - dFzMinus) > 0.004 * thrust,
          "the control differential must survive the budget");
}

} // namespace

int main() {
    TestConvergesAndBudgets();
    TestRotorFeelsTheVanes();
    TestControlResponseSurvives();

    if (failures == 0) { std::cout << "PASS: TestRotorVaneCoupling\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestRotorVaneCoupling\n";
    return 1;
}
