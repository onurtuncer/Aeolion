// TestVaneCascade.cpp -- the cascade momentum closure for duct-jet vane
// loads (Solver::SolveVaneCascade) and its role as the coupled driver's
// default vane closure.
//
// The closure exists because the vane lattice at hover is bistable once
// deflected: mutually-inflating strip induction let a single 8-degree
// vane carry more side force than the whole net thrust on one solution
// branch. The cascade bounds every force by its sector's mass flow BY
// CONSTRUCTION, so what must hold here is exactly that bookkeeping: the
// undeflected cruciform cannot recover more torque than the jet's
// angular-momentum flux; a dead-air strip carries nothing; deflection
// makes bounded, sign-correct control forces; and the coupled solve
// under the closure converges with no runaway anywhere in it.
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/RotorVaneCoupling.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/VaneCascade.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <cmath>
#include <functional>
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
constexpr double ExitRadius = 0.156;
constexpr double VaneChord = 0.075;

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

// A prescribed jet: uniform axial core with solid-body swirl, zero
// outside the jet radius (dead air), so the momentum bookkeeping can be
// checked against a hand integral.
std::function<Solver::Vec3(const Solver::Vec3&)> MakeJet(double axial, double tipSwirl,
                                                         double jetRadius) {
    return [=](const Math::Vec3& point) {
        const double radius = std::hypot(point.y, point.z);
        if (radius > jetRadius) return Math::Vec3(0.0, 0.0, 0.0);
        const Math::Vec3 tangent(0.0, -point.z, point.y);
        const double swirl = tipSwirl * radius / jetRadius;
        return Math::Vec3(axial, 0.0, 0.0) +
               (radius > 1e-9 ? tangent * (swirl / radius) : Math::Vec3());
    };
}

Solver::ViscousCoupledResult SolveCascade(const std::vector<double>& deflectionsDeg,
                                          const std::function<Solver::Vec3(const Solver::Vec3&)>& jet) {
    const auto surfaces = MakeVaneSet();
    const auto panels = PanelBuilder::BuildDuctVanes(surfaces, ExitRadius, 0.5 * VaneChord,
                                                     VaneChord, deflectionsDeg, jet);
    const auto strips = PanelBuilder::BuildDuctVaneStrips(surfaces, ExitRadius, VaneChord,
                                                          deflectionsDeg);
    Solver::FreestreamConditions fc;
    fc.Vinf = HoverFloorSpeed;
    fc.alphaDeg = 0.0;
    fc.betaDeg = 0.0;
    fc.rho = Rho;
    fc.RefPoint = Solver::Vec3(0.0, 0.0, 0.0);
    return Solver::SolveVaneCascade(panels, strips, fc, jet, Solver::AnalyticSectionModel{});
}

// The jet's angular-momentum flux through the vane annulus, integrated
// from the same prescribed field the closure sees.
double JetFlux(double axial, double tipSwirl, double jetRadius) {
    const double inner = 0.4 * ExitRadius;
    const double outer = std::min(ExitRadius, jetRadius);
    double flux = 0.0;
    const int bands = 400;
    for (int i = 0; i < bands; ++i) {
        const double r = inner + (outer - inner) * (i + 0.5) / bands;
        const double dr = (outer - inner) / bands;
        const double w = tipSwirl * r / jetRadius;
        flux += Rho * axial * (2.0 * std::numbers::pi * r * dr) * w * r;
    }
    return flux;
}

void TestMomentumBoundByConstruction() {
    const double axial = 12.0, tipSwirl = 4.0;
    const auto run = SolveCascade({}, MakeJet(axial, tipSwirl, ExitRadius));
    const double flux = JetFlux(axial, tipSwirl, ExitRadius);
    std::cout << "cascade undeflected: Mx=" << run.Base.Mx << " N*m  (jet flux " << flux
              << " N*m)  Fx=" << run.Base.Di << " N  Fy=" << run.Base.Y << "  Fz=" << run.Base.L
              << "\n";
    CHECK(run.Converged, "the cascade closure is a direct solve and always converges");
    CHECK(run.Base.Mx > 0.0, "vanes straightening a +x-wise swirl must recover counter-torque");
    CHECK(run.Base.Mx <= 1.05 * flux,
          "recovered torque may not exceed the jet's angular-momentum flux -- THE bound");
    CHECK(run.Base.Mx > 0.5 * flux,
          "a full-chord cruciform should recover the bulk of a modest swirl");
    const double thrustScale = Rho * axial * axial * ExitRadius * ExitRadius;
    CHECK(std::fabs(run.Base.Y) < 0.01 * thrustScale + 1e-9, "cruciform symmetry must cancel Fy");
    CHECK(std::fabs(run.Base.L) < 0.01 * thrustScale + 1e-9, "cruciform symmetry must cancel Fz");
}

void TestDeadAirCarriesNothing() {
    // Contract the jet to 70% of the exit radius: the outboard strips sit
    // in still air and must carry ~no force -- the jet-edge pathology
    // that plagued the lattice closes itself here, because force follows
    // mass flow.
    const double axial = 12.0, tipSwirl = 4.0;
    const auto full = SolveCascade({}, MakeJet(axial, tipSwirl, ExitRadius));
    const auto contracted = SolveCascade({}, MakeJet(axial, tipSwirl, 0.7 * ExitRadius));
    std::cout << "contracted jet: Mx " << full.Base.Mx << " -> " << contracted.Base.Mx << " N*m\n";
    CHECK(contracted.Base.Mx < full.Base.Mx,
          "a contracted jet wets less vane span and must recover less torque");
    CHECK(contracted.Base.Mx > 0.0, "but the wetted strips must still work");
    for (std::size_t i = 0; i < contracted.Strips.size(); ++i) {
        const auto& strip = contracted.Strips[i];
        const double radius = std::hypot(strip.Mid.y, strip.Mid.z);
        if (radius > 0.75 * ExitRadius)
            CHECK(strip.Force.Norm() < 1e-9,
                  "a strip in dead air outside the jet must carry no force");
    }
}

void TestBoundedControlResponse() {
    // Deflect vane 0 both ways in a swirl-rich jet -- the exact regime
    // where the lattice bifurcated into runaway branches. The response
    // must flip sign with the command and stay bounded by the sector
    // momentum the jet can deliver.
    const double axial = 9.0, tipSwirl = 6.0;
    const auto jet = MakeJet(axial, tipSwirl, ExitRadius);
    const auto neutral = SolveCascade({}, jet);
    const auto plus = SolveCascade({8.0, 0.0, 0.0, 0.0}, jet);
    const auto minus = SolveCascade({-8.0, 0.0, 0.0, 0.0}, jet);

    const double dFzPlus = plus.Base.L - neutral.Base.L;
    const double dFzMinus = minus.Base.L - neutral.Base.L;
    std::cout << "cascade control response: dFz(+8)=" << dFzPlus << " N  dFz(-8)=" << dFzMinus
              << " N\n";
    CHECK(dFzPlus * dFzMinus < 0.0, "opposite commands must pull opposite ways");
    CHECK(std::fabs(dFzPlus) > 1e-4 && std::fabs(dFzMinus) > 1e-4,
          "the response must be real, not numerically dead");

    // The bound that the lattice violated: one vane's whole-force change
    // cannot exceed the momentum flux of its own sector.
    const double sectorFlux = Rho * axial * (std::numbers::pi * ExitRadius * ExitRadius / 4.0) *
                              std::hypot(axial, tipSwirl);
    CHECK(std::fabs(dFzPlus) < sectorFlux && std::fabs(dFzMinus) < sectorFlux,
          "one vane's response must stay inside its sector's momentum flux");
}

void TestCoupledSolveUnderCascade() {
    // The production path: the two-way rotor-vane solve with the cascade
    // closure (the default), on the reference hover configuration. No
    // budget factor is needed -- the closure conserves angular momentum
    // per sector -- and no state anywhere may run away.
    Geometry::Propeller prop;
    prop.BladeCount = 2;
    prop.Radius = 0.15;
    prop.HubRadius = 0.02;
    for (int i = 0; i < 12; ++i) {
        const double r = prop.HubRadius + (prop.Radius - prop.HubRadius) * (i + 0.5) / 13.0;
        const double eta = r / prop.Radius;
        prop.Stations.push_back({r, 0.025 * (1.0 - 0.3 * eta), 30.0 - 20.0 * eta});
    }
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

    Solver::FreestreamConditions vaneFc = fc;
    vaneFc.p = 0.0;
    Solver::ReferenceGeometry vaneRef;
    vaneRef.Area = 1.0;
    vaneRef.Span = 2.0 * shroudInner;
    vaneRef.Chord = shroudChord;

    const auto solveJoint = [&](const std::vector<double>& deflections) {
        const Solver::VaneBuilder vaneBuilder =
            [&](const std::function<Solver::Vec3(const Solver::Vec3&)>& flow)
            -> std::pair<std::vector<Solver::Panel>, std::vector<Solver::StripSection>> {
            return {PanelBuilder::BuildDuctVanes(surfaces, shroudInner, 0.5 * shroudChord,
                                                 shroudChord, deflections, flow),
                    PanelBuilder::BuildDuctVaneStrips(surfaces, shroudInner, shroudChord,
                                                      deflections)};
        };
        const Solver::RotorBuilder rotorBuilder =
            [&](const std::function<double(double)>& inflow) {
                return PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega, inflow);
            };
        return Solver::SolveRotorVaneCoupled(
            panels, strips, fc, ref, trail, Solver::AnalyticSectionModel{}, duct, vaneBuilder,
            vaneFc, vaneRef, Solver::DefaultTrailSpanFactor * 2.0 * shroudInner,
            Solver::AnalyticSectionModel{}, 0.0, Rho, {}, rotorBuilder);
    };

    const auto neutral = solveJoint({});
    const double thrust = -(neutral.Rotor.Base.Di + neutral.Vanes.Base.Di);
    const double jetTorque = -(neutral.Rotor.InducedMoment.x + neutral.Rotor.ProfileMoment.x);
    std::cout << "coupled cascade hover: outer=" << neutral.OuterIterations
              << " residual=" << neutral.Residual << "  T=" << thrust
              << " N  jet Q=" << jetTorque << " N*m  vane Mx=" << neutral.Vanes.Base.Mx
              << " N*m\n";
    CHECK(neutral.Converged, "the coupled solve under the cascade closure must converge");
    CHECK(neutral.SwirlFactor == 1.0, "the cascade needs no swirl-budget factor");
    CHECK(neutral.Vanes.Base.Mx > 0.0, "the vanes must recover swirl");
    CHECK(neutral.Vanes.Base.Mx <= 1.1 * jetTorque,
          "recovered torque bounded by the jet's flux, no budget iteration required");
    CHECK(thrust > 0.0, "the system must thrust");
    CHECK(std::fabs(neutral.Vanes.Base.Y) < 0.05 * thrust + 1e-6, "no phantom side force");
    CHECK(std::fabs(neutral.Vanes.Base.L) < 0.05 * thrust + 1e-6, "no phantom vertical force");

    const auto plus = solveJoint({8.0, 0.0, 0.0, 0.0});
    const auto minus = solveJoint({-8.0, 0.0, 0.0, 0.0});
    CHECK(plus.Converged && minus.Converged, "deflected coupled solves must converge");
    const double dFzPlus = plus.Vanes.Base.L - neutral.Vanes.Base.L;
    const double dFzMinus = minus.Vanes.Base.L - neutral.Vanes.Base.L;
    std::cout << "coupled cascade response: dFz(+8)=" << dFzPlus << " N  dFz(-8)=" << dFzMinus
              << " N  (T=" << thrust << " N)\n";
    CHECK(dFzPlus * dFzMinus < 0.0, "opposite commands must pull opposite ways");
    // The runaway that motivated all this: a deflected vane's side force
    // exceeded the NET THRUST. Under the cascade it must stay a modest
    // fraction of it.
    CHECK(std::fabs(plus.Vanes.Base.L) < 0.5 * thrust && std::fabs(minus.Vanes.Base.L) < 0.5 * thrust,
          "a single 8-degree vane may not out-force the propulsor -- bounded by construction");
}

} // namespace

int main() {
    TestMomentumBoundByConstruction();
    TestDeadAirCarriesNothing();
    TestBoundedControlResponse();
    TestCoupledSolveUnderCascade();

    if (failures == 0) { std::cout << "PASS: TestVaneCascade\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestVaneCascade\n";
    return 1;
}
