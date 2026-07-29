// TestPropellerVanes.cpp -- the downstream duct-jet control vanes: meshed
// from contract-shaped ControlSurfaces, solved in the STATIC frame against
// the azimuthal-mean slipstream of the converged rotor+duct solution, and
// integrated into the propulsive wrench.
//
// The physics pinned here is the reason thrust-vectoring vanes exist: at
// essentially zero airspeed, a deflected vane in the propwash must make
// real side force and a real control moment (the propwash IS its dynamic
// pressure), the sign must follow the deflection, an undeflected cruciform
// must make ~none, and the vanes sitting in the swirl must recover part of
// the rotor's torque as counter-torque even undeflected.
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/SectionBoundaryLayer.h"
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

// The contract's cruciform, in miniature: four all-moving vanes on the
// cardinal hinge axes (CONTRACT frame, x-forward/z-down -- the builder
// owns the conversion), radial 0.4..1.0 of the exit radius, full duct
// chord.
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

struct VaneRun {
    Solver::ViscousCoupledResult Rotor;
    Solver::ViscousCoupledResult Vanes;
};

VaneRun Solve(const std::vector<double>& deflectionsDeg) {
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
    const double trail = Solver::DefaultTrailSpanFactor * 2.0 * prop.Radius;

    VaneRun run;
    run.Rotor = Solver::SolveViscousCoupled(panels, strips, fc, ref, trail,
                                            Solver::AnalyticSectionModel{}, {}, duct);

    const auto surfaces = MakeVaneSet();
    const auto vanePanels = PanelBuilder::BuildDuctVanes(surfaces, shroudInner, 0.5 * shroudChord,
                                                         shroudChord, deflectionsDeg);
    const auto vaneStrips =
        PanelBuilder::BuildDuctVaneStrips(surfaces, shroudInner, shroudChord, deflectionsDeg);
    const auto slipstream = Solver::SlipstreamField(run.Rotor, 0.0, Rho, duct);

    Solver::FreestreamConditions vaneFc = fc;
    vaneFc.p = 0.0; // static frame: the vanes do not rotate
    Solver::ReferenceGeometry vaneRef;
    for (const auto& panel : vanePanels) vaneRef.Area += panel.PlanformArea;
    vaneRef.Span = 2.0 * shroudInner;
    vaneRef.Chord = shroudChord;

    // The vanes get the PROPER viscous treatment -- the transpiration-coupled
    // boundary-layer section solver (flat plates: no camber callable), the
    // same model the blades default to. Its interior is piecewise (discrete
    // transition/separation stations), hence the coarser tolerance, as in
    // TestViscousCoupling.
    Solver::ViscousCouplingOptions vaneOptions;
    vaneOptions.Tolerance = 2e-2;
    vaneOptions.MaxIterations = 400;
    run.Vanes = Solver::SolveViscousCoupled(vanePanels, vaneStrips, vaneFc, vaneRef,
                                            Solver::DefaultTrailSpanFactor * 2.0 * shroudInner,
                                            Solver::BoundaryLayerSectionModel{}, vaneOptions, {},
                                            slipstream);
    std::cout << "  [vane solve: iters=" << run.Vanes.Iterations
              << " residual=" << run.Vanes.MaxResidual << "]\n";
    return run;
}

void TestMeshAndAlignment() {
    const auto surfaces = MakeVaneSet();
    const auto panels = PanelBuilder::BuildDuctVanes(surfaces, 0.156, 0.0375, 0.075, {});
    const auto strips = PanelBuilder::BuildDuctVaneStrips(surfaces, 0.156, 0.075, {});
    CHECK(panels.size() == 4u * PanelBuilder::DefaultVaneRadialPanels,
          "four vanes at one Weissinger row per radial strip");
    CHECK(strips.size() == panels.size(), "vane strips must align one-to-one with vane panels");
    for (const auto& panel : panels) {
        CHECK(std::fabs(panel.Normal.Norm() - 1.0) < 1e-9, "vane normals must be unit");
        CHECK(panel.TrailDirA.x == 1.0, "vane wakes must trail downstream");
    }
}

void TestUndeflectedCruciform() {
    const auto run = Solve({});
    for (std::size_t i = 0; i < run.Vanes.Strips.size() && i < 5; ++i) {
        const auto& st = run.Vanes.Strips[i];
        std::cout << "  strip " << i << ": r=" << std::hypot(st.Mid.y, st.Mid.z)
                  << " alphaEff=" << st.alphaEffDeg << " cl=" << st.cl << " Vrel=" << st.Vrel
                  << " R=" << st.Residual << "\n";
    }
    CHECK(run.Vanes.Converged, "the undeflected vane solve must converge");

    const double thrust = -run.Rotor.Base.Di;
    std::cout << "undeflected: vane drag=" << run.Vanes.Base.Di << " N  Fy=" << run.Vanes.Base.Y
              << "  Fz=" << run.Vanes.Base.L << "  Mx=" << run.Vanes.Base.Mx
              << " (rotor Mx=" << run.Rotor.Base.Mx << ")  My=" << run.Vanes.Base.My
              << "  Mz=" << run.Vanes.Base.Mz << "\n";

    // A symmetric cruciform at zero deflection: side forces and control
    // moments cancel...
    CHECK(std::fabs(run.Vanes.Base.Y) < 0.02 * thrust + 1e-6, "cruciform symmetry must cancel Fy");
    CHECK(std::fabs(run.Vanes.Base.L) < 0.02 * thrust + 1e-6, "cruciform symmetry must cancel Fz");

    // ...but the swirl does not: every vane sees the same tangential wash,
    // and their reaction opposes the rotor's torque on the airframe. The
    // rotor's Mx is negative (it absorbs shaft torque); the vanes' must
    // push back the other way.
    CHECK(run.Vanes.Base.Mx * run.Rotor.Base.Mx < 0.0,
          "undeflected vanes in the swirl must produce counter-torque");
}

void TestDeflectionMakesControlWrench() {
    // Deflect the +y-hinged vane (index 0) alone. The control RESPONSE is
    // measured about the neutral state: at hover the vane sits in tens of
    // degrees of swirl, so its zero-force point is nowhere near zero
    // deflection and the raw forces at +-delta are NOT antisymmetric --
    // but the response to opposite commands must still pull opposite ways,
    // and the differential must be a real fraction of the thrust.
    const auto neutral = Solve({});
    const auto plus = Solve({8.0, 0.0, 0.0, 0.0});
    const auto minus = Solve({-8.0, 0.0, 0.0, 0.0});
    CHECK(plus.Vanes.Converged && minus.Vanes.Converged, "deflected vane solves must converge");

    const double thrust = -plus.Rotor.Base.Di;
    const double dFzPlus = plus.Vanes.Base.L - neutral.Vanes.Base.L;
    const double dFzMinus = minus.Vanes.Base.L - neutral.Vanes.Base.L;
    const double dMyPlus = plus.Vanes.Base.My - neutral.Vanes.Base.My;
    const double dMyMinus = minus.Vanes.Base.My - neutral.Vanes.Base.My;
    std::cout << "vane0 response about neutral: dFz(+8)=" << dFzPlus << " N  dFz(-8)=" << dFzMinus
              << " N  dMy(+8)=" << dMyPlus << " N*m  dMy(-8)=" << dMyMinus
              << " N*m  (rotor T=" << thrust << " N)\n";

    CHECK(dFzPlus * dFzMinus < 0.0, "opposite commands must pull the side force opposite ways");
    CHECK(std::fabs(dFzPlus - dFzMinus) > 0.005 * thrust,
          "the side-force differential must be a real fraction of the thrust -- this is the\n"
          "      hover control authority the vanes exist for");
    CHECK(dMyPlus * dMyMinus < 0.0, "opposite commands must flip the control moment");
    CHECK(std::fabs(dMyPlus - dMyMinus) > 1e-4, "the moment differential must be real");

    // The wrench debit: commanding at least one way costs drag over neutral
    // (the other way may unload the swirl-biased vane and cost less).
    CHECK(std::max(plus.Vanes.Base.Di, minus.Vanes.Base.Di) > neutral.Vanes.Base.Di,
          "deflecting against the swirl must cost vane drag");
}

} // namespace

int main() {
    TestMeshAndAlignment();
    TestUndeflectedCruciform();
    TestDeflectionMakesControlWrench();

    if (failures == 0) { std::cout << "PASS: TestPropellerVanes\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestPropellerVanes\n";
    return 1;
}
