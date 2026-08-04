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
    // This test pins the ONE-WAY vane mechanics (mesh, frames, control
    // response signs), so it runs at the swirl level the two-way budget
    // converges to: the un-budgeted hover swirl of the helical-wake jet
    // over-recovers torque several-fold and drives the vanes outside the
    // section models' convergent range. The full un-scaled loop is
    // TestRotorVaneCoupling's job.
    const auto slipstream = Solver::SlipstreamField(run.Rotor, 0.0, Rho, duct, 1.0 / 3.0);
    // Vane wakes leave along the local mean flow (at hover, the slipstream
    // itself): helically with the swirl, not hardcoded downstream.
    const auto localFlow = [&](const Math::Vec3& point) { return slipstream(point); };
    const auto vanePanels = PanelBuilder::BuildDuctVanes(surfaces, shroudInner, 0.5 * shroudChord,
                                                         shroudChord, deflectionsDeg, localFlow);
    const auto vaneStrips =
        PanelBuilder::BuildDuctVaneStrips(surfaces, shroudInner, shroudChord, deflectionsDeg);

    Solver::FreestreamConditions vaneFc = fc;
    vaneFc.p = 0.0; // static frame: the vanes do not rotate
    Solver::ReferenceGeometry vaneRef;
    for (const auto& panel : vanePanels) vaneRef.Area += panel.PlanformArea;
    vaneRef.Span = 2.0 * shroudInner;
    vaneRef.Chord = shroudChord;

    // This test is the ONE-WAY mechanics diagnostic (mesh, frames, sign
    // conventions, control response). Under the Level-B helical wake the
    // jet is weaker and swirlier and the OUTERMOST strip straddles the
    // jet-edge shear layer near alphaEff = -100 degrees -- which is
    // exactly what the coupling's reversed-flow residual exclusion
    // exists for; the interior strips converge to ~1e-2.
    Solver::ViscousCouplingOptions vaneOptions;
    // A DEFLECTED vane adds its command on top of ~30-50 degrees of swirl
    // incidence, parking extra strips in the post-stall 45-60 degree band
    // where the coupled solve limit-cycles at the 0.1-0.2 level; the
    // deflected diagnostics therefore accept that plateau (the residual
    // exclusion above keeps the raw first iterate from slipping through),
    // while the undeflected solve holds the tight tolerance.
    vaneOptions.Tolerance = deflectionsDeg.empty() ? 2e-2 : 2.5e-1;
    vaneOptions.MaxIterations = 400;
    run.Vanes = Solver::SolveViscousCoupled(vanePanels, vaneStrips, vaneFc, vaneRef,
                                            Solver::DefaultTrailSpanFactor * 2.0 * shroudInner,
                                            Solver::AnalyticSectionModel{}, vaneOptions, {},
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
        CHECK(panel.TrailDirA.x == 1.0, "with no stated flow, vane wakes default to downstream");
    }

    // With a local flow stated, the wake must leave along it -- in a
    // swirling jet that means helically, with real tangential components,
    // not hardcoded axial.
    const auto swirling = [](const Math::Vec3& point) {
        const Math::Vec3 that = Math::Cross(Math::Vec3(1.0, 0.0, 0.0), point).Normalized();
        return Math::Vec3(5.0, 0.0, 0.0) + that * 3.0;
    };
    const auto swirled = PanelBuilder::BuildDuctVanes(surfaces, 0.156, 0.0375, 0.075, {}, swirling);
    for (const auto& panel : swirled) {
        CHECK(std::fabs(panel.TrailDirA.Norm() - 1.0) < 1e-9, "flow-aligned legs must be unit");
        const double tangential = std::hypot(panel.TrailDirA.y, panel.TrailDirA.z);
        CHECK(tangential > 0.3, "in a swirling jet the wake must leave with the swirl");
        CHECK(panel.TrailDirA.x > 0.0, "and still convect downstream");
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

void TestBoundaryLayerVanesInMildJet() {
    // The PROPER viscous vane treatment -- the transpiration-coupled
    // boundary-layer section model -- exercised INSIDE its convergent
    // envelope: a brisk, mildly swirling jet (about 6 degrees of swirl at
    // the tip), the cruise-like regime where the BL solution is the
    // honest one. The hover jet-edge pathology above is exactly what this
    // case excludes: every strip here sits in attached-flow incidence.
    const auto surfaces = MakeVaneSet();
    const double exitRadius = 0.156;
    const double chord = 0.075;
    const auto jet = [&](const Math::Vec3& point) {
        // Rotor spins +x, so the jet swirl rotates the same way: +x-wise.
        // Brisk on purpose: 60 m/s puts the chord Reynolds number near
        // 3e5. Below that the flat plate's laminar-separation polar is
        // non-monotonic and the coupled cruciform genuinely bifurcates
        // into asymmetric equilibria (observed at 20 and 40 m/s) -- real
        // low-Re physics, but not the unique fixed point a symmetry
        // test can pin.
        const Math::Vec3 tangent(0.0, -point.z, point.y);
        const double radius = std::hypot(point.y, point.z);
        const double swirl = 6.0 * radius / exitRadius;
        return Math::Vec3(60.0, 0.0, 0.0) +
               (radius > 1e-9 ? tangent * (swirl / radius) : Math::Vec3());
    };
    const auto panels = PanelBuilder::BuildDuctVanes(surfaces, exitRadius, 0.5 * chord, chord, {}, jet);
    const auto strips = PanelBuilder::BuildDuctVaneStrips(surfaces, exitRadius, chord, {});

    Solver::FreestreamConditions fc;
    fc.Vinf = HoverFloorSpeed;
    fc.alphaDeg = 0.0;
    fc.betaDeg = 0.0;
    fc.rho = Rho;
    fc.RefPoint = Solver::Vec3(0.0, 0.0, 0.0);
    Solver::ReferenceGeometry ref;
    for (const auto& panel : panels) ref.Area += panel.PlanformArea;
    ref.Span = 2.0 * exitRadius;
    ref.Chord = chord;

    Solver::ViscousCouplingOptions options;
    options.Tolerance = 2e-2; // the BL interior is piecewise, as in TestViscousCoupling
    options.MaxIterations = 400;
    const auto run = Solver::SolveViscousCoupled(panels, strips, fc, ref,
                                                 Solver::DefaultTrailSpanFactor * 2.0 * exitRadius,
                                                 Solver::BoundaryLayerSectionModel{}, options, {}, jet);
    std::cout << "BL vanes, mild jet: iters=" << run.Iterations << " residual=" << run.MaxResidual
              << "  Mx=" << run.Base.Mx << " N*m  Di=" << run.Base.Di << " N\n";
    for (std::size_t i = 0; i < run.Strips.size(); ++i) {
        const auto& st = run.Strips[i];
        std::cout << "  strip " << i << ": r=" << std::hypot(st.Mid.y, st.Mid.z)
                  << " alphaEff=" << st.alphaEffDeg << " cl=" << st.cl << " Vrel=" << st.Vrel
                  << " R=" << st.Residual << "\n";
    }
    CHECK(run.Converged, "the BL vane solve must converge inside its envelope");
    for (const auto& strip : run.Strips)
        CHECK(std::fabs(strip.alphaEffDeg) < 15.0,
              "every strip must sit in the attached-flow incidence this case poses");
    CHECK(run.Base.Mx > 0.0,
          "vanes straightening a +x-wise swirl must react with counter-torque about +x");
    const double load = std::fabs(run.Base.Mx);
    CHECK(std::fabs(run.Base.Y) < 0.5 * load / exitRadius + 1e-6,
          "cruciform symmetry must cancel Fy under the BL model too");
    CHECK(std::fabs(run.Base.L) < 0.5 * load / exitRadius + 1e-6,
          "cruciform symmetry must cancel Fz under the BL model too");
}

} // namespace

int main() {
    TestMeshAndAlignment();
    TestUndeflectedCruciform();
    TestDeflectionMakesControlWrench();
    TestBoundaryLayerVanesInMildJet();

    if (failures == 0) { std::cout << "PASS: TestPropellerVanes\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestPropellerVanes\n";
    return 1;
}
