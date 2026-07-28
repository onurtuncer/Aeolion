// TestPropellerDuct.cpp -- the ducted-propeller interaction: the shroud
// meshed by BuildPropellerDuct() sharing the viscous-coupled solve with the
// blade lattice.
//
// Mesh checks are the same physics-not-arithmetic style as TestDuctPanels:
// a closed ring must enclose exactly its annular volume and point every
// normal outward. Interaction checks pin what a shroud must DO in this
// model: the coupled solve still converges; the ducted system's in-plane
// forces still cancel by symmetry; and the duct carries a nonzero share of
// the axial force through its lip-suction pressure field -- the
// interaction being the point of panelling it at all.
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

void TestShroudMeshCloses() {
    const double rInner = 0.156, rOuter = 0.18, chord = 0.075;
    const auto panels = PanelBuilder::BuildPropellerDuct(rInner, rOuter, chord);
    CHECK(!panels.empty(), "the shroud mesh must not be empty");

    // Divergence theorem: volume = (1/3) * sum (cp . n) A over a closed
    // surface; the ring's exact volume is pi (rO^2 - rI^2) chord.
    double volume = 0.0;
    Solver::Vec3 areaSum(0, 0, 0);
    for (const auto& panel : panels) {
        volume += Solver::Dot(panel.ControlPoint, panel.Normal) * panel.Area / 3.0;
        areaSum = areaSum + panel.Normal * panel.Area;
        CHECK(std::fabs(panel.Normal.Norm() - 1.0) < 1e-9, "every shroud normal must be unit");
    }
    const double exact = std::numbers::pi * (rOuter * rOuter - rInner * rInner) * chord;
    CHECK(std::fabs(volume - exact) < 0.02 * exact,
          "the shroud must enclose its annular volume (closure + outward normals)");
    CHECK(areaSum.Norm() < 1e-6, "a closed surface's area-weighted normals must sum to zero");
}

void TestDuctedInteraction() {
    const auto prop = MakeTestProp();
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega);
    const auto strips = PanelBuilder::BuildPropellerStrips(prop);
    const auto duct = PanelBuilder::BuildPropellerDuct(1.04 * prop.Radius, 1.2 * prop.Radius,
                                                       0.5 * prop.Radius);

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

    const auto open = Solver::SolveViscousCoupled(panels, strips, fc, ref, trail,
                                                  Solver::AnalyticSectionModel{});
    const auto ducted = Solver::SolveViscousCoupled(panels, strips, fc, ref, trail,
                                                    Solver::AnalyticSectionModel{}, {}, duct);

    const double openThrust = -open.Base.Di;
    const double ductedThrust = -ducted.Base.Di;
    const double ductThrust = -ducted.SourceForce.x;
    std::cout << "hover: open T=" << openThrust << " N   ducted total T=" << ductedThrust
              << " N (duct share " << ductThrust << " N)  iters=" << ducted.Iterations << "\n";

    CHECK(ducted.Converged, "the ducted coupled solve must converge");
    CHECK(std::isfinite(ductedThrust) && std::isfinite(ductThrust), "ducted loads must be finite");
    CHECK(ductThrust != 0.0, "the duct must carry a share of the axial force -- that IS the interaction");
    CHECK(std::fabs(ducted.Base.L) < 0.05 * std::fabs(ductedThrust) + 1e-6,
          "axisymmetry must cancel the ducted system's vertical force");
    CHECK(std::fabs(ducted.Base.Y) < 0.05 * std::fabs(ductedThrust) + 1e-6,
          "axisymmetry must cancel the ducted system's side force");

    // The duct must also talk back to the blades: the blade circulation
    // distribution with the shroud present must differ measurably from the
    // open prop's.
    double maxShift = 0.0;
    for (std::size_t i = 0; i < open.Base.gamma.size(); ++i)
        maxShift = std::max(maxShift, std::fabs(open.Base.gamma[i] - ducted.Base.gamma[i]));
    std::cout << "max blade-circulation shift from the duct: " << maxShift << "\n";
    CHECK(maxShift > 1e-4, "the shroud must change the blade loading, not just add its own force");
}

} // namespace

int main() {
    TestShroudMeshCloses();
    TestDuctedInteraction();

    if (failures == 0) { std::cout << "PASS: TestPropellerDuct\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestPropellerDuct\n";
    return 1;
}
