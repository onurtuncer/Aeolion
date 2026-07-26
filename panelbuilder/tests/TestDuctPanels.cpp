// TestDuctPanels.cpp -- validates the duct panelling against physics that
// holds regardless of the specific chord/diameters, the same way
// TestBodyPanels.cpp validates the fuselage:
//
//   - Closure. The outer wall, inner bore wall, and two end caps together
//     enclose exactly the annular ring's material volume, and every panel
//     normal faces away from that material.
//   - d'Alembert. A closed body in steady potential flow feels no net force,
//     regardless of shape -- sensitive to the frame conversion, the winding,
//     and the solve all at once.
//   - Frame. The contract states x forward; the solver is x aft, so the
//     duct's leading edge must end up ahead of its trailing edge in solver
//     axes.
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/Solver.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using Aeolion::Geometry::DuctGeometry;
using Aeolion::Geometry::HandoffContract;
using Aeolion::Lattice::SourcePanel;
using Aeolion::Math::Vec3;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

namespace PB = Aeolion::PanelBuilder;

// A minimal contract carrying only a duct (plus the wing every contract must
// have). Center is off the origin on purpose, so a test that only passes at
// (0,0,0) would be caught.
HandoffContract DuctOnlyContract(double chord, double innerDiameter, double outerDiameter,
                                 Vec3 center = Vec3(-0.4, 0.0, 0.0)) {
    HandoffContract contract;
    Aeolion::Geometry::PlanformStation root;
    root.Eta = 0.0; root.Chord = 0.2;
    Aeolion::Geometry::PlanformStation tip = root;
    tip.Eta = 1.0;
    contract.Span = 1.0;
    contract.Stations = {root, tip};
    contract.Mesh.SpanwisePanelsPerSection = 2;
    contract.Mesh.ChordwisePanels = 1;

    DuctGeometry duct;
    duct.IsStated = true;
    duct.Chord = chord;
    duct.InnerDiameter = innerDiameter;
    duct.OuterDiameter = outerDiameter;
    duct.Center = center;
    contract.Duct = duct;
    return contract;
}

Aeolion::Solver::SolveResult SolveDuct(const std::vector<SourcePanel>& duct, double alphaDeg, double Vinf) {
    Aeolion::Solver::PanelSystem system{{}, duct};
    const auto prepared = Aeolion::Solver::Prepare(system, 1.0);

    Aeolion::Solver::FreestreamConditions fc;
    fc.Vinf = Vinf; fc.alphaDeg = alphaDeg; fc.rho = 1.225;
    Aeolion::Solver::ReferenceGeometry ref;
    ref.Area = 1.0; ref.Span = 1.0; ref.Chord = 1.0;
    return Aeolion::Solver::SolveWithSystem(prepared, fc, ref);
}

void TestFrameConversion() {
    const double chord = 0.1, innerD = 0.2, outerD = 0.24;
    const Vec3 center(-0.4, 0.0, 0.0);
    const auto contract = DuctOnlyContract(chord, innerD, outerD, center);
    const auto duct = PB::LatticeBuilder(contract).BuildDuct();
    CHECK(!duct.empty(), "a duct contract should produce duct panels");

    double minX = 1e30, maxX = -1e30;
    for (const auto& panel : duct)
        for (const Vec3& corner : panel.Corners) {
            minX = std::min(minX, corner.x);
            maxX = std::max(maxX, corner.x);
        }
    // Contract leading edge is at center.x + chord/2 (forward); the solver is
    // x aft, so the leading edge maps to the SMALLER solver x.
    const double expectedLead = -(center.x + chord * 0.5);
    const double expectedTrail = -(center.x - chord * 0.5);
    CHECK(std::fabs(minX - expectedLead) < 1e-9,
          "the leading edge should sit at x = " + std::to_string(expectedLead) + " in solver axes, got " +
              std::to_string(minX));
    CHECK(std::fabs(maxX - expectedTrail) < 1e-9,
          "the trailing edge should sit at x = " + std::to_string(expectedTrail) + " in solver axes, got " +
              std::to_string(maxX));
    CHECK(minX < maxX, "the leading edge must be ahead of the trailing edge in solver axes");
}

void TestClosureAndOutwardNormals() {
    const double chord = 0.1, innerD = 0.2, outerD = 0.24;
    const Vec3 center(-0.4, 0.05, -0.02); // off-axis, to catch a hardcoded (0,0,0) assumption
    const auto contract = DuctOnlyContract(chord, innerD, outerD, center);
    PB::LatticeOptions options;
    options.DuctCircumferentialPanels = 32;
    options.DuctAxialPanels = 5;
    const auto duct = PB::LatticeBuilder(contract, options).BuildDuct();

    const int sectors = options.DuctCircumferentialPanels;
    const int axial = options.DuctAxialPanels;
    CHECK(duct.size() == static_cast<std::size_t>(2 * axial * sectors + 2 * sectors),
          "panel count should be two cylinders' worth plus two capped annuli");

    // Every normal must face away from the ring's own material: outward on
    // the outer wall, inward (toward the axis) on the bore wall, and
    // axially off the two caps.
    // The duct axis is offset by (center.y, -center.z) in solver axes (y is
    // unchanged by the frame conversion, z flips sign -- see BuildDuct()), so
    // the radial direction off that axis is the control point's offset from
    // (center.y, -center.z), not from the origin.
    int wrongFacing = 0;
    for (const auto& panel : duct) {
        const Vec3 radial(0.0, panel.ControlPoint.y - center.y, panel.ControlPoint.z + center.z);
        if (panel.Surface == PB::DuctOuterSurfaceName) {
            if (Aeolion::Math::Dot(panel.Normal, radial.Normalized()) < -1e-9) ++wrongFacing;
        } else if (panel.Surface == PB::DuctInnerSurfaceName) {
            if (Aeolion::Math::Dot(panel.Normal, radial.Normalized()) > 1e-9) ++wrongFacing;
        } else if (panel.Surface == PB::DuctLeadingCapSurfaceName) {
            if (panel.Normal.x > -1e-9) ++wrongFacing; // must face forward (-x, solver axes)
        } else if (panel.Surface == PB::DuctTrailingCapSurfaceName) {
            if (panel.Normal.x < 1e-9) ++wrongFacing; // must face aft (+x, solver axes)
        }
    }
    CHECK(wrongFacing == 0, std::to_string(wrongFacing) + " duct panels face the wrong way");

    // Divergence theorem: the enclosed volume must be the annular ring's
    // material volume, pi * (rOuter^2 - rInner^2) * chord.
    double enclosed = 0.0;
    for (const auto& panel : duct)
        enclosed += Aeolion::Math::Dot(panel.ControlPoint, panel.Normal) * panel.Area / 3.0;

    const double rInner = innerD * 0.5, rOuter = outerD * 0.5;
    const double exact = std::numbers::pi * (rOuter * rOuter - rInner * rInner) * chord;
    std::cout << "duct: " << duct.size() << " panels, enclosed volume " << enclosed << " vs exact " << exact
              << " (" << 100.0 * std::fabs(enclosed - exact) / exact << "% off)\n";
    CHECK(std::fabs(enclosed - exact) < 0.02 * exact,
          "the panelled duct should enclose the annulus's volume within 2%");
}

void TestDAlembert() {
    const auto contract = DuctOnlyContract(0.09135, 0.209, 0.225);
    PB::LatticeOptions options;
    options.DuctCircumferentialPanels = 24;
    options.DuctAxialPanels = 4;
    const auto duct = PB::LatticeBuilder(contract, options).BuildDuct();

    const double Vinf = 20.0;
    const double dynamicPressure = 0.5 * 1.225 * Vinf * Vinf;
    const double reference = dynamicPressure * std::numbers::pi * 0.225 * 0.225;

    for (double alphaDeg : {0.0, 6.0, 12.0}) {
        const auto result = SolveDuct(duct, alphaDeg, Vinf);
        const double netForce = std::sqrt(result.L * result.L + result.Di * result.Di + result.Y * result.Y);
        std::cout << "duct alpha=" << alphaDeg << ": |F|/(q*A) = " << netForce / reference << "\n";
        CHECK(netForce < 0.05 * reference,
              "a closed duct must feel essentially no net force (d'Alembert) at alpha=" +
                  std::to_string(alphaDeg) + ", got |F|/(qA) = " + std::to_string(netForce / reference));
    }
}

void TestFlagAndAbsentDuct() {
    const auto contract = DuctOnlyContract(0.1, 0.2, 0.24);

    PB::LatticeOptions off;
    off.IncludeDuct = false;
    CHECK(PB::LatticeBuilder(contract, off).BuildDuct().empty(),
          "IncludeDuct=false must produce no duct panels");

    HandoffContract noDuct = contract;
    noDuct.Duct = DuctGeometry{};
    CHECK(PB::LatticeBuilder(noDuct).BuildDuct().empty(),
          "a contract without a duct must produce no duct panels");
}

void TestFixtureDuctParses() {
    const auto contract = Aeolion::Geometry::LoadHandoff(std::string(AEOLION_TEST_DATA_DIR) +
                                                        "/AeolionGeometryHandoff-1.8.0.json");
    const auto duct = PB::LatticeBuilder(contract).BuildDuct();
    CHECK(!duct.empty(), "the 1.8.0 fixture's duct block should produce duct panels");

    // The fixture's duct sits well aft of the fuselage nose and is much
    // smaller in radius than the propeller disk it surrounds -- a loose
    // sanity bound that would catch a gross unit or frame mistake.
    for (const auto& panel : duct)
        CHECK(panel.ControlPoint.x > 0.0, "the fixture's duct must sit aft of the nose in solver axes");
}

} // namespace

int main() {
    TestFrameConversion();
    TestClosureAndOutwardNormals();
    TestDAlembert();
    TestFlagAndAbsentDuct();
    TestFixtureDuctParses();

    if (failures == 0) { std::cout << "PASS: TestDuctPanels\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestDuctPanels\n";
    return 1;
}
