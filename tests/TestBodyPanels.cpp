// TestBodyPanels.cpp -- validates the fuselage panelling against physics
// that holds regardless of the body's shape, so the checks do not encode
// the builder's own arithmetic:
//
//   - Closure. A capped body of revolution must enclose its stated volume
//     and its panel normals must all point outward. A single inverted
//     winding would flip that panel's pressure contribution and is
//     otherwise very hard to see.
//   - d'Alembert. A closed body in steady potential flow feels NO net
//     force, at any incidence. This is exact, it is independent of the
//     body's shape, and it is sensitive to the frame conversion, the
//     winding, the kernel sign and the solve all at once.
//   - Munk moment. The same body DOES feel a pitching moment at incidence,
//     and it is destabilizing (nose-up for nose-up incidence). That is the
//     dominant thing a fuselage does to an airframe's stability and the
//     reason for modelling it at all.
//   - Frame. The contract states x forward; the solver is x aft. The nose
//     must end up AHEAD of the tail in solver axes.
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/Solver.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using Aeolion::Geometry::BodyGeometry;
using Aeolion::Geometry::BodyStation;
using Aeolion::Geometry::HandoffContract;
using Aeolion::Lattice::SourcePanel;
using Aeolion::Math::Vec3;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

namespace PB = Aeolion::PanelBuilder;

// A closed ellipsoid of revolution, stated the way the contract states a
// body: x forward, so the nose is at 0 and the tail at -length.
HandoffContract EllipsoidContract(double length, double maxRadius, int axialStations) {
    HandoffContract contract;
    // A minimal but valid wing, since a contract always carries one.
    Aeolion::Geometry::PlanformStation root;
    root.Eta = 0.0; root.Chord = 0.2;
    Aeolion::Geometry::PlanformStation tip = root;
    tip.Eta = 1.0;
    contract.Span = 1.0;
    contract.Stations = {root, tip};
    contract.Mesh.SpanwisePanelsPerSection = 2;
    contract.Mesh.ChordwisePanels = 1;

    BodyGeometry body;
    body.Length = length;
    for (int i = 0; i <= axialStations; ++i) {
        const double t = static_cast<double>(i) / axialStations; // 0 at nose, 1 at tail
        BodyStation station;
        station.x = -t * length;
        // Ellipse closing to a point at BOTH ends.
        const double u = 2.0 * t - 1.0;
        station.Radius = maxRadius * std::sqrt(std::max(0.0, 1.0 - u * u));
        body.Stations.push_back(station);
    }
    contract.Body = body;
    return contract;
}

Aeolion::Solver::SolveResult SolveBody(const std::vector<SourcePanel>& body, double alphaDeg, double Vinf) {
    Aeolion::Solver::PanelSystem system{{}, body};
    const auto prepared = Aeolion::Solver::Prepare(system, 1.0);

    Aeolion::Solver::FreestreamConditions fc;
    fc.Vinf = Vinf; fc.alphaDeg = alphaDeg; fc.rho = 1.225;
    Aeolion::Solver::ReferenceGeometry ref;
    ref.Area = 1.0; ref.Span = 1.0; ref.Chord = 1.0;
    return Aeolion::Solver::SolveWithSystem(prepared, fc, ref);
}

void TestFrameConversion() {
    const auto contract = EllipsoidContract(2.0, 0.3, 16);
    const auto body = PB::LatticeBuilder(contract).BuildBody();
    CHECK(!body.empty(), "an ellipsoid contract should produce body panels");

    double minX = 1e30, maxX = -1e30;
    for (const auto& panel : body) {
        for (const Vec3& corner : panel.Corners) {
            minX = std::min(minX, corner.x);
            maxX = std::max(maxX, corner.x);
        }
    }
    // Contract nose is x = 0 and tail x = -length; the solver is x aft, so
    // the body must span 0 .. +length.
    CHECK(std::fabs(minX) < 1e-9, "the nose should sit at x = 0 in solver axes, got " + std::to_string(minX));
    CHECK(std::fabs(maxX - 2.0) < 1e-9,
          "the tail should sit at x = +length in solver axes, got " + std::to_string(maxX));
}

void TestClosureAndOutwardNormals() {
    const double length = 2.0;
    const double maxRadius = 0.3;
    const auto contract = EllipsoidContract(length, maxRadius, 40);
    PB::LatticeOptions options;
    options.BodyCircumferentialPanels = 32;
    const auto body = PB::LatticeBuilder(contract, options).BuildBody();

    // Every normal must point away from the body axis (or aft, on a base).
    // The axis runs along x, so "outward" means a non-negative radial
    // component -- an inverted winding shows up immediately.
    int inwardFacing = 0;
    for (const auto& panel : body) {
        const Vec3 radial(0.0, panel.ControlPoint.y, panel.ControlPoint.z);
        if (radial.Norm() < 1e-12) continue;
        if (Aeolion::Math::Dot(panel.Normal, radial.Normalized()) < -1e-9) ++inwardFacing;
    }
    CHECK(inwardFacing == 0, std::to_string(inwardFacing) + " body panels have inward-facing normals");

    // Divergence theorem: the volume enclosed is (1/3) * sum over panels of
    // (r . n) * area, which must match the ellipsoid's 4/3 pi a b^2.
    double enclosed = 0.0;
    for (const auto& panel : body)
        enclosed += Aeolion::Math::Dot(panel.ControlPoint, panel.Normal) * panel.Area / 3.0;

    const double exact = 4.0 / 3.0 * std::numbers::pi * (length * 0.5) * maxRadius * maxRadius;
    // The centroid of the enclosed volume is offset from the origin here, but
    // the divergence identity is origin-independent for a CLOSED surface, so
    // any discrepancy is discretization, not placement.
    std::cout << "ellipsoid: " << body.size() << " panels, enclosed volume " << enclosed << " vs exact "
              << exact << " (" << 100.0 * std::fabs(enclosed - exact) / exact << "% off)\n";
    CHECK(std::fabs(enclosed - exact) < 0.02 * exact,
          "the panelled body should enclose the ellipsoid's volume within 2%");
}

void TestDAlembertAndMunkMoment() {
    const auto contract = EllipsoidContract(2.0, 0.3, 32);
    PB::LatticeOptions options;
    options.BodyCircumferentialPanels = 24;
    const auto body = PB::LatticeBuilder(contract, options).BuildBody();

    const double Vinf = 20.0;
    const double dynamicPressure = 0.5 * 1.225 * Vinf * Vinf;
    // Scale forces on something with the body's own size, not the dummy
    // reference area, so the "near zero" below means near zero physically.
    const double reference = dynamicPressure * std::numbers::pi * 0.3 * 0.3;

    for (double alphaDeg : {0.0, 6.0, 12.0}) {
        const auto result = SolveBody(body, alphaDeg, Vinf);
        const double netForce = std::sqrt(result.L * result.L + result.Di * result.Di + result.Y * result.Y);

        std::cout << "ellipsoid alpha=" << alphaDeg << ": |F|/(q*A) = " << netForce / reference
                  << ", Cm-ish My/(q*A*L) = " << result.My / (reference * 2.0) << "\n";

        // d'Alembert: a closed body in potential flow feels no net force.
        CHECK(netForce < 0.05 * reference,
              "a closed body must feel essentially no net force (d'Alembert) at alpha=" +
                  std::to_string(alphaDeg) + ", got |F|/(qA) = " + std::to_string(netForce / reference));
    }

    // Munk moment: zero at zero incidence by symmetry, and destabilizing
    // (nose-up) once there is incidence.
    const auto atZero = SolveBody(body, 0.0, Vinf);
    const auto atAlpha = SolveBody(body, 10.0, Vinf);
    CHECK(std::fabs(atZero.My) < 0.02 * reference * 2.0,
          "a body of revolution at zero incidence should feel no pitching moment");
    CHECK(std::fabs(atAlpha.My) > 0.05 * reference * 2.0,
          "a body at incidence MUST feel a Munk moment, got My = " + std::to_string(atAlpha.My));
    CHECK(atAlpha.My > 0.0,
          "the Munk moment should be nose-up (destabilizing) for nose-up incidence, got My = " +
              std::to_string(atAlpha.My));

    // Slender-body theory says the moment goes as sin(2*alpha), not as alpha.
    // That is a sharp, shape-independent signature: if the pressure
    // integration or the crossflow were wrong, the moment would still grow
    // with incidence but would not follow this curve.
    const double atSix = SolveBody(body, 6.0, Vinf).My;
    const double atTwelve = SolveBody(body, 12.0, Vinf).My;
    const double predicted = std::sin(2.0 * 12.0 * std::numbers::pi / 180.0) /
                             std::sin(2.0 * 6.0 * std::numbers::pi / 180.0);
    const double observed = atTwelve / atSix;
    std::cout << "Munk moment ratio My(12deg)/My(6deg) = " << observed << " vs sin(2a) prediction "
              << predicted << "\n";
    CHECK(std::fabs(observed - predicted) < 0.02 * predicted,
          "the Munk moment should follow sin(2*alpha): got ratio " + std::to_string(observed) +
              ", expected " + std::to_string(predicted));
}

void TestFlagAndAbsentBody() {
    const auto contract = EllipsoidContract(2.0, 0.3, 12);

    PB::LatticeOptions off;
    off.IncludeBody = false;
    CHECK(PB::LatticeBuilder(contract, off).BuildBody().empty(),
          "IncludeBody=false must produce no body panels");

    // A contract with no body block at all: same degenerate result, no
    // special-casing needed by the caller.
    HandoffContract wingOnly = contract;
    wingOnly.Body = BodyGeometry{};
    CHECK(PB::LatticeBuilder(wingOnly).BuildBody().empty(),
          "a contract without a body must produce no body panels");
}

// The real airframe: an open base, which is the case the ellipsoid does not
// cover.
void TestFixtureBodyIsCapped() {
    const auto contract = Aeolion::Geometry::LoadHandoff(std::string(AEOLION_TEST_DATA_DIR) +
                                                        "/AeolionGeometryHandoff-1.4.0.json");
    const auto body = PB::LatticeBuilder(contract).BuildBody();
    CHECK(!body.empty(), "the 1.4.0 fixture should produce body panels");

    int basePanels = 0;
    for (const auto& panel : body)
        if (panel.Surface == PB::BodyBaseSurfaceName) ++basePanels;
    CHECK(basePanels > 0, "the fixture's open tail must be capped, so the body encloses a volume");

    // Capped, so the divergence identity should give a sensible positive
    // volume rather than the nonsense an open surface produces.
    double enclosed = 0.0;
    for (const auto& panel : body)
        enclosed += Aeolion::Math::Dot(panel.ControlPoint, panel.Normal) * panel.Area / 3.0;
    std::cout << "fixture body: " << body.size() << " panels (" << basePanels
              << " base), enclosed volume " << enclosed << " m^3\n";
    CHECK(enclosed > 0.0, "a capped body must enclose a positive volume");
    // Sanity bound: it cannot exceed the cylinder that contains it.
    const double bound = std::numbers::pi * 0.0491 * 0.0491 * 0.49099;
    CHECK(enclosed < bound, "enclosed volume must be less than the bounding cylinder");
}

} // namespace

int main() {
    TestFrameConversion();
    TestClosureAndOutwardNormals();
    TestDAlembertAndMunkMoment();
    TestFlagAndAbsentBody();
    TestFixtureBodyIsCapped();

    if (failures == 0) { std::cout << "PASS: TestBodyPanels\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestBodyPanels\n";
    return 1;
}
