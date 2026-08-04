// TestSectionBoundaryLayer.cpp -- the 2-D transpiration-coupled
// boundary-layer section solver (Solver::BoundaryLayerSectionModel),
// checked against thin-airfoil theory and boundary-layer scalings rather
// than against itself:
//
//   - flat plate at zero incidence: no lift, and a drag on the order of
//     the flat-plate friction values the correlations are fits OF;
//   - lift slope: below but near 2*pi per radian (the boundary layer can
//     only DECAMBER a section, never add lift);
//   - positive camber lifts at zero incidence;
//   - drag falls as Reynolds number rises;
//   - the transpiration feedback reduces lift relative to the same
//     lumped-vortex solve with the boundary layer silenced.
#include "Aeolion/Solver/SectionBoundaryLayer.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>

using namespace Aeolion;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

Solver::StripSection FlatStrip() {
    Solver::StripSection strip;
    strip.Chord = 0.1;
    strip.Width = 0.01;
    strip.Eta = 0.5;
    return strip;
}

// Parabolic camber z/c = 4*eps*psi*(1-psi): slope = 4*eps*(1-2*psi).
Solver::BoundaryLayerSectionModel CamberedModel(double eps) {
    Solver::BoundaryLayerSectionModel model;
    model.CamberSlope = [eps](double /*eta*/, double psi) { return 4.0 * eps * (1.0 - 2.0 * psi); };
    return model;
}

void TestFlatPlate() {
    const Solver::BoundaryLayerSectionModel model; // no camber
    const auto strip = FlatStrip();

    const auto zero = model(strip, 0.0, 5e5, 0.0);
    std::cout << "flat, alpha=0, Re=5e5: cl=" << zero.cl << "  cd=" << zero.cd << "\n";
    CHECK(std::fabs(zero.cl) < 0.01, "a flat plate at zero incidence must not lift");
    CHECK(zero.cd > 0.002 && zero.cd < 0.03,
          "flat-plate drag must land in the friction-drag range at Re=5e5");

    const auto lifted = model(strip, 4.0, 5e5, 0.0);
    const double thinAirfoil = 2.0 * std::numbers::pi * Math::DegToRad(4.0);
    std::cout << "flat, alpha=4: cl=" << lifted.cl << "  (2*pi*alpha=" << thinAirfoil << ")\n";
    CHECK(lifted.cl > 0.6 * thinAirfoil && lifted.cl < thinAirfoil,
          "viscous lift must sit below but near the thin-airfoil slope");

    const auto more = model(strip, 8.0, 5e5, 0.0);
    CHECK(more.cl > lifted.cl, "lift must grow with incidence below stall");
    CHECK(more.cd > lifted.cd, "drag must grow with incidence");
}

void TestCamberLifts() {
    const auto model = CamberedModel(0.02); // 2% parabolic camber
    const auto strip = FlatStrip();
    const auto c = model(strip, 0.0, 5e5, 0.0);
    // Inviscid thin-airfoil: cl = 4*pi*eps = 0.251 at alpha = 0.
    std::cout << "2% camber, alpha=0: cl=" << c.cl << "\n";
    CHECK(c.cl > 0.1 && c.cl < 4.0 * std::numbers::pi * 0.02 * 1.05,
          "positive camber must lift at zero incidence, below the inviscid value");
}

void TestReynoldsTrend() {
    const Solver::BoundaryLayerSectionModel model;
    const auto strip = FlatStrip();
    const double cdLow = model(strip, 2.0, 1e5, 0.0).cd;
    const double cdHigh = model(strip, 2.0, 2e6, 0.0).cd;
    std::cout << "cd(Re=1e5)=" << cdLow << "  cd(Re=2e6)=" << cdHigh << "\n";
    CHECK(cdLow > cdHigh, "drag must fall as Reynolds number rises");
}

void TestTranspirationDecambers() {
    // The same solver with the boundary layer silenced (one iteration, no
    // transpiration applied yet) is the inviscid lumped-vortex answer; the
    // converged coupled lift must sit below it.
    Solver::BoundaryLayerSectionModel coupled;
    Solver::BoundaryLayerSectionModel inviscid;
    inviscid.MaxIterations = 1; // first pass solves with zero transpiration
    const auto strip = FlatStrip();
    const double clCoupled = coupled(strip, 6.0, 3e5, 0.0).cl;
    const double clInviscid = inviscid(strip, 6.0, 3e5, 0.0).cl;
    std::cout << "alpha=6: cl inviscid pass=" << clInviscid << "  coupled=" << clCoupled << "\n";
    CHECK(clCoupled < clInviscid,
          "the displacement transpiration must decamber (reduce lift), never add it");
    CHECK(clCoupled > 0.5 * clInviscid, "the viscous decrement must stay a correction, not a collapse");
}

} // namespace

int main() {
    TestFlatPlate();
    TestCamberLifts();
    TestReynoldsTrend();
    TestTranspirationDecambers();

    if (failures == 0) { std::cout << "PASS: TestSectionBoundaryLayer\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestSectionBoundaryLayer\n";
    return 1;
}
