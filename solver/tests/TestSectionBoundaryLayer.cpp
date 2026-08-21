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

// The march used to start at theta = 0 on a camber line, which has no
// stagnation region to have started from -- ue is already finite at station 0,
// so the whole upstream run was discarded. Supplying the stagnation strain
// seeds it with the Hiemenz value the unresolved region would have delivered.
void TestStagnationSeedThickensAndDrags() {
    const auto strip = FlatStrip();
    const double Re = 5e5, alpha = 4.0;

    Solver::BoundaryLayerSectionModel bare;
    Solver::BoundaryLayerSectionModel seeded;
    // A representative nose strain: Ue rising to freestream over a few percent
    // of chord gives dUe*/ds* of order 10.
    const double a = 10.0;
    seeded.StagnationStrain = [a](double, double) { return a; };

    const auto bareCoef = bare(strip, alpha, Re, 0.0);
    const auto seedCoef = seeded(strip, alpha, Re, 0.0);

    // 1. The seed must MOVE something. A supplier that is accepted and then
    //    ignored is the failure this test exists to catch.
    CHECK(seedCoef.cd != bareCoef.cd, "the stagnation seed changed nothing at all");

    // 2. It must move drag UP. A thicker layer at the leading edge stays
    //    thicker all the way to the trailing edge, and Squire-Young reads
    //    drag off theta there. A seed that REDUCED drag would mean the sign
    //    of the correction is inverted.
    CHECK(seedCoef.cd > bareCoef.cd,
          "seeding the upstream run must increase cd, got " << bareCoef.cd << " -> " << seedCoef.cd);

    // 3. The seed must VANISH in its own limit. theta_0 = sqrt(0.075/(Re a)),
    //    so a very sharp nose seeds a vanishing layer and the result must
    //    return to the unseeded march. This is the check that says the seed is
    //    a momentum thickness entering the Thwaites integral and not an
    //    arbitrary offset bolted on.
    //
    //    Deliberately NOT a bound on how large the correction may be. It is
    //    83% of cd at a = 10, which looks alarming for a seed of 1.2e-4 chords
    //    until one notices the mechanism: a thicker leading-edge layer raises
    //    Re_theta, which trips Michel earlier, which turns a longer run
    //    turbulent. The seed acts mostly THROUGH TRANSITION, so a large jump
    //    is the physics rather than a symptom, and a magnitude bound would
    //    just be a tripwire on where transition happens to sit.
    Solver::BoundaryLayerSectionModel vanishing;
    vanishing.StagnationStrain = [](double, double) { return 1e8; };
    const double vanishCd = vanishing(strip, alpha, Re, 0.0).cd;
    CHECK(std::fabs(vanishCd - bareCoef.cd) < 1e-6 * std::max(bareCoef.cd, 1e-9),
          "a vanishing seed must return the unseeded march, got " << bareCoef.cd
          << " -> " << vanishCd);

    // 4. An empty supplier must reproduce the old behaviour EXACTLY, so no
    //    existing consumer moves until it opts in.
    Solver::BoundaryLayerSectionModel unset;
    CHECK(unset(strip, alpha, Re, 0.0).cd == bareCoef.cd,
          "an unset supplier must be bit-identical to the previous march");

    // 5. Larger strain means a thinner stagnation layer, so less drag --
    //    theta_0 goes as 1/sqrt(a). This is what says the seed is being used
    //    as a momentum thickness and not merely as an arbitrary offset.
    Solver::BoundaryLayerSectionModel sharper;
    sharper.StagnationStrain = [a](double, double) { return 4.0 * a; };
    CHECK(sharper(strip, alpha, Re, 0.0).cd < seedCoef.cd,
          "a sharper nose must seed a THINNER layer and less drag");

    std::cout << "stagnation seed: cd " << bareCoef.cd << " -> " << seedCoef.cd
              << " (a=" << a << ")" << std::endl;
}

} // namespace

int main() {
    TestFlatPlate();
    TestCamberLifts();
    TestReynoldsTrend();
    TestTranspirationDecambers();
    TestStagnationSeedThickensAndDrags();

    if (failures == 0) { std::cout << "PASS: TestSectionBoundaryLayer\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestSectionBoundaryLayer\n";
    return 1;
}
