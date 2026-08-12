// TestDiscreteVortexSection.cpp -- the LESP-modulated discrete-vortex
// section (Solver/DiscreteVortexSection.h), pinned against the classical
// results that bound an unsteady thin-airfoil wake model: the steady limit,
// Wagner's growth, Kelvin's theorem to machine precision, the LESP
// criterion's on/off behaviour, and the bluff-plate band at 90 degrees.

#include "Aeolion/Solver/DiscreteVortexSection.h"

#include <cmath>
#include <iostream>
#include <numbers>

using namespace Aeolion;
namespace S = Aeolion::Solver;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

constexpr double Pi = std::numbers::pi;

void TestSteadyLimitAndWagner() {
    // Attached flow: a huge LESP threshold suppresses the LEV entirely, and
    // the run must relax onto the steady thin-airfoil value from below,
    // through Wagner's growth.
    S::DiscreteVortexOptions options;
    options.LespCritical = 100.0;
    options.Duration = 30.0;
    options.KeepHistory = true;

    const double alphaDeg = 4.0;
    const auto run = S::SolveDiscreteVortexSection(alphaDeg, options);
    CHECK(run.Valid, "the attached run completes");
    CHECK(run.LevShed == 0, "no leading-edge vortex below the threshold");

    // The impulse loads are exact against Wagner, so the expectation over
    // the averaging window is NOT the steady value but the window mean of
    // Wagner's function (Jones' two-exponential fit), phi ~ 0.938 over
    // s = 15..30. Pinning against that catches both a lift deficit AND a
    // formulation that quietly rides above Wagner's growth.
    const double steady = 2.0 * Pi * std::sin(alphaDeg * Pi / 180.0);
    double wagnerMean = 0.0;
    int wagnerCount = 0;
    for (double s = 15.0; s < 30.0; s += 0.1) {
        wagnerMean += 1.0 - 0.165 * std::exp(-0.045 * 2.0 * s) - 0.335 * std::exp(-0.3 * 2.0 * s);
        ++wagnerCount;
    }
    wagnerMean /= wagnerCount;
    CHECK(std::fabs(run.MeanCl - wagnerMean * steady) < 0.04 * steady,
          "mean cl tracks the Wagner-window mean within 4%");
    CHECK(run.RmsCl < 0.06 * steady, "the attached wake is quiet to Wagner's own drift");

    // Wagner: the circulatory lift starts near HALF the steady value and
    // grows monotonically toward it. Sample at one and at fifteen chords.
    double clAt1 = 0.0, clAt15 = 0.0;
    for (const auto& s : run.History) {
        if (std::fabs(s.t - 1.0) < 0.5 * options.TimeStep) clAt1 = s.cl;
        if (std::fabs(s.t - 15.0) < 0.5 * options.TimeStep) clAt15 = s.cl;
    }
    CHECK(clAt1 > 0.45 * steady && clAt1 < 0.85 * steady,
          "cl at one chord of travel sits in Wagner's band");
    CHECK(clAt15 > clAt1, "the lift grows toward the steady value");
    CHECK(run.KelvinResidual < 1e-10, "Kelvin's theorem holds to machine precision");

    std::cout << "attached alpha=4: mean cl " << run.MeanCl << " (steady " << steady
              << "), cl(s=1)/steady " << clAt1 / steady << "\n";
}

void TestLespCriterionSwitches() {
    // The same incidence with a realistic threshold must shed from the
    // leading edge; A0 must be capped at the threshold while it does.
    S::DiscreteVortexOptions options;
    options.LespCritical = 0.2;
    options.Duration = 20.0;
    options.KeepHistory = true;

    const auto run = S::SolveDiscreteVortexSection(25.0, options);
    CHECK(run.Valid, "the separated run completes");
    CHECK(run.LevShed > 0, "alpha = 25 exceeds a 0.2 threshold and sheds LEVs");
    CHECK(run.KelvinResidual < 1e-10, "Kelvin holds with both edges shedding");

    double worstA0 = 0.0;
    for (const auto& s : run.History)
        if (s.LevActive) worstA0 = std::max(worstA0, std::fabs(s.A0));
    CHECK(worstA0 <= options.LespCritical + 1e-6,
          "A0 is capped at the critical LESP while the LEV sheds");

    // The shed flow must carry unsteadiness the quasi-steady tiers cannot.
    CHECK(run.RmsCl > 0.02, "deep separation fluctuates");
}

void TestDeepStallBand() {
    // The normal plate, pinned to the TWO-dimensional expectation: a 2-D
    // vortex street has no spanwise breakup and famously overpredicts the
    // measured plate drag (~2) by 60-70% -- classical 2-D simulations sit
    // near 3.3. A mean inside [2.2, 4.0] says the street is doing what a
    // 2-D street does; a mean near the measured 2 would actually be
    // SUSPICIOUS here (see the header's dimension-bias note).
    S::DiscreteVortexOptions options;
    options.LespCritical = 0.2;
    options.Duration = 25.0;

    const auto run = S::SolveDiscreteVortexSection(90.0, options);
    CHECK(run.Valid, "the normal-plate run completes");
    CHECK(run.MeanCd > 2.2 && run.MeanCd < 4.0,
          "mean cd(90) lands in the 2-D street's band");
    CHECK(std::fabs(run.MeanCl) < 0.35, "mean cl(90) is near zero by symmetry");
    CHECK(run.RmsCd > 0.05, "the plate wake sheds");
    std::cout << "plate alpha=90: mean cd " << run.MeanCd << " rms cd " << run.RmsCd
              << " mean cl " << run.MeanCl << " (LEVs " << run.LevShed << ", TEVs "
              << run.TevShed << ")\n";
}

void TestSymmetry() {
    // Odd symmetry of the mean lift under alpha -> -alpha. Shedding is
    // chaotic, so the tolerance is loose but the sign is not negotiable.
    S::DiscreteVortexOptions options;
    options.LespCritical = 0.2;
    options.Duration = 20.0;
    const auto plus = S::SolveDiscreteVortexSection(30.0, options);
    const auto minus = S::SolveDiscreteVortexSection(-30.0, options);
    CHECK(plus.Valid && minus.Valid, "both signs complete");
    CHECK(plus.MeanCl > 0.0 && minus.MeanCl < 0.0, "mean lift carries the sign of alpha");
    // The shedding is chaotic and twenty convective times hold only a few
    // cycles, so the mirror tolerances are those of a short average, not
    // of the formulation.
    CHECK(std::fabs(plus.MeanCl + minus.MeanCl) < 0.2 * std::fabs(plus.MeanCl) + 0.08,
          "mean lift is odd in alpha to shedding scatter");
    CHECK(std::fabs(plus.MeanCd - minus.MeanCd) < 0.25 * plus.MeanCd + 0.1,
          "mean drag is even in alpha to shedding scatter");
    std::cout << "symmetry: cl " << plus.MeanCl << " / " << minus.MeanCl << ", cd "
              << plus.MeanCd << " / " << minus.MeanCd << "\n";
}

} // namespace

int main() {
    TestSteadyLimitAndWagner();
    TestLespCriterionSwitches();
    TestDeepStallBand();
    TestSymmetry();

    if (failures) {
        std::cerr << failures << " check(s) failed in TestDiscreteVortexSection\n";
        return 1;
    }
    std::cout << "TestDiscreteVortexSection: all checks passed\n";
    return 0;
}
