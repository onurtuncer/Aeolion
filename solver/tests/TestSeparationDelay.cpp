// Regression test for the contract the B1 separation measurement rests on.
//
// B1 answers "how far does the aft fan delay separation" by reading the
// separation point off the coupled solve directly, instead of differencing
// two limit-cycle means. That measurement is only meaningful because
// MakeSeparationFunction has four properties, none of which the separation
// tables themselves enforce -- they are functions of (eta, alpha) alone and
// have never heard of the fan. Every property below is one whose loss would
// corrupt the measured map SILENTLY, without any magnitude looking wrong.
#include "Aeolion/Solver/SeparationTables.h"
#include <cassert>
#include <cmath>
#include <cstdio>

int main() {
    using namespace Aeolion::Solver;

    // A table shaped like the real ones: attached near the trailing edge at
    // low incidence, walking forward, reaching the leading edge by ~30 deg.
    StripSeparationTable t;
    t.AlphaDeg = {0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0};
    t.Psi      = {0.99, 0.97, 0.94, 0.88, 0.76, 0.55, 0.28, 0.06};
    const auto f = MakeSeparationFunction({t, t}, {0.25, 0.75});

    // 1. MONOTONE in incidence. This is what makes "the fan lowers local
    //    incidence" and "the fan moves separation aft" the SAME statement.
    //    Lose it and a rise in f stops meaning a delay.
    double prev = f(0.25, 0.0);
    for (double a = 0.5; a <= 40.0; a += 0.5) {
        const double v = f(0.25, a);
        assert(v <= prev + 1e-12 && "separation point moved AFT with rising incidence");
        prev = v;
    }

    // 2. SIGN CONVENTION: f = 1 attached, f = 0 separated at the LE, so a fan
    //    that delays separation RAISES f. A silent flip here would inverting
    //    the paper's headline table without changing one magnitude -- exactly
    //    the failure a magnitude check cannot catch.
    assert(f(0.25, 13.0) > f(0.25, 14.0) && "lowering incidence must RAISE f");

    // 3. CLAMPED, and extrapolating to the plate limit rather than freezing.
    //    The measured map reads exactly 0 on the worst strip from alpha 20 up,
    //    at every thrust setting including Tc = 8; that saturation is why
    //    differencing two saturated values is not a measurement, and why B1
    //    reports the span mean past that point. Freezing f above zero instead
    //    would deny the model its plate limit and manufacture a fan effect
    //    where the metric has simply run out of range.
    assert(f(0.25, 60.0) == 0.0 && "f must reach the plate limit, not freeze");
    for (double a = 0.0; a <= 120.0; a += 3.0) {
        const double v = f(0.25, a);
        assert(v >= 0.0 && v <= 1.0 && "f left [0,1]");
    }

    // 4. SYMMETRIC in incidence -- the table is keyed on |alpha| from the
    //    zero-lift line. The map is swept to negative alpha, so a sign-
    //    sensitive lookup would corrupt those rows only.
    assert(std::fabs(f(0.25, -14.0) - f(0.25, 14.0)) < 1e-12 && "lookup not symmetric");

    std::printf("TestSeparationDelay: monotone, sign, plate limit, symmetry -- OK\n");
    return 0;
}
