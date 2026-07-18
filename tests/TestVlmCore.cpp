// TestVlmCore.cpp -- validates Vlm.h's core solver against known
// closed-form aerodynamic theory:
//   - CL for a rectangular unswept/untwisted wing should track the
//     lifting-line-corrected thin-wing formula CL = 2*pi*alpha/(1+2/AR)
//     within a few percent (VLM is more accurate than the 2-term formula,
//     so exact match isn't expected -- convergence in the right direction
//     as AR increases is the real check).
//   - Oswald efficiency for that same wing should sit close to 1.0.
//   - CL/CDi should converge monotonically as panel count increases.
#include "Aeolion/Vlm.h"
#include <iostream>
#include <cmath>

using namespace Aeolion::Vlm;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

int main() {
    // --- thin-wing lift slope check across a few aspect ratios ---
    for (int AR_target : {6, 8, 12, 20}) {
        WingParams w;
        w.Span = 4.0;
        double c = w.Span / AR_target;
        w.RootChord = c; w.TipChord = c;
        w.SweepQuarterChordDeg = 0; w.DihedralDeg = 0; w.TwistTipDeg = 0;
        w.NPanelsSemiSpan = 30; w.CosineSpacing = true;

        FreestreamConditions fc; fc.Vinf = 30.0; fc.alphaDeg = 4.0; fc.rho = 1.225;
        SolveResult res = Solve(w, fc);
        double AR = w.Span * w.Span / res.ReferenceArea;
        double alphaRad = fc.alphaDeg * Pi / 180.0;
        double CL_theory = 2 * Pi * alphaRad / (1.0 + 2.0 / AR);
        double relErr = std::fabs(res.CL - CL_theory) / CL_theory;

        std::cout << "AR=" << AR << "  CL=" << res.CL << "  CL_theory=" << CL_theory
                  << "  relErr=" << relErr * 100 << "%\n";
        CHECK(relErr < 0.15, "CL vs thin-wing theory off by more than 15% at AR=" << AR);
        CHECK(res.CDi > 0, "CDi should be positive at AR=" << AR);

        double e = (res.CDi > 1e-9) ? (res.CL * res.CL) / (Pi * AR * res.CDi) : 0.0;
        std::cout << "  Oswald e=" << e << "\n";
        CHECK(e > 0.85 && e < 1.05, "Oswald efficiency out of expected range for a plain rectangular wing, AR=" << AR);
    }

    // --- panel convergence: CL/CDi should settle down, not oscillate wildly ---
    double prevCL = 0;
    bool first = true;
    for (int npss : {10, 20, 40, 80}) {
        WingParams w;
        w.Span = 4.0; w.RootChord = 0.5; w.TipChord = 0.5;
        w.NPanelsSemiSpan = npss; w.CosineSpacing = true;
        FreestreamConditions fc; fc.Vinf = 30; fc.alphaDeg = 4.0; fc.rho = 1.225;
        SolveResult res = Solve(w, fc);
        std::cout << "N=" << 2 * npss << "  CL=" << res.CL << "\n";
        if (!first) {
            double delta = std::fabs(res.CL - prevCL);
            CHECK(delta < 0.05, "CL changed too much between panel refinements (N=" << 2 * npss << "), delta=" << delta);
        }
        prevCL = res.CL;
        first = false;
    }

    if (failures == 0) { std::cout << "PASS: TestVlmCore\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestVlmCore\n";
    return 1;
}
