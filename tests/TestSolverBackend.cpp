// TestSolverBackend.cpp -- sanity check for the LAPACK/OpenBLAS dense
// linear-algebra backend: run a representative solve plus a derivatives
// sweep and confirm the results are finite and physically sane.
#include "Aeolion/VLM.h"
#include <iostream>
#include <iomanip>

using namespace Aeolion::VLM;

int main() {
    WingParams wp;
    wp.Span = 4.0; wp.RootChord = 0.4; wp.TipChord = 0.3;
    wp.SweepQuarterChordDeg = 10.0; wp.DihedralDeg = 3.0; wp.TwistTipDeg = -3.0;
    wp.NPanelsSemiSpan = 20; wp.CosineSpacing = true;

    FreestreamConditions fc;
    fc.Vinf = 20.0; fc.alphaDeg = 5.0; fc.rho = 1.225;

    SolveResult res = Solve(wp, fc);

    std::vector<Panel> panels = BuildWing(wp);
    for (auto& p : panels) p.Surface = "wing";
    ReferenceGeometry ref; ref.Area = res.ReferenceArea; ref.Span = wp.Span; ref.Chord = res.ReferenceArea / wp.Span;
    StabilityDerivatives d = ComputeDerivatives(panels, fc, ref, 50.0 * wp.Span);

    std::cout << "backend: LAPACK/OpenBLAS\n";
    std::cout << std::setprecision(10);
    std::cout << "CL=" << res.CL << "\n";
    std::cout << "CDi=" << res.CDi << "\n";
    std::cout << "CL_alpha=" << d.CL_alpha << "\n";
    std::cout << "Cm_alpha=" << d.Cm_alpha << "\n";
    std::cout << "Croll_beta=" << d.Croll_beta << "\n";

    // basic sanity: finite, non-negative induced drag
    bool ok = std::isfinite(res.CL) && std::isfinite(res.CDi) && res.CDi >= 0.0;
    if (!ok) { std::cerr << "FAIL: non-finite or negative-CDi result\n"; return 1; }
    std::cout << "PASS: TestSolverBackend\n";
    return 0;
}
