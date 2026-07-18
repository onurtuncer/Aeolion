// TestBoundaryLayer.cpp -- validates BoundaryLayer.h on a zero-
// pressure-gradient flat plate against known closed-form results:
//   - laminar Thwaites should match the exact Blasius solution (Thwaites
//     is calibrated to reproduce zero-pressure-gradient flow exactly).
//   - turbulent Head's method should converge to H~1.3-1.4 (the known
//     equilibrium shape factor) and agree with standard flat-plate Cf
//     correlations to within normal cross-method scatter (~15-20%).
#include "Aeolion/BoundaryLayer.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace Aeolion;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

int main() {
    double Vinf = 20.0, nu = 1.5e-5, L = 1.0;
    int N = 2000;
    std::vector<double> s(N + 1), Ue(N + 1);
    for (int i = 0; i <= N; ++i) { s[i] = L * i / N; Ue[i] = Vinf; }
    double Re_L = Vinf * L / nu;

    auto lam = BoundaryLayer::AnalyzeSurface(s, Ue, Vinf, nu, L, false);
    double cd_blasius = 1.328 / std::sqrt(Re_L);
    double lamRelErr = std::fabs(lam.cd - cd_blasius) / cd_blasius;
    std::cout << "Re_L=" << Re_L << "  laminar cd=" << lam.cd << "  Blasius=" << cd_blasius
              << "  relErr=" << lamRelErr * 100 << "%\n";
    CHECK(lamRelErr < 0.03, "Thwaites laminar result should match Blasius within ~3%");
    CHECK(lam.sTransition < 0, "zero-pressure-gradient flat plate at this Re should stay laminar to the end");

    auto turb = BoundaryLayer::AnalyzeSurface(s, Ue, Vinf, nu, L, true);
    double cd_schlichting = 0.455 / std::pow(std::log10(Re_L), 2.58);
    double cd_prandtl = 0.074 / std::pow(Re_L, 0.2);
    double errSch = std::fabs(turb.cd - cd_schlichting) / cd_schlichting;
    double errPr = std::fabs(turb.cd - cd_prandtl) / cd_prandtl;
    std::cout << "turbulent cd=" << turb.cd << "  Schlichting=" << cd_schlichting << " (err " << errSch * 100 << "%)"
              << "  Prandtl-1/7=" << cd_prandtl << " (err " << errPr * 100 << "%)  H_te=" << turb.H_te << "\n";
    CHECK(errSch < 0.25, "turbulent cd should be within ~25% of the Schlichting correlation (independent empirical closures, some scatter expected)");
    CHECK(errPr < 0.25, "turbulent cd should be within ~25% of the Prandtl 1/7-power correlation");
    CHECK(turb.H_te > 1.2 && turb.H_te < 1.6, "turbulent flat-plate shape factor should sit near the known ~1.3-1.4 equilibrium value");

    if (failures == 0) { std::cout << "PASS: TestBoundaryLayer\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestBoundaryLayer\n";
    return 1;
}
