// TestPanelMethod.cpp -- validates AirfoilPanel.h against exact
// potential-flow results:
//   - source-panel-only solve around a circle should exactly reproduce
//     Cp(theta) = 1 - 4*sin(theta)^2 (textbook closed form).
//   - a symmetric NACA section at alpha=0 must give Cl=0 exactly.
//   - Cl_alpha for a thin symmetric section should sit near 2*pi
//     (thin-airfoil theory), converging as panel count increases.
#include "Aeolion/AirfoilPanel.h"
#include <iostream>
#include <cmath>

using namespace Aeolion::Airfoil;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

static double CylinderCpError(int N) {
    std::vector<Vec2> loop;
    for (int i = 0; i <= N; ++i) {
        double th = -2 * Pi * i / N; // clockwise, matching BuildPanels' outward-normal convention
        loop.push_back({std::cos(th), std::sin(th)});
    }
    auto panels = BuildPanels(loop);
    int n = (int)panels.size();

    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    std::vector<double> b(n, 0.0);
    Vec2 Vinf{1.0, 0.0};
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            auto infl = ComputePanelInfluence(panels[i].Mid, panels[j]);
            A[i][j] = Dot2(infl.FromSource, panels[i].Normal);
        }
        b[i] = -Dot2(Vinf, panels[i].Normal);
    }
    // simple Gaussian elimination (source-only system, small N, fine for a test)
    for (int col = 0; col < n; ++col) {
        int piv = col; double best = std::fabs(A[col][col]);
        for (int r = col + 1; r < n; ++r) if (std::fabs(A[r][col]) > best) { best = std::fabs(A[r][col]); piv = r; }
        if (piv != col) { std::swap(A[piv], A[col]); std::swap(b[piv], b[col]); }
        for (int r = col + 1; r < n; ++r) {
            double f = A[r][col] / A[col][col];
            if (f == 0.0) continue;
            for (int c = col; c < n; ++c) A[r][c] -= f * A[col][c];
            b[r] -= f * b[col];
        }
    }
    std::vector<double> sigma(n, 0.0);
    for (int r = n - 1; r >= 0; --r) {
        double s = b[r];
        for (int c = r + 1; c < n; ++c) s -= A[r][c] * sigma[c];
        sigma[r] = (std::fabs(A[r][r]) > 1e-13) ? s / A[r][r] : 0.0;
    }

    double maxErr = 0;
    for (int i = 0; i < n; ++i) {
        Vec2 v = Vinf;
        for (int j = 0; j < n; ++j) {
            auto infl = ComputePanelInfluence(panels[i].Mid, panels[j]);
            v = v + infl.FromSource * sigma[j];
        }
        double vt = Dot2(v, panels[i].Tangent);
        double cp = 1.0 - vt * vt;
        double th = std::atan2(panels[i].Mid.y, panels[i].Mid.x);
        double cpExact = 1.0 - 4.0 * std::sin(th) * std::sin(th);
        maxErr = std::max(maxErr, std::fabs(cp - cpExact));
    }
    return maxErr;
}

int main() {
    double cylErr = CylinderCpError(120);
    std::cout << "cylinder max |Cp - Cp_exact| = " << cylErr << "\n";
    CHECK(cylErr < 1e-6, "source-panel cylinder solution should match the exact potential-flow result almost exactly");

    auto loop0012 = Naca4("0012", 160);
    auto panels0012 = BuildPanels(loop0012);
    SolveResult symZero = Solve(panels0012, 0.0);
    std::cout << "NACA0012 alpha=0: Cl=" << symZero.Cl << "\n";
    CHECK(std::fabs(symZero.Cl) < 1e-9, "symmetric airfoil at alpha=0 must give exactly zero lift");

    double a1 = 2.0 * Pi / 180.0, a2 = -2.0 * Pi / 180.0;
    SolveResult r1 = Solve(panels0012, a1), r2 = Solve(panels0012, a2);
    double clAlpha = (r1.Cl - r2.Cl) / (a1 - a2);
    std::cout << "NACA0012 Cl_alpha=" << clAlpha << " /rad (thin-airfoil theory: " << 2 * Pi << ")\n";
    CHECK(clAlpha > 2 * Pi * 0.95 && clAlpha < 2 * Pi * 1.25,
          "Cl_alpha should be close to (and typically slightly above) 2*pi for a thin-ish symmetric section");

    auto loop2412 = Naca4("2412", 160);
    auto panels2412 = BuildPanels(loop2412);
    SolveResult cambered = Solve(panels2412, 0.0);
    std::cout << "NACA2412 alpha=0: Cl=" << cambered.Cl << "\n";
    CHECK(cambered.Cl > 0.15 && cambered.Cl < 0.4, "cambered section should show positive lift at alpha=0, roughly 0.2-0.3");

    if (failures == 0) { std::cout << "PASS: TestPanelMethod\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestPanelMethod\n";
    return 1;
}
