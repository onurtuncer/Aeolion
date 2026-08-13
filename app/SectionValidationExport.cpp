// SectionValidationExport.cpp -- section-level validation of the anchored
// post-stall model (Solver/PostStallSection.h) against the Sheldahl-Klimas
// NACA 0015 table at Re = 3.6e5: the declared acceptance test of the
// post-separation study's quasi-steady tier (Part II, "Validation and the
// remaining tiers").
//
// The comparison is END TO END through the same machinery the airframe
// solve uses -- nothing is fed from the data being compared against:
//
//   1. The NACA 0015 half-thickness law (closed trailing edge) is fitted
//      with the schema's own CST parameterization (N1 = 1/2, N2 = 1,
//      Bernstein order 7) by least squares; the fit residual and the
//      recovered leading-edge radius are reported so the geometry error is
//      on the record (~2e-4 of chord; r_LE within a percent of 0.0248c).
//   2. The suction-side separation point f(alpha) comes from the SAME
//      pipeline as the airframe's: Hess-Smith on the contour
//      (SolveSectionContour), then the attachment-point march
//      (MarchSurfaceRun) at the table's Reynolds number.
//   3. PostStallSectionModel(f) produces the polar to 90 deg. AspectRatio
//      is set to the Viterna fit's AR = 50 edge -- the two-dimensional
//      limit of the ceiling, CdMax = 2.01 (see PostStallSection.h).
//
// Output: JSON polar (alpha, f, cl, cd, cm) plus the fit diagnostics; the
// figure/prose comparison against the Sandia-distributed table lives in
// papers/journal-of-aircraft-poststall/figures/render-validation-figure.py.
//
// Provenance caveat, stated where the numbers are made: Sheldahl & Klimas
// (SAND80-2114) measured this section through stall; the deep post-stall
// range of the distributed 360-degree tables is their synthesis. The deep
// range of the comparison is therefore model-against-community-standard,
// not model-against-raw-measurement, and the paper says so.

#include "Aeolion/Geometry/AirfoilSection.h"
#include "Aeolion/Geometry/SectionContour.h"
#include "Aeolion/Solver/AttachmentBoundaryLayer.h"
#include "Aeolion/Solver/DiscreteVortexSection.h"
#include "Aeolion/Solver/PostStallSection.h"
#include "Aeolion/Solver/SectionPanelMethod.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace Aeolion;
namespace S = Aeolion::Solver;

namespace {

constexpr double DefaultReynolds = 3.6e5; // the Sheldahl-Klimas block compared against
constexpr int CstOrder = 7;               // Bernstein order of the thickness fit
constexpr int FitNodes = 60;              // cosine-spaced psi nodes for the least squares
constexpr double SeparationTableMaxDeg = 40.0;
constexpr double SeparationTableStepDeg = 1.0;
constexpr double PolarMaxDeg = 90.0;
constexpr double PolarStepDeg = 1.0;

// NACA 00xx half-thickness, closed-trailing-edge form (the -0.1036
// variant), t = maximum thickness as a chord fraction.
double NacaHalfThickness(double t, double psi) {
    const double r = std::sqrt(std::max(psi, 0.0));
    return 5.0 * t *
           (0.2969 * r - 0.1260 * psi - 0.3516 * psi * psi + 0.2843 * psi * psi * psi -
            0.1036 * psi * psi * psi * psi);
}

double Bernstein(int order, int i, double x) {
    double comb = 1.0;
    for (int k = 1; k <= i; ++k) comb *= static_cast<double>(order - k + 1) / k;
    return comb * std::pow(x, i) * std::pow(1.0 - x, order - i);
}

// Least-squares CST fit of the half-thickness law: z = psi^N1 (1-psi)^N2 *
// sum_i A_i B_i(psi). Solved via the (small, well-conditioned at order 7)
// normal equations.
std::vector<double> FitCst(double thickness, double& maxResidual) {
    const int n = CstOrder + 1;
    std::vector<double> normal(n * n, 0.0), rhs(n, 0.0);
    std::vector<double> psis(FitNodes);
    for (int k = 0; k < FitNodes; ++k)
        psis[k] = 0.5 * (1.0 - std::cos(std::numbers::pi * (k + 0.5) / FitNodes));

    for (double psi : psis) {
        const double classFn = std::pow(psi, Geometry::CstN1) * std::pow(1.0 - psi, Geometry::CstN2);
        const double target = NacaHalfThickness(thickness, psi);
        for (int i = 0; i < n; ++i) {
            const double bi = classFn * Bernstein(CstOrder, i, psi);
            rhs[i] += bi * target;
            for (int j = 0; j < n; ++j)
                normal[i * n + j] += bi * classFn * Bernstein(CstOrder, j, psi);
        }
    }
    // Gaussian elimination with partial pivoting.
    std::vector<double> a(normal);
    std::vector<double> x(rhs);
    for (int col = 0; col < n; ++col) {
        int best = col;
        for (int row = col + 1; row < n; ++row)
            if (std::fabs(a[row * n + col]) > std::fabs(a[best * n + col])) best = row;
        for (int k = 0; k < n; ++k) std::swap(a[col * n + k], a[best * n + k]);
        std::swap(x[col], x[best]);
        for (int row = col + 1; row < n; ++row) {
            const double f = a[row * n + col] / a[col * n + col];
            for (int k = col; k < n; ++k) a[row * n + k] -= f * a[col * n + k];
            x[row] -= f * x[col];
        }
    }
    for (int row = n; row-- > 0;) {
        for (int k = row + 1; k < n; ++k) x[row] -= a[row * n + k] * x[k];
        x[row] /= a[row * n + row];
    }

    maxResidual = 0.0;
    for (double psi : psis) {
        const double classFn = std::pow(psi, Geometry::CstN1) * std::pow(1.0 - psi, Geometry::CstN2);
        double z = 0.0;
        for (int i = 0; i < n; ++i) z += x[i] * classFn * Bernstein(CstOrder, i, psi);
        maxResidual = std::max(maxResidual, std::fabs(z - NacaHalfThickness(thickness, psi)));
    }
    return x;
}

} // namespace

int main(int argc, char** argv) {
    const std::string outPath = (argc > 1) ? argv[1] : "section-validation.json";
    const double Re = (argc > 2) ? std::atof(argv[2]) : DefaultReynolds;

    // --- the section, through the schema's own parameterization ------------
    double fitResidual = 0.0;
    Geometry::AirfoilSection section;
    section.Eta = 0.0;
    section.CoefficientsUpper = FitCst(0.15, fitResidual);
    section.CoefficientsLower = section.CoefficientsUpper;
    for (double& c : section.CoefficientsLower) c = -c;

    const Geometry::SectionContour contour = Geometry::BuildSectionContour(section);
    if (!contour.Valid()) {
        std::cerr << "contour construction failed\n";
        return 1;
    }
    std::cout << "NACA 0015 CST fit: max residual " << fitResidual << " c, r_LE/c = "
              << contour.LeadingEdgeRadius << " (exact 0.0248)\n";

    // --- f(alpha) from the same pipeline as the airframe's ------------------
    std::vector<double> tabAlpha, tabPsi;
    for (double aDeg = 0.0; aDeg <= SeparationTableMaxDeg + 1e-9; aDeg += SeparationTableStepDeg) {
        const S::SectionSolution sol = S::SolveSectionContour(contour, Math::DegToRad(aDeg));
        if (!sol.Valid || !sol.StagnationFound) continue;
        const S::SurfaceMarch march = S::MarchSurfaceRun(sol.UpperRun(), Re);
        if (!march.Valid) continue;
        tabAlpha.push_back(aDeg);
        tabPsi.push_back(march.SeparationPsi);
    }
    if (tabAlpha.size() < 4) {
        std::cerr << "separation table too sparse (" << tabAlpha.size() << " points)\n";
        return 1;
    }

    S::PostStallSectionModel model;
    model.AspectRatio = S::ViternaMaxAspectRatio; // the fit's 2-D edge, CdMax = 2.01
    model.SeparationPoint = [tabAlpha, tabPsi](double, double aDeg) {
        const double a = std::fabs(aDeg);
        if (a <= tabAlpha.front()) return std::clamp(tabPsi.front(), 0.0, 1.0);
        std::size_t hi = 1;
        while (hi + 1 < tabAlpha.size() && tabAlpha[hi] < a) ++hi;
        const double a0 = tabAlpha[hi - 1], a1 = tabAlpha[hi];
        const double frac = (a1 - a0 > 1e-9) ? (a - a0) / (a1 - a0) : 0.0; // >1 extrapolates
        return std::clamp(tabPsi[hi - 1] + frac * (tabPsi[hi] - tabPsi[hi - 1]), 0.0, 1.0);
    };

    S::StripSection strip;
    strip.Eta = 0.0;
    strip.Alpha0Deg = 0.0;
    strip.Chord = 1.0;
    strip.Width = 1.0;

    const double stallDeg = model.StallAngleDeg(strip.Eta);
    std::cout << "emergent stall: " << stallDeg << " deg from zero lift\n";

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot open " << outPath << "\n";
        return 1;
    }
    out.precision(9);
    out << "{\n\"meta\":{\"section\":\"NACA0015\",\"Re\":" << Re << ",\"fitResidual\":"
        << fitResidual << ",\"leadingEdgeRadius\":" << contour.LeadingEdgeRadius
        << ",\"aspectRatio\":" << model.AspectRatio << ",\"cdMax\":"
        << S::ViternaCdMax(model.AspectRatio) << ",\"stallDeg\":" << stallDeg
        << ",\"clAlphaPerRad\":" << model.ClAlphaPerRad << "},\n\"polar\":[\n";
    bool first = true;
    for (double aDeg = 0.0; aDeg <= PolarMaxDeg + 1e-9; aDeg += PolarStepDeg) {
        const S::SectionCoefficients c = model(strip, aDeg, Re, 0.0);
        const double f = model.SeparationPoint(strip.Eta, aDeg);
        if (!first) out << ",\n";
        first = false;
        out << " [" << aDeg << ',' << f << ',' << c.cl << ',' << c.cd << ',' << c.cm << ']';
    }
    out << "\n],\n";

    // --- the 2-D unsteady cross-check at the same attitudes ---------------------
    // Mean and RMS loads of the LESP-modulated discrete-vortex section
    // (tier 3's two-dimensional half) every ten degrees. Flat plate: the
    // cross-check's question is posed at deep incidence, where camber is
    // secondary and its 2-D coherence bias is the stated caveat.
    out << "\"dvm\":[\n";
    bool firstDvm = true;
    for (double aDeg = 10.0; aDeg <= 90.0 + 1e-9; aDeg += 10.0) {
        const S::DiscreteVortexResult run = S::SolveDiscreteVortexSection(aDeg);
        if (!run.Valid) continue;
        if (!firstDvm) out << ",\n";
        firstDvm = false;
        out << " [" << aDeg << ',' << run.MeanCl << ',' << run.MeanCd << ',' << run.RmsCl << ','
            << run.RmsCd << ',' << run.LevShed << ']';
        std::cout << "dvm alpha=" << aDeg << ": mean cl " << run.MeanCl << " cd " << run.MeanCd
                  << "  rms cl " << run.RmsCl << " cd " << run.RmsCd << "\n";
    }
    out << "\n]\n}\n";
    std::cout << "wrote " << outPath << " (" << tabAlpha.size() << " separation-table points)\n";
    return 0;
}
