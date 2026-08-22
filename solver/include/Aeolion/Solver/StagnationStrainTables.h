// Solver/StagnationStrainTables.h
//
// The per-section stagnation-strain tables that seed the boundary-layer
// march in SectionBoundaryLayer.h: dU_e/ds at the attachment point, per
// chord and per unit oncoming speed, as a function of (eta, alpha).
//
// --- why this is a table and not a call ---------------------------------
// BoundaryLayerSectionModel is a CAMBER LINE. It has no thickness, so no
// stagnation region exists inside it, and a march starting at station 0 with
// theta = 0 discards the whole upstream run. The momentum thickness that run
// would have delivered is the Hiemenz value theta_0 = sqrt(0.075 nu / a), and
// `a` is a property of the thickness-resolved flow around the nose. It cannot
// be recovered from a camber line at any price.
//
// It CAN be computed, by the Hess-Smith solve in SectionPanelMethod.h -- but
// not per call. BoundaryLayerSectionModel iterates to convergence inside the
// viscous coupling, which itself iterates up to a thousand times per flight
// condition, so a section panel solve on that path would be paid tens of
// thousands of times per condition to produce a number that depends only on
// (eta, alpha). Hence a table, built once, exactly as BuildSeparationTables
// handles separation for the same reason.
//
// --- why alpha is signed here and not in the separation tables ----------
// The separation tables key on |alpha| from the zero-lift line, which is
// sound because that is the variable the section model is posed in. The
// stagnation strain is not that kind of quantity: it is a property of where
// the attachment point sits on a CAMBERED nose, and a cambered nose is not
// symmetric. At +5 and -5 degrees the stagnation point sits on opposite
// surfaces at different curvatures, so the strains differ. Tabulate signed.
#pragma once

#include "Aeolion/Geometry/AirfoilSection.h"
#include "Aeolion/Geometry/SectionContour.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Solver/SectionPanelMethod.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

namespace Aeolion::Solver {

/** One section's stagnation strain against incidence, ascending in alpha. */
struct StagnationStrainTable {
    std::vector<double> AlphaDeg; ///< Section incidence about its own chord line, signed.
    std::vector<double> Strain;   ///< dU_e/ds, per chord per unit oncoming speed.
};

/**
 * Range and step of the tabulation.
 *
 * The march this feeds only runs inside BlendEndDeg of the zero-lift line --
 * past that the section model returns its saturated analytic polar and never
 * marches at all. Twenty degrees each way covers that envelope about any
 * plausible zero-lift angle with room to spare, and the strain is smooth in
 * incidence, so a one-degree step is comfortably finer than it needs to be
 * for a quantity that enters as a square root.
 */
inline constexpr double StrainTableMaxAlphaDeg = 20.0;
inline constexpr double StrainTableStepDeg = 1.0;

/**
 * Build one table per section. Each entry is a Hess-Smith solve on the
 * section's own closed contour, so the cost is (sections x alphas) panel
 * solves, paid once.
 */
[[nodiscard]] inline std::vector<StagnationStrainTable> BuildStagnationStrainTables(
    const std::vector<Geometry::AirfoilSection>& sections,
    int pointsPerSurface = Geometry::DefaultContourPointsPerSurface) {
    std::vector<StagnationStrainTable> tables(sections.size());
    for (std::size_t i = 0; i < sections.size(); ++i) {
        const Geometry::SectionContour contour =
            Geometry::BuildSectionContour(sections[i], pointsPerSurface);
        if (!contour.Valid()) continue;
        for (double a = -StrainTableMaxAlphaDeg; a <= StrainTableMaxAlphaDeg + 1e-9;
             a += StrainTableStepDeg) {
            const SectionSolution solution = SolveSectionContour(contour, Math::DegToRad(a));
            tables[i].AlphaDeg.push_back(a);
            tables[i].Strain.push_back(solution.StagnationStrain);
        }
    }
    return tables;
}

/**
 * Nearest section by span fraction, linear interpolation in incidence,
 * CLAMPED at both ends rather than extrapolated.
 *
 * Clamping is the right edge rule here, and it is not the one the separation
 * tables use. A separation point must keep walking forward past the table
 * edge or the model is denied its plate limit; a stagnation strain must not,
 * because a linear extrapolation of it eventually goes NEGATIVE, and a
 * negative strain has no Hiemenz solution -- theta_0 = sqrt(0.075 nu / a)
 * would take the square root of a negative number. Beyond the table the
 * march is not running anyway, so the clamp costs nothing real.
 *
 * Returns 0 for an empty table, which the consumer reads as "no seed" and
 * which reproduces the unseeded march exactly.
 */
[[nodiscard]] inline std::function<double(double eta, double alphaDeg)>
MakeStagnationStrainFunction(std::vector<StagnationStrainTable> tables,
                             std::vector<double> etas) {
    return [tables = std::move(tables), etas = std::move(etas)](double eta, double aDeg) {
        if (tables.empty() || etas.empty()) return 0.0;
        std::size_t best = 0;
        double bestDist = 1e30;
        for (std::size_t i = 0; i < etas.size() && i < tables.size(); ++i) {
            const double d = std::fabs(etas[i] - eta);
            if (d < bestDist) {
                bestDist = d;
                best = i;
            }
        }
        const StagnationStrainTable& t = tables[best];
        const std::size_t n = t.AlphaDeg.size();
        if (n == 0) return 0.0;
        if (aDeg <= t.AlphaDeg.front()) return std::max(t.Strain.front(), 0.0);
        if (aDeg >= t.AlphaDeg.back()) return std::max(t.Strain.back(), 0.0);
        std::size_t hi = 1;
        while (hi + 1 < n && t.AlphaDeg[hi] < aDeg) ++hi;
        const double a0 = t.AlphaDeg[hi - 1], a1 = t.AlphaDeg[hi];
        const double frac = (a1 - a0 > 1e-12) ? (aDeg - a0) / (a1 - a0) : 0.0;
        return std::max(t.Strain[hi - 1] + frac * (t.Strain[hi] - t.Strain[hi - 1]), 0.0);
    };
}

} // namespace Aeolion::Solver
