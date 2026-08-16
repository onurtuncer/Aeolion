// Solver/SeparationTables.h
//
// The per-strip separation tables the anchored post-stall section model
// (PostStallSection.h) is driven by: f(eta, alpha_local) -- the
// suction-side separation point at a strip's LOCAL incidence measured
// from its own zero-lift line.
//
// The tables are built ONCE from an inviscid alpha sweep at beta = 0.
// Two things justify that, and both are measured rather than assumed:
// the Phase-0 map established that sideslip to 30 degrees only rescales
// the loads through cos(beta), so a one-dimensional table suffices; and
// keying by LOCAL incidence rather than by the aircraft attitude the
// sweep happened to run at is what lets one table serve every attitude,
// since local incidence is the variable the section model is posed in.
//
// This was driver-local code in app/PostStallSweepExport.cpp until a
// second consumer (the DAVE-ML aero map) needed the same tables. It is
// moved here VERBATIM -- same constants, same interpolation, same
// extrapolation rule -- so the post-stall study's numbers are unchanged
// by the extraction.
#pragma once

#include "Aeolion/Geometry/AirfoilSection.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Solver/AttachmentBoundaryLayer.h" // StationSeparation, SurveySeparation
#include "Aeolion/Solver/AttachmentLine.h"
#include "Aeolion/Solver/PostStallSection.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace Aeolion::Solver {

/** One strip's separation-point table, ascending in local incidence. */
struct StripSeparationTable {
    std::vector<double> AlphaDeg; ///< Local incidence from zero lift, ascending.
    std::vector<double> Psi;      ///< Suction-side separation point x/c at that incidence.
};

inline constexpr double SeparationTableMaxAlphaDeg = 32.0;
inline constexpr double SeparationTableStepDeg = 2.0;

/**
 * Build one table per strip from an inviscid alpha sweep on an already
 * prepared system. `fc` supplies the flight condition; its alpha and beta
 * are overwritten by the sweep.
 */
[[nodiscard]] inline std::vector<StripSeparationTable> BuildSeparationTables(
    const PreparedSystem& prepared, const std::vector<Panel>& wing,
    const std::vector<StripSection>& strips,
    const std::vector<Geometry::AirfoilSection>& sections, FreestreamConditions fc,
    const ReferenceGeometry& ref) {
    std::vector<StripSeparationTable> tables(strips.size());
    fc.betaDeg = 0.0;
    for (double alpha = 0.0; alpha <= SeparationTableMaxAlphaDeg + 1e-9;
         alpha += SeparationTableStepDeg) {
        fc.alphaDeg = alpha;
        const SolveResult solved = SolveWithSystem(prepared, fc, ref);
        const FlowField field = MakeFlowField(prepared, fc, solved.gamma, solved.sigma);
        const AttachmentLine line = ComputeAttachmentLine(field, wing, strips, sections);
        const SeparationSurvey survey = SurveySeparation(line);
        for (const StationSeparation& entry : survey.Stations) {
            if (!entry.Resolved || entry.Strip < 0) continue;
            const std::size_t i = static_cast<std::size_t>(entry.Strip);
            if (i >= strips.size()) continue;
            const Vec3 v = field.BoundMidpointVelocity(static_cast<int>(i));
            const double localDeg =
                Math::RadToDeg(std::atan2(Dot(v, strips[i].LiftDir), Dot(v, strips[i].ChordDir))) -
                strips[i].EffectiveAlpha0Deg();
            tables[i].AlphaDeg.push_back(localDeg);
            tables[i].Psi.push_back(entry.Upper.SeparationPsi);
        }
    }
    return tables;
}

/**
 * Nearest strip by span fraction, linear interpolation in local
 * incidence, LINEAR EXTRAPOLATION beyond the table's last point (clamped
 * to [0, 1]): the computed separation point keeps walking forward past
 * the table edge, and holding it constant instead would freeze f above
 * zero and deny the model its plate limit.
 */
[[nodiscard]] inline SeparationPointFunction MakeSeparationFunction(
    std::vector<StripSeparationTable> tables, std::vector<double> etas) {
    return [tables = std::move(tables), etas = std::move(etas)](double eta, double aDeg) {
        std::size_t best = 0;
        double bestDist = 1e30;
        for (std::size_t i = 0; i < etas.size(); ++i) {
            const double d = std::fabs(etas[i] - eta);
            if (d < bestDist) {
                bestDist = d;
                best = i;
            }
        }
        const StripSeparationTable& t = tables[best];
        const std::size_t n = t.AlphaDeg.size();
        if (n == 0) return 1.0;
        const double a = std::fabs(aDeg);
        if (a <= t.AlphaDeg.front()) return std::clamp(t.Psi.front(), 0.0, 1.0);
        std::size_t hi = 1;
        while (hi + 1 < n && t.AlphaDeg[hi] < a) ++hi;
        const double a0 = t.AlphaDeg[hi - 1], a1 = t.AlphaDeg[hi];
        const double frac = (a1 - a0 > 1e-9) ? (a - a0) / (a1 - a0) : 0.0; // >1 extrapolates
        return std::clamp(t.Psi[hi - 1] + frac * (t.Psi[hi] - t.Psi[hi - 1]), 0.0, 1.0);
    };
}

} // namespace Aeolion::Solver
