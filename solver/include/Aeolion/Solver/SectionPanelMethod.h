// Solver/SectionPanelMethod.h
//
// A 2-D Hess-Smith panel solve on a section's REAL contour: constant-strength
// source panels carrying the thickness, one constant vortex strength over
// the whole surface carrying the circulation, and a Kutta condition at the
// trailing edge closing the system.
//
// --- what this is for --------------------------------------------------------
// The lattice cannot locate a wing's stagnation point. It is a camber sheet
// of zero thickness; the flow never stops on it. What the lattice DOES know
// is what each strip sees -- the local incidence, induced flow and all --
// and that is exactly the input a section problem needs. So the division of
// labour is: the 3-D lattice owns the induced field, and this owns the
// section, in the same way ViscousCoupling.h already lets a section model
// own the lift curve.
//
// Two things come out, and the second is the one that matters downstream:
//
//   the stagnation point,  where the surface velocity changes sign, and
//   U_e(s) on both surfaces, measured from that point.
//
// A stagnation point on its own would be half an answer. An integral
// boundary layer marches from the attachment point along the edge velocity
// distribution, and needs the strain rate dU_e/ds at the stagnation point
// to start at all: the Hiemenz similarity solution gives the initial
// momentum thickness there as theta_0^2 = 0.075 nu / (dU_e/ds), which is
// finite and NOT the zero a flat-plate start assumes. Both surfaces' runs
// therefore come out of the same solve that located the point.
//
// --- formulation --------------------------------------------------------------
// Classical Hess-Smith (Moran, "An Introduction to Theoretical and
// Computational Aerodynamics", 4.10; Katz & Plotkin 11.4). N panels carry
// N unknown source strengths plus one shared vortex strength, against N
// flow-tangency equations plus one Kutta equation, so the system is
// (N+1)x(N+1) and is solved by the same LAPACK-backed dense factorization
// the 3-D lattice uses.
//
// Sources alone cannot lift and a distributed vortex alone cannot make a
// closed body's thickness; the pair does both, and the single shared vortex
// strength -- rather than one per panel -- is what keeps the system square
// with exactly one Kutta condition.
//
// The contour must be the CLOCKWISE loop Geometry/SectionContour.h produces,
// because the outward normal is taken as (-sin theta, cos theta) with no
// per-panel sign test, and because the Kutta condition is written on the
// first and last panels on the assumption that they are the two adjacent to
// the trailing edge.
//
// Incompressible and inviscid, like everything else in the toolkit: this
// gives the EDGE velocity a boundary layer is driven by, not the boundary
// layer itself.

#pragma once

#include "Aeolion/Geometry/SectionContour.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Solver/Solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

namespace Aeolion::Solver {

// --- guards -------------------------------------------------------------------
inline constexpr double SectionPanelLengthFloor = 1e-12; ///< A collapsed panel contributes nothing.
inline constexpr double SectionRadiusFloor = 1e-14;      ///< Keeps the influence logarithm finite.
inline constexpr int    MinSectionPanels = 4;

/**
 * The boundary layer's view of one surface: arc length measured FROM the
 * stagnation point, and the edge speed along it.
 *
 * "Upper" and "lower" name the trailing edge each run ends at, not every
 * point it passes through -- the upper run starts at the stagnation point,
 * which at positive incidence lies on the lower surface, wraps around the
 * nose, and only then runs aft. That wrap is not an artefact to be trimmed
 * off: it is real boundary-layer run length, and dropping it is the usual
 * way a coupled method ends up with a too-thin suction-side layer.
 */
struct SurfaceRun {
    std::vector<double> S;  ///< Arc from the stagnation point, increasing [chord fractions].
    std::vector<double> Ue; ///< Edge speed, positive [multiples of the oncoming speed].
    std::vector<double> Psi, Zeta; ///< Where each station sits on the section.

    [[nodiscard]] std::size_t Count() const { return S.size(); }
    [[nodiscard]] double Length() const { return S.empty() ? 0.0 : S.back(); }
};

/** Everything a Hess-Smith solve of one section produces. */
struct SectionSolution {
    bool Valid = false;

    // --- per panel (all aligned, one entry per panel, in contour order) ----
    std::vector<double> S;    ///< Arc length of the panel midpoint from the contour's start [chord].
    std::vector<double> Ue;   ///< SIGNED surface velocity along the traversal direction [Vinf].
    std::vector<double> Cp;   ///< 1 - Ue^2.
    std::vector<double> Psi, Zeta; ///< Panel midpoints, chord fractions.

    double cl = 0.0;    ///< Section lift coefficient from the total bound circulation.
    double Gamma = 0.0; ///< Vortex strength per unit length (the shared unknown).
    double Perimeter = 0.0;

    // --- the stagnation point ------------------------------------------------
    bool StagnationFound = false;
    double StagnationArc = 0.0;  ///< Its arc position along the contour [chord].
    double StagnationPsi = 0.0;  ///< x/c of the stagnation point.
    double StagnationZeta = 0.0; ///< z/c.

    /**
     * Whether the stagnation point sits on the lower surface. At incidence
     * above the zero-lift angle it does; below, it moves onto the upper
     * surface. The sign of this is the sign of the section's loading, which
     * is why it is stated rather than assumed.
     */
    bool StagnationOnLower = true;

    /**
     * Where the point actually lands, for calibration against intuition.
     * Measured from the leading edge along the surface, the offset follows
     *
     *     s_stag / c  ~  sqrt(2 r_LE / c) * alpha_e,
     *
     * with alpha_e the incidence from the zero-lift line. The SQUARE ROOT
     * of the nose radius, not the nose radius itself: an airfoil nose does
     * not sit in a uniform stream at angle alpha the way an isolated
     * cylinder does, it sits inside the leading-edge singularity of the
     * outer thin-airfoil solution, where u ~ alpha sqrt(c/s). Matching that
     * against the parabolic nose's own sqrt(s/r_LE) is what produces the
     * square root. The cylinder reading gives r_LE * alpha, which is an
     * order of magnitude too small on a typical section -- 0.001 c against
     * a true 0.013 c for a NACA 0012 at four degrees.
     *
     * That gap is the whole argument for solving the section rather than
     * correlating it, and it is asserted directly in
     * TestSectionPanelMethod.cpp across a twelvefold range of nose radius.
     */

    /**
     * dU_e/ds at the stagnation point, per chord and per unit oncoming
     * speed. Positive by construction. This is the Hiemenz strain rate that
     * sets a boundary layer's initial momentum thickness,
     *
     *     theta_0^2 = 0.075 nu / (dU_e/ds),
     *
     * so dimensionalize as (StagnationStrain * Vinf / chord) [1/s] before
     * using it.
     */
    double StagnationStrain = 0.0;

    [[nodiscard]] SurfaceRun UpperRun() const;
    [[nodiscard]] SurfaceRun LowerRun() const;
};

namespace Detail {

/** Per-panel geometry of a contour: endpoints, angle, length, midpoint. */
struct SectionPanels {
    int Count = 0;
    std::vector<double> X0, Z0, X1, Z1; ///< Endpoints.
    std::vector<double> Xm, Zm;         ///< Midpoints (control points).
    std::vector<double> Theta, Length, Arc; ///< Angle, length, arc of the midpoint from the start.
    double Perimeter = 0.0;
};

[[nodiscard]] inline SectionPanels MakeSectionPanels(const Geometry::SectionContour& contour) {
    SectionPanels panels;
    const std::size_t points = contour.Count();
    if (points < 2) return panels;

    double arc = 0.0;
    for (std::size_t k = 0; k + 1 < points; ++k) {
        const double dx = contour.Psi[k + 1] - contour.Psi[k];
        const double dz = contour.Zeta[k + 1] - contour.Zeta[k];
        const double length = std::hypot(dx, dz);
        if (length < SectionPanelLengthFloor) continue; // a repeated point carries no panel

        panels.X0.push_back(contour.Psi[k]);
        panels.Z0.push_back(contour.Zeta[k]);
        panels.X1.push_back(contour.Psi[k + 1]);
        panels.Z1.push_back(contour.Zeta[k + 1]);
        panels.Xm.push_back(Math::Half * (contour.Psi[k] + contour.Psi[k + 1]));
        panels.Zm.push_back(Math::Half * (contour.Zeta[k] + contour.Zeta[k + 1]));
        panels.Theta.push_back(std::atan2(dz, dx));
        panels.Length.push_back(length);
        panels.Arc.push_back(arc + Math::Half * length);
        arc += length;
    }
    panels.Count = static_cast<int>(panels.Length.size());
    panels.Perimeter = arc;
    return panels;
}

} // namespace Detail

/**
 * Solve one section contour at incidence `alphaRad`, measured from the
 * section's own chord line (the psi axis), positive nose-up.
 *
 * Everything is nondimensional: the contour is in chord fractions and the
 * oncoming speed is unity, so Ue comes back as a multiple of the oncoming
 * speed and arc lengths as fractions of chord.
 */
[[nodiscard]] inline SectionSolution SolveSectionContour(const Geometry::SectionContour& contour,
                                                          double alphaRad) {
    SectionSolution solution;
    if (!contour.Valid()) return solution;

    const Detail::SectionPanels panels = Detail::MakeSectionPanels(contour);
    const int n = panels.Count;
    if (n < MinSectionPanels) return solution;

    const std::size_t nn = static_cast<std::size_t>(n);
    const double cosAlpha = std::cos(alphaRad), sinAlpha = std::sin(alphaRad);
    const double inv2Pi = 1.0 / (2.0 * std::numbers::pi);

    // --- influence coefficients -----------------------------------------------
    // A: normal influence, B: tangential influence, both of a unit source on
    // panel j at the midpoint of panel i. The single shared VORTEX's
    // influences follow from the same numbers by the standard duality --
    // rotating a source's velocity by ninety degrees turns it into a
    // vortex's -- which is why they are summed from A and B rather than
    // computed again:
    //
    //     vortex normal at i     = -sum_j B_ij
    //     vortex tangential at i = +sum_j A_ij
    std::vector<double> A(nn * nn), B(nn * nn);
    std::vector<double> vortexNormal(nn, 0.0), vortexTangential(nn, 0.0);

    for (int i = 0; i < n; ++i) {
        const std::size_t ii = static_cast<std::size_t>(i);
        for (int j = 0; j < n; ++j) {
            const std::size_t jj = static_cast<std::size_t>(j);
            const double sinDiff = std::sin(panels.Theta[ii] - panels.Theta[jj]);
            const double cosDiff = std::cos(panels.Theta[ii] - panels.Theta[jj]);

            double logRatio = 0.0;
            double subtended = std::numbers::pi; // a panel subtends pi at its own midpoint
            if (i != j) {
                const double dx0 = panels.Xm[ii] - panels.X0[jj], dz0 = panels.Zm[ii] - panels.Z0[jj];
                const double dx1 = panels.Xm[ii] - panels.X1[jj], dz1 = panels.Zm[ii] - panels.Z1[jj];
                const double r0 = std::max(std::hypot(dx0, dz0), SectionRadiusFloor);
                const double r1 = std::max(std::hypot(dx1, dz1), SectionRadiusFloor);
                logRatio = std::log(r1 / r0);
                subtended = std::atan2(dz1 * dx0 - dx1 * dz0, dx1 * dx0 + dz1 * dz0);
            }

            A[ii * nn + jj] = inv2Pi * (sinDiff * logRatio + cosDiff * subtended);
            B[ii * nn + jj] = inv2Pi * (sinDiff * subtended - cosDiff * logRatio);
        }
        for (int j = 0; j < n; ++j) {
            vortexNormal[ii] -= B[ii * nn + static_cast<std::size_t>(j)];
            vortexTangential[ii] += A[ii * nn + static_cast<std::size_t>(j)];
        }
    }

    // --- assemble the (N+1) x (N+1) system --------------------------------------
    const std::size_t m = nn + 1;
    std::vector<double> matrix(m * m, 0.0);
    std::vector<double> rhs(m, 0.0);

    // N flow-tangency rows: -Vinf . n_i = sin(theta_i - alpha).
    for (int i = 0; i < n; ++i) {
        const std::size_t ii = static_cast<std::size_t>(i);
        for (int j = 0; j < n; ++j) matrix[ii * m + static_cast<std::size_t>(j)] = A[ii * nn + static_cast<std::size_t>(j)];
        matrix[ii * m + nn] = vortexNormal[ii];
        rhs[ii] = std::sin(panels.Theta[ii] - alphaRad);
    }

    // The Kutta condition: equal and opposite tangential speed on the two
    // panels flanking the trailing edge, V_t(first) + V_t(last) = 0. Their
    // traversal directions are opposite there, so equal MAGNITUDE is what
    // the sum being zero states.
    const std::size_t first = 0, last = nn - 1;
    for (int j = 0; j < n; ++j) {
        const std::size_t jj = static_cast<std::size_t>(j);
        matrix[nn * m + jj] = B[first * nn + jj] + B[last * nn + jj];
    }
    matrix[nn * m + nn] = vortexTangential[first] + vortexTangential[last];
    rhs[nn] = -(std::cos(panels.Theta[first] - alphaRad) + std::cos(panels.Theta[last] - alphaRad));

    const ConstDenseMatrixView view(matrix.data(), static_cast<int>(m), static_cast<int>(m));
    const LUFactorization factorization = LuFactorize(view);
    const std::vector<double> strengths = LuSolve(factorization, rhs);
    if (strengths.size() != m) return solution;

    const double gamma = strengths[nn];

    // --- surface velocity -----------------------------------------------------
    solution.S.resize(nn);
    solution.Ue.resize(nn);
    solution.Cp.resize(nn);
    solution.Psi.resize(nn);
    solution.Zeta.resize(nn);

    for (int i = 0; i < n; ++i) {
        const std::size_t ii = static_cast<std::size_t>(i);
        double tangential = std::cos(panels.Theta[ii] - alphaRad) + vortexTangential[ii] * gamma;
        for (int j = 0; j < n; ++j)
            tangential += B[ii * nn + static_cast<std::size_t>(j)] * strengths[static_cast<std::size_t>(j)];
        solution.S[ii] = panels.Arc[ii];
        solution.Ue[ii] = tangential;
        solution.Cp[ii] = 1.0 - tangential * tangential;
        solution.Psi[ii] = panels.Xm[ii];
        solution.Zeta[ii] = panels.Zm[ii];
    }

    // Total bound circulation is the shared strength times the perimeter it
    // is distributed over; cl = 2 Gamma on a unit chord at unit speed. The
    // contour runs clockwise while positive lift is a clockwise circulation
    // in these axes, which is where the sign comes from.
    solution.Gamma = gamma;
    solution.Perimeter = panels.Perimeter;
    solution.cl = 2.0 * gamma * panels.Perimeter;

    // --- the stagnation point --------------------------------------------------
    // Traversal runs trailing edge -> lower -> leading edge -> upper ->
    // trailing edge. On the lower surface the flow runs aft, AGAINST that
    // direction, so Ue is negative there; past the stagnation point the flow
    // runs around the nose WITH the traversal and Ue is positive. So the
    // stagnation point is the first negative-to-positive crossing.
    for (int i = 0; i + 1 < n; ++i) {
        const std::size_t ii = static_cast<std::size_t>(i);
        const double here = solution.Ue[ii], next = solution.Ue[ii + 1];
        if (!(here <= 0.0 && next > 0.0)) continue;

        const double span = solution.S[ii + 1] - solution.S[ii];
        const double fraction = (next - here != 0.0) ? -here / (next - here) : 0.0;
        solution.StagnationFound = true;
        solution.StagnationArc = solution.S[ii] + fraction * span;
        solution.StagnationPsi = solution.Psi[ii] + fraction * (solution.Psi[ii + 1] - solution.Psi[ii]);
        solution.StagnationZeta = solution.Zeta[ii] + fraction * (solution.Zeta[ii + 1] - solution.Zeta[ii]);
        // The lower stretch is the first half of the traversal, up to the
        // leading edge -- which is where psi stops decreasing.
        solution.StagnationOnLower = (solution.Psi[ii + 1] <= solution.Psi[ii]);
        solution.StagnationStrain = (span > SectionPanelLengthFloor) ? (next - here) / span : 0.0;
        break;
    }

    solution.Valid = true;
    return solution;
}

// --- the two boundary-layer runs ------------------------------------------------
inline SurfaceRun SectionSolution::UpperRun() const {
    SurfaceRun run;
    if (!Valid || !StagnationFound) return run;
    for (std::size_t k = 0; k < S.size(); ++k) {
        if (S[k] <= StagnationArc) continue;
        run.S.push_back(S[k] - StagnationArc);
        run.Ue.push_back(std::fabs(Ue[k]));
        run.Psi.push_back(Psi[k]);
        run.Zeta.push_back(Zeta[k]);
    }
    return run;
}

inline SurfaceRun SectionSolution::LowerRun() const {
    SurfaceRun run;
    if (!Valid || !StagnationFound) return run;
    // Walk BACK toward the contour's start, so arc from the stagnation point
    // still comes out increasing.
    for (std::size_t k = S.size(); k-- > 0;) {
        if (S[k] >= StagnationArc) continue;
        run.S.push_back(StagnationArc - S[k]);
        run.Ue.push_back(std::fabs(Ue[k]));
        run.Psi.push_back(Psi[k]);
        run.Zeta.push_back(Zeta[k]);
    }
    return run;
}

/**
 * Initial momentum thickness of a laminar boundary layer starting at a
 * stagnation point, from the Hiemenz similarity solution:
 *
 *     theta_0 = sqrt(0.075 nu / (dU_e/ds)).
 *
 * `strain` is SectionSolution::StagnationStrain, `chord` and `speed` are
 * what the section was nondimensionalized by, and `kinematicViscosity` is
 * in m^2/s. Returns metres. This is the number that a march starting from
 * theta = 0 at the leading edge gets wrong, and it matters most exactly
 * where the layer is thinnest and the transition correlations are most
 * sensitive.
 */
[[nodiscard]] inline double StagnationMomentumThickness(double strain, double chord, double speed,
                                                        double kinematicViscosity) {
    const double rate = strain * speed / std::max(chord, Math::Tiny); // [1/s]
    if (!(rate > Math::Tiny)) return 0.0;
    return std::sqrt(0.075 * kinematicViscosity / rate);
}

} // namespace Aeolion::Solver
