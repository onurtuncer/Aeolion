// Solver/DiscreteVortexSection.h
//
// Tier 3 of the post-separation study, two-dimensional half: an unsteady
// discrete-vortex section with leading-edge shedding modulated by the
// leading-edge suction parameter (LESP) -- the LDVM of Ramesh et al.
// (J. Fluid Mech. 751, 2014) in its fixed-incidence form. This is the
// UNSTEADY CROSS-CHECK of the quasi-steady tiers: it does not replace the
// anchored section model in the coupling, it measures what the coupling's
// limit-cycle means leave out -- the mean and fluctuating loads of the
// actual shedding flow at a given deep incidence.
//
// The formulation is large-angle unsteady thin-airfoil theory:
//
//   - The plate's bound vorticity is the Glauert series
//         gamma(theta) = 2 U [ A0 (1+cos)/sin + sum_n An sin(n theta) ],
//     with the time-varying coefficients projected each step from the
//     downwash of everything that is NOT the bound sheet (freestream +
//     every free vortex):  A0 = (1/pi) int (w/U) dtheta,
//     An = -(2/pi) int (w/U) cos(n theta) dtheta. No influence matrix is
//     ever assembled or solved.
//   - One trailing-edge vortex is shed every step; its strength satisfies
//     Kelvin's theorem. A0 -- the LESP -- measures the leading-edge
//     suction peak, and while |A0| exceeds the critical value a
//     leading-edge vortex is shed as well, with the pair's strengths
//     satisfying Kelvin AND |A0| = LESP_crit simultaneously. Both
//     conditions are LINEAR in the unknown strengths (the A-projections
//     of a unit vortex are just quadratures), so each step closes with a
//     1x1 or 2x2 linear solve -- no inner iteration, and the Kelvin
//     residual is machine zero by construction (asserted in the result).
//   - Free vortices are Vatistas-core regularized (core ~ 1.3 U dt, the
//     LDVM's own choice) and convect with the full local velocity,
//     including the bound sheet's, by forward Euler. Vortices far
//     downstream are dropped; their circulation is retired from the
//     Kelvin ledger explicitly so the bookkeeping stays exact.
//
// Loads come from the IMPULSE theorem, not from the linearized A-series
// force formulas: F = rho d/dt [ sum_k Gamma_k (z_k, -x_k) + (0, -Bx) ],
// with Bx the bound sheet's own first moment (closed-form in the A's)
// and the sum over every tracked free vortex, in the lifting sign
// convention (Gamma > 0 clockwise, L = rho U Gamma). The A-series normal
// force is a small-disturbance formula and was measured collapsing at
// deep incidence (cd(90) less than half the bluff-plate level) while the
// impulse form is exact for tracked vorticity at any angle -- the
// attached limit (Wagner, 2 pi sin alpha) falls out of the same
// expression through the starting vortex's recession. Vortices retired
// past the far cutoff keep contributing their impulse RATE analytically,
// as circulation drifting with the freestream. The quarter-chord moment
// keeps the steady series term only; its apparent-mass counterpart is
// omitted, which biases the INSTANTANEOUS cm during strong shedding but
// not the converged means this tier exists to report.
//
// Stated limitations: flat plate (no camber -- the cross-check question is
// posed at deep incidence where camber is secondary); fixed incidence
// (no pitching terms, U constant); Euler convection; no viscous decay of
// the free vortices. The LESP_crit default (0.2) is the literature's
// Re ~ 1e5-1e6 range, and it is a PARAMETER of the cross-check, not a fitted
// constant of the model under test. And one bias is inherent to the
// dimension, not the implementation: near 90 degrees a two-dimensional
// vortex street stays perfectly coherent -- no spanwise breakup -- and
// overpredicts the mean drag of a normal plate by roughly 60-70% (the
// classical 2-D-simulation result, ~3.3 against the measured ~2). The
// cross-check is therefore trusted for the stall break and mid post-stall
// range and for FLUCTUATION content; its own 90-degree overshoot is the
// standing argument for the three-dimensional particle tier.

#pragma once

#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

namespace Aeolion::Solver {

// --- discretization defaults -------------------------------------------------
inline constexpr double DvmDefaultTimeStep = 0.02;    ///< Convective, c/U.
inline constexpr double DvmDefaultDuration = 25.0;    ///< Convective times.
inline constexpr double DvmDefaultLespCritical = 0.2; ///< Ramesh's Re ~ 1e5..1e6 band.
inline constexpr int    DvmFourierTerms = 3;          ///< A0..A2 (loads need no more).
inline constexpr int    DvmChordStations = 64;        ///< Quadrature stations on the chord (the
                                                      ///< newest TEV sits close behind the TE and
                                                      ///< its downwash peak needs resolving).
inline constexpr double DvmCoreFactor = 1.3;          ///< Core radius = factor * U * dt.
inline constexpr double DvmWakeCutoffChords = 9.0;    ///< Drop vortices this far downstream.
inline constexpr double DvmShedFraction = 1.0 / 3.0;  ///< New vortex at this fraction toward the last.

struct DiscreteVortexOptions {
    double TimeStep = DvmDefaultTimeStep;
    double Duration = DvmDefaultDuration;
    double LespCritical = DvmDefaultLespCritical;
    double AverageFraction = 0.5; ///< Statistics over this trailing fraction of the run.
    bool KeepHistory = false;     ///< Retain the per-step load trace.
};

struct DiscreteVortexSample {
    double t = 0.0; ///< Convective time.
    double cl = 0.0, cd = 0.0, cm = 0.0;
    double A0 = 0.0;
    bool LevActive = false;
};

struct DiscreteVortexResult {
    bool Valid = false;
    double MeanCl = 0.0, MeanCd = 0.0, MeanCm = 0.0, MeanCn = 0.0;
    double RmsCl = 0.0, RmsCd = 0.0;
    int TevShed = 0, LevShed = 0;
    /** max |total circulation| / (U c) over the run -- machine zero by construction. */
    double KelvinResidual = 0.0;
    std::vector<DiscreteVortexSample> History;
};

namespace Detail {

struct FreeVortex {
    double x = 0.0, z = 0.0, gamma = 0.0;
};

/** Vatistas (n = 2) regularized point-vortex velocity. */
inline void VortexVelocity(double dx, double dz, double gamma, double core2, double& u,
                           double& w) {
    const double r2 = dx * dx + dz * dz;
    const double denom = 2.0 * std::numbers::pi * std::sqrt(r2 * r2 + core2 * core2);
    if (!(denom > Math::Tiny)) {
        u = 0.0;
        w = 0.0;
        return;
    }
    u = gamma * dz / denom;
    w = -gamma * dx / denom;
}

} // namespace Detail

/**
 * Fixed-incidence LDVM run on a unit-chord flat plate, impulsively started
 * in a unit freestream at `alphaDeg`. Returns converged mean/RMS loads
 * over the trailing AverageFraction of the run.
 */
[[nodiscard]] inline DiscreteVortexResult SolveDiscreteVortexSection(
    double alphaDeg, const DiscreteVortexOptions& options = {}) {
    DiscreteVortexResult res;
    const double dt = options.TimeStep;
    const int steps = static_cast<int>(options.Duration / dt);
    if (!(dt > 0.0) || steps < 8) return res;

    const double alpha = Math::DegToRad(alphaDeg);
    const double ca = std::cos(alpha), sa = std::sin(alpha);
    // Plate along x in [0, 1]; freestream (ca, sa); plate normal +z.
    const double core = DvmCoreFactor * dt;
    const double core2 = core * core;

    // Chord quadrature stations (cosine): theta midpoint rule.
    const int m = DvmChordStations;
    std::vector<double> theta(m), xs(m), dtheta(m);
    for (int j = 0; j < m; ++j) {
        theta[j] = std::numbers::pi * (j + 0.5) / m;
        xs[j] = 0.5 * (1.0 - std::cos(theta[j]));
        dtheta[j] = std::numbers::pi / m;
    }

    // A-projections of a downwash distribution sampled at the stations.
    const auto project = [&](const std::vector<double>& w, double coeff[DvmFourierTerms]) {
        for (int n = 0; n < DvmFourierTerms; ++n) coeff[n] = 0.0;
        for (int j = 0; j < m; ++j) {
            coeff[0] += w[j] * dtheta[j] / std::numbers::pi;
            for (int n = 1; n < DvmFourierTerms; ++n)
                coeff[n] -= 2.0 * w[j] * std::cos(n * theta[j]) * dtheta[j] / std::numbers::pi;
        }
    };

    std::vector<Detail::FreeVortex> wake;
    wake.reserve(2 * steps);
    int lastTev = -1, lastLev = -1; // indices into wake, for placement
    bool levActivePrev = false;     // placement continuity holds only across CONSECUTIVE sheds
    double retired = 0.0;           // circulation dropped past the cutoff
    double A[DvmFourierTerms] = {0, 0, 0};

    // Bound-sheet discrete equivalent for convecting the wake: vortices of
    // strength gamma(theta_j) (dx/dtheta) dtheta at the stations.
    std::vector<double> boundGamma(m, 0.0);
    const auto rebuildBound = [&]() {
        for (int j = 0; j < m; ++j) {
            const double st = std::sin(theta[j]);
            const double g = 2.0 * (A[0] * (1.0 + std::cos(theta[j])) / std::max(st, 1e-9) +
                                    A[1] * st + A[2] * std::sin(2.0 * theta[j]));
            boundGamma[j] = g * 0.5 * st * dtheta[j]; // gamma * dx, dx = sin/2 dtheta
        }
    };

    // Downwash (z-component, positive up) at the stations from freestream
    // plus the CURRENT wake; unit-vortex columns are added separately.
    std::vector<double> wBase(m), wUnitTev(m), wUnitLev(m);
    const auto sampleWake = [&](std::vector<double>& w) {
        for (int j = 0; j < m; ++j) {
            double sum = sa; // freestream normal component
            for (const Detail::FreeVortex& v : wake) {
                double u, wz;
                Detail::VortexVelocity(xs[j] - v.x, -v.z, v.gamma, core2, u, wz);
                sum += wz;
            }
            w[j] = sum;
        }
    };
    const auto sampleUnit = [&](double vx, double vz, std::vector<double>& w) {
        for (int j = 0; j < m; ++j) {
            double u, wz;
            Detail::VortexVelocity(xs[j] - vx, -vz, 1.0, core2, u, wz);
            w[j] = wz;
        }
    };

    const double avgStart = options.Duration * (1.0 - options.AverageFraction);
    double sumCl = 0.0, sumCd = 0.0, sumCm = 0.0, sumCn = 0.0, sumCl2 = 0.0, sumCd2 = 0.0;
    int avgCount = 0;
    double pxPrev = 0.0, pzPrev = 0.0; // vorticity impulse at the previous step

    for (int step = 0; step < steps; ++step) {
        const double t = (step + 1) * dt;

        // --- placements -----------------------------------------------------
        double tevX, tevZ;
        if (lastTev >= 0) {
            tevX = 1.0 + DvmShedFraction * (wake[lastTev].x - 1.0);
            tevZ = DvmShedFraction * wake[lastTev].z;
        } else {
            tevX = 1.0 + 0.5 * dt * ca;
            tevZ = 0.5 * dt * sa;
        }

        // --- solve the new strengths ---------------------------------------
        sampleWake(wBase);
        sampleUnit(tevX, tevZ, wUnitTev);
        double Abase[DvmFourierTerms], Atev[DvmFourierTerms];
        project(wBase, Abase);
        project(wUnitTev, Atev);

        // Kelvin: Gamma_b + sum(wake) + retired + Gamma_tev (+ Gamma_lev) = 0,
        // with Gamma_b = pi c U (A0 + A1/2) linear in every strength.
        double wakeSum = retired;
        for (const Detail::FreeVortex& v : wake) wakeSum += v.gamma;
        const auto gb = [&](const double a[DvmFourierTerms]) {
            return std::numbers::pi * (a[0] + 0.5 * a[1]);
        };

        double gTev = 0.0, gLev = 0.0;
        bool levActive = false;
        {
            // TEV alone first.
            const double k = gb(Atev) + 1.0;
            gTev = -(gb(Abase) + wakeSum) / k;
            const double a0 = Abase[0] + gTev * Atev[0];
            if (std::fabs(a0) > options.LespCritical) {
                levActive = true;
                double levX, levZ;
                if (levActivePrev && lastLev >= 0) {
                    levX = DvmShedFraction * wake[lastLev].x;
                    levZ = DvmShedFraction * wake[lastLev].z;
                } else {
                    // Off the leading edge, along the local upstream normal side.
                    levX = -0.5 * dt * ca;
                    levZ = 0.5 * dt * ((a0 > 0.0) ? sa + 0.5 : sa - 0.5);
                }
                sampleUnit(levX, levZ, wUnitLev);
                double Alev[DvmFourierTerms];
                project(wUnitLev, Alev);

                // 2x2 linear system: Kelvin and A0 = sign * crit.
                const double target = (a0 > 0.0) ? options.LespCritical : -options.LespCritical;
                const double a11 = gb(Atev) + 1.0, a12 = gb(Alev) + 1.0;
                const double a21 = Atev[0], a22 = Alev[0];
                const double b1 = -(gb(Abase) + wakeSum);
                const double b2 = target - Abase[0];
                const double det = a11 * a22 - a12 * a21;
                if (std::fabs(det) > Math::Tiny) {
                    gTev = (b1 * a22 - b2 * a12) / det;
                    gLev = (a11 * b2 - a21 * b1) / det;
                    wake.push_back({levX, levZ, gLev});
                    lastLev = static_cast<int>(wake.size()) - 1;
                    ++res.LevShed;
                } else {
                    levActive = false;
                }
            }
        }
        wake.push_back({tevX, tevZ, gTev});
        lastTev = static_cast<int>(wake.size()) - 1;
        ++res.TevShed;

        // --- final coefficients and the Kelvin ledger ------------------------
        sampleWake(wBase); // wake now includes the new vortices
        project(wBase, A);
        double ledger = gb(A) + retired;
        for (const Detail::FreeVortex& v : wake) ledger += v.gamma;
        res.KelvinResidual = std::max(res.KelvinResidual, std::fabs(ledger));

        // --- loads (impulse) --------------------------------------------------
        // Bound sheet's first moment, closed-form from the Glauert series:
        // integral gamma x dx = (pi U c^2 / 4)(A0 + A1 - A2/2).
        const double bx = 0.25 * std::numbers::pi * (A[0] + A[1] - 0.5 * A[2]);
        double px = 0.0, pz = -bx;
        for (const Detail::FreeVortex& v : wake) {
            px += v.gamma * v.z;
            pz -= v.gamma * v.x;
        }
        double fx = 0.0, fz = 0.0;
        if (step > 0) {
            fx = (px - pxPrev) / dt + retired * sa;  // retired circulation drifts
            fz = (pz - pzPrev) / dt - retired * ca;  // with the freestream
        }
        pxPrev = px;
        pzPrev = pz;
        const double cn = 2.0 * fz;                  // plate-normal force coefficient
        const double cl = 2.0 * (fz * ca - fx * sa); // wind axes
        const double cd = 2.0 * (fx * ca + fz * sa);
        const double cm = -0.25 * std::numbers::pi * (A[1] - A[2]);

        if (options.KeepHistory)
            res.History.push_back({t, cl, cd, cm, A[0], levActive});
        if (t >= avgStart) {
            sumCl += cl;
            sumCd += cd;
            sumCm += cm;
            sumCn += cn;
            sumCl2 += cl * cl;
            sumCd2 += cd * cd;
            ++avgCount;
        }
        levActivePrev = levActive;

        // --- convect --------------------------------------------------------
        rebuildBound();
        std::vector<Detail::FreeVortex> next = wake;
        for (std::size_t i = 0; i < wake.size(); ++i) {
            double u = ca, w = sa;
            for (std::size_t k = 0; k < wake.size(); ++k) {
                if (k == i) continue;
                double du, dw;
                Detail::VortexVelocity(wake[i].x - wake[k].x, wake[i].z - wake[k].z,
                                       wake[k].gamma, core2, du, dw);
                u += du;
                w += dw;
            }
            for (int j = 0; j < m; ++j) {
                double du, dw;
                Detail::VortexVelocity(wake[i].x - xs[j], wake[i].z, boundGamma[j], core2, du,
                                       dw);
                u += du;
                w += dw;
            }
            next[i].x += u * dt;
            next[i].z += w * dt;
        }
        wake = std::move(next);

        // --- retire the far wake -------------------------------------------
        // Indices shift; the shed-placement anchors are re-found by identity
        // of position, cheaper bookkeeping than stable handles at this size.
        const double tevXAnchor = wake[lastTev].x, tevZAnchor = wake[lastTev].z;
        const double levXAnchor = (lastLev >= 0) ? wake[lastLev].x : 0.0;
        const double levZAnchor = (lastLev >= 0) ? wake[lastLev].z : 0.0;
        std::vector<Detail::FreeVortex> kept;
        kept.reserve(wake.size());
        for (const Detail::FreeVortex& v : wake) {
            if (v.x > DvmWakeCutoffChords) {
                retired += v.gamma;
                // The vortex leaves the tracked impulse sum: remove it from
                // the PREVIOUS baseline too, or the next difference reads
                // its disappearance as a force spike of Gamma x/dt -- the
                // bug that showed up as steady negative drag at moderate
                // incidence, where retirement runs at one vortex per step.
                // From here on its impulse RATE is carried analytically by
                // the freestream-drift term on `retired`.
                pxPrev -= v.gamma * v.z;
                pzPrev += v.gamma * v.x;
                continue;
            }
            kept.push_back(v);
        }
        wake = std::move(kept);
        lastTev = lastLev = -1;
        for (std::size_t i = 0; i < wake.size(); ++i) {
            if (wake[i].x == tevXAnchor && wake[i].z == tevZAnchor) lastTev = static_cast<int>(i);
            if (wake[i].x == levXAnchor && wake[i].z == levZAnchor) lastLev = static_cast<int>(i);
        }
    }

    if (avgCount == 0) return res;
    res.Valid = true;
    res.MeanCl = sumCl / avgCount;
    res.MeanCd = sumCd / avgCount;
    res.MeanCm = sumCm / avgCount;
    res.MeanCn = sumCn / avgCount;
    res.RmsCl = std::sqrt(std::max(sumCl2 / avgCount - res.MeanCl * res.MeanCl, 0.0));
    res.RmsCd = std::sqrt(std::max(sumCd2 / avgCount - res.MeanCd * res.MeanCd, 0.0));
    return res;
}

} // namespace Aeolion::Solver
