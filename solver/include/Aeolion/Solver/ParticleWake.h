// Solver/ParticleWake.h
//
// Tier 3 of the post-separation study, three-dimensional half: an unsteady
// single-row ring lattice whose wake is carried by vortex PARTICLES, shed
// from the trailing edge everywhere and from the leading edge of separated
// strips. Like its two-dimensional sibling (DiscreteVortexSection.h) this
// is a CROSS-CHECK, not a coupling: run at a fixed attitude it reports the
// mean and fluctuating loads of the actual three-dimensional shedding
// flow, the thing the quasi-steady map's limit-cycle means average away
// and the 2-D cross-check overpredicts for want of spanwise breakup.
//
// The pieces:
//
//   - Each strip carries a RING: leading segment on the bound quarter-chord
//     line, chordwise side legs, closing segment at the trailing edge;
//     flow tangency at the three-quarter-chord control point. The ring
//     influence matrix is geometry-only and factored once.
//   - The wake leaves through a one-step BUFFER ring per strip (the
//     standard unsteady-lattice device): a filament quad from the trailing
//     edge to one convection step downstream, carrying the strip's
//     PREVIOUS circulation. Its near segment cancels the bound ring's
//     closing segment -- without it a full Gamma of spurious spanwise
//     vorticity sits on the trailing edge, and the bound solve was
//     measured converging to 42% of the steady VLM's circulation before
//     this was understood. Each step the old buffer converts to particles
//     as four segment-particles (closed loops in, closed loops out, so
//     conservation needs no bookkeeping at all): the +/- spanwise pairs of
//     successive steps annihilate wherever the circulation is steady and
//     leave exactly the shed sheet where it is not, and the side pieces
//     superpose into the trailing sheet.
//   - A strip whose instantaneous incidence is past its separation
//     boundary (the same f(eta, alpha) tables the anchored model uses, or
//     a plain incidence threshold when none are supplied) sheds from the
//     LEADING edge as well: a particle pair carrying the boundary-layer
//     vorticity flux (1/2 V_loc^2 per unit span, the classical shear-layer
//     rate), spanwise-oriented, one sign at the leading edge and its exact
//     negative folded into that strip's trailing-edge shed -- total
//     circulation untouched, counter-rotating layers as a bluff section
//     sheds them. The pair's sign convention (leading edge carries the
//     bound sense) is pinned by the tests, not argued from first
//     principles here.
//   - Particles convect with the full local velocity (freestream, rings,
//     particles; algebraic-core regularized) and STRETCH: d(alpha)/dt =
//     (alpha . grad) u with the analytic kernel gradient, magnitude-capped
//     as a stability guard. Stretching is not optional equipment here --
//     spanwise breakup of the street is the very physics this tier exists
//     to add over the 2-D cross-check.
//   - Loads come from the impulse theorem, as in 2-D:
//         F = -rho d/dt [ sum Gamma (ring + buffer vector areas)
//                         + 1/2 sum x_p x alpha_p ],
//     exact for tracked vorticity at ANY incidence. Two alternatives were
//     tried and rejected with their reasons recorded: the naive
//     four-segment buffer conversion made the impulse a difference of
//     huge cancelling +/- pair terms (force RMS thirty times the mean
//     while the bound circulation was verified correct through the same
//     run), and LOCAL unsteady Kutta-Joukowski forces on the bound line
//     are smooth but structurally cannot produce bluff-plate drag -- at
//     90 degrees a line force perpendicular to the local velocity
//     averages to zero, and the pressure content lives chordwise. The
//     resolution is the MERGED conversion above: impulse exactness with
//     a wake of small-strength particles. Retired particles leave the
//     impulse baseline as they leave the sum (the 2-D lesson) and drift
//     analytically with the freestream after.
//
// Stated limitations: single-row lattice (the map's own resolution), no
// body/duct sources (the comparison target is the map's WING-ONLY forces),
// direct N^2 summation (OpenMP-parallel, no treecode), and an eddy
// viscosity that is a stated model parameter rather than resolved physics.
// Convection is second-order midpoint; diffusion is core spreading with
// far-wake merging; long-run consistency is kept by Pedrizzetti
// relaxation; convergence is a claim the batch-mean confidence interval
// makes, not the run length.

#pragma once

#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Solver/PostStallSection.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numbers>
#include <vector>

namespace Aeolion::Solver {

// --- defaults ---------------------------------------------------------------
inline constexpr double PwDefaultTimeStep = 0.08;   ///< Convective, ref-chord/V.
inline constexpr double PwDefaultDuration = 18.0;   ///< Convective times.
inline constexpr double PwCoreFactor = 1.3;         ///< Core = factor * V * dt.
inline constexpr double PwCutoffSpans = 6.0;        ///< Retire particles this far downstream.
inline constexpr double PwLeSheddingF = 0.35;       ///< Shed the LE when f falls below this...
inline constexpr double PwLeSheddingIncidenceDeg = 20.0; ///< ...or past this incidence, without tables.
inline constexpr double PwStretchCap = 4.0;         ///< |alpha| growth limit, multiples of birth strength.
/**
 * Ceiling on the local speed entering the leading-edge flux, as a multiple
 * of the freestream. The flux is quadratic in the local speed, and the
 * local speed includes the particles the flux itself shed -- a positive
 * feedback that was measured running away at deep incidence on the real
 * (non-uniform, trimmed) wing while the uniform test fixture stayed
 * subcritical. A plate edge's potential-flow speedup is order two, so the
 * ceiling states known physics rather than tuning the answer.
 */
inline constexpr double PwLeFluxSpeedCap = 2.5;
inline constexpr std::size_t PwMaxParticles = 30000;///< Runaway guard; exceeding it invalidates the run.

// --- Phase A: the mean-capable upgrades -------------------------------------
// An INVISCID particle street at this resolution has no dissipation: the
// enstrophy piles up at the core scale and the fluctuation-to-mean ratio
// GROWS with the averaging window (measured 3-40 across the configuration
// attitudes). Core spreading with an eddy viscosity is the physical
// regularizer; its coefficient is THE dissipation parameter of this tier,
// dimensionally nu_t = coeff * Vinf * ref-chord, its sensitivity reported
// with the runs rather than hidden. Zero disables diffusion (and with it
// merging), recovering the inviscid tier exactly.
inline constexpr double PwTurbulentViscosityCoeff = 1e-3;
inline constexpr int    PwRelaxEvery = 5;      ///< Pedrizzetti realignment cadence [steps].
inline constexpr double PwRelaxFraction = 0.3; ///< Blend toward the local vorticity direction.
inline constexpr int    PwMergeEvery = 10;     ///< Far-wake merge cadence [steps].
inline constexpr double PwMergeDistanceFactor = 0.5; ///< Merge when closer than this x core.
inline constexpr double PwMergeZoneChords = 1.5;     ///< Merging only this far behind the TE.
inline constexpr double PwBatchConvectiveTimes = 5.0;///< Statistics batch length ~ a shedding period.

struct ParticleWakeOptions {
    double TimeStep = PwDefaultTimeStep;
    double Duration = PwDefaultDuration;
    double AverageFraction = 0.5;
    /** f(eta, |alpha - alpha0|) as the anchored model consumes it; empty
     *  falls back to the plain incidence threshold. */
    SeparationPointFunction SeparationPoint;
    bool KeepHistory = false;
    /** nu_t = this x Vinf x ref chord; zero = inviscid (no spreading, no merge). */
    double TurbulentViscosityCoeff = PwTurbulentViscosityCoeff;
};

struct ParticleWakeSample {
    double t = 0.0;
    double CL = 0.0, CD = 0.0, CY = 0.0, CN = 0.0;
    int Particles = 0;
};

struct ParticleWakeResult {
    bool Valid = false;
    double MeanCL = 0.0, MeanCD = 0.0, MeanCY = 0.0, MeanCN = 0.0;
    double RmsCL = 0.0, RmsCN = 0.0;
    /**
     * The Kutta-Joukowski reading of the same window, 2 sum(Gamma w)/(V S):
     * what the bound circulation alone claims the lift is. Attached, it
     * must agree with MeanCL -- disagreement separates a bound-solve
     * defect from an impulse-accounting one, which is exactly the
     * diagnostic that shaped this header.
     */
    double MeanCirculationCL = 0.0;
    /**
     * Batch-mean 95% half-width on MeanCN (batches of PwBatchConvectiveTimes;
     * ~2 sigma / sqrt(batches)). "Converged" is a claim THIS number makes,
     * not the run length.
     */
    double MeanCN_CI = 0.0;
    int Batches = 0;
    int Merged = 0; ///< Particles absorbed by far-wake merging.
    int MaxParticles = 0;
    int SheddingStrips = 0; ///< Strips that shed from the leading edge at the last step.
    std::vector<ParticleWakeSample> History;
};

namespace Detail {

struct WakeParticle {
    Vec3 X{0, 0, 0};
    Vec3 Alpha{0, 0, 0};   ///< Vector strength, circulation x length.
    double Birth = 0.0;    ///< |Alpha| at creation, for the stretch cap.
    double Core2 = 0.0;    ///< Squared core radius; grows by core spreading.
};

/** Regularized particle-induced velocity, v = (alpha x r) / (4 pi rho^3). */
inline Vec3 ParticleVelocity(const Vec3& at, const WakeParticle& p, double core2) {
    const Vec3 r = at - p.X;
    const double rho2 = Dot(r, r) + core2;
    const double inv = 1.0 / (4.0 * std::numbers::pi * rho2 * std::sqrt(rho2));
    return Cross(p.Alpha, r) * inv;
}

/** The kernel's velocity gradient at `at`, for the stretching term. */
inline void ParticleVelocityGradient(const Vec3& at, const WakeParticle& p, double core2,
                                     Vec3 grad[3]) {
    const Vec3 r = at - p.X;
    const double rho2 = Dot(r, r) + core2;
    const double inv3 = 1.0 / (4.0 * std::numbers::pi * rho2 * std::sqrt(rho2));
    const double inv5 = 3.0 * inv3 / rho2;
    // v_i = eps_ijk a_j r_k inv3;  dv_i/dx_l = eps_ijl a_j inv3 - v_i^hat r_l inv5-ish
    const Vec3 axr = Cross(p.Alpha, r);
    const Vec3 unit[3] = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    for (int l = 0; l < 3; ++l) {
        const Vec3 depsilon = Cross(p.Alpha, unit[l]) * inv3;
        const double rl = (l == 0) ? r.x : (l == 1) ? r.y : r.z;
        grad[l] = depsilon - axr * (inv5 * rl);
    }
}

/** Robust finite-segment Biot-Savart (unit circulation). */
inline Vec3 SegmentVelocityUnit(const Vec3& at, const Vec3& p1, const Vec3& p2, double core2) {
    const Vec3 r1 = at - p1, r2 = at - p2;
    const double l1 = r1.Norm(), l2 = r2.Norm();
    const Vec3 cr = Cross(r1, r2);
    const double denom = 4.0 * std::numbers::pi *
                         (l1 * l2 * (l1 * l2 + Dot(r1, r2)) + core2 * (l1 + l2) * (l1 + l2));
    if (!(denom > Math::Tiny)) return {0, 0, 0};
    return cr * ((l1 + l2) / denom);
}

} // namespace Detail

/**
 * Fixed-attitude particle-wake cross-check on a single-row lattice. `strips`
 * must align one-to-one with `panels` and carry TRUE-chord frames (the
 * same contract as the coupled solve). Loads are wing-only, wind-axis,
 * normalized by `ref`.
 */
[[nodiscard]] inline ParticleWakeResult SolveParticleWake(
    const std::vector<Panel>& panels, const std::vector<StripSection>& strips,
    const FreestreamConditions& fc, const ReferenceGeometry& ref,
    const ParticleWakeOptions& options = {}) {
    ParticleWakeResult res;
    const std::size_t n = panels.size();
    if (n == 0 || strips.size() != n || !(ref.Chord > 0.0) || !(fc.Vinf > 0.0)) return res;

    const double dtPhys = options.TimeStep * ref.Chord / fc.Vinf;
    const int steps = static_cast<int>(options.Duration * ref.Chord / (fc.Vinf * dtPhys));
    const double core = PwCoreFactor * fc.Vinf * dtPhys;
    const double core2 = core * core;
    const Vec3 Vinf = FreestreamVelocity(fc);
    const double cutoff = PwCutoffSpans * ref.Span;

    // --- ring geometry -------------------------------------------------------
    struct Ring {
        Vec3 C[4];     ///< A, B, TE_B, TE_A: bound segment first, loop closed.
        Vec3 SVec;     ///< 1/2 closed-loop integral of x x dl (vector area).
        Vec3 Mid, Te;  ///< Bound midpoint and trailing-edge midpoint.
        Vec3 Span;     ///< B - A.
    };
    std::vector<Ring> rings(n);
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 aft = strips[i].ChordDir * (0.75 * strips[i].Chord);
        Ring& r = rings[i];
        r.C[0] = panels[i].A;
        r.C[1] = panels[i].B;
        r.C[2] = panels[i].B + aft;
        r.C[3] = panels[i].A + aft;
        r.Mid = (panels[i].A + panels[i].B) * 0.5;
        r.Te = r.Mid + aft;
        r.Span = panels[i].B - panels[i].A;
        r.SVec = Vec3(0, 0, 0);
        for (int k = 0; k < 4; ++k) r.SVec = r.SVec + Cross(r.C[k], r.C[(k + 1) % 4]) * 0.5;
    }
    const auto ringVelocityUnit = [&](const Vec3& at, const Ring& r) {
        Vec3 v(0, 0, 0);
        for (int k = 0; k < 4; ++k)
            v = v + Detail::SegmentVelocityUnit(at, r.C[k], r.C[(k + 1) % 4], 1e-12);
        return v;
    };

    // Influence matrix (geometry-only), factored once.
    std::vector<double> a(n * n);
    DenseMatrixView view(a.data(), static_cast<int>(n), static_cast<int>(n));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            view[static_cast<int>(i), static_cast<int>(j)] =
                Dot(ringVelocityUnit(panels[i].ControlPoint, rings[j]), panels[i].Normal);
    const LUFactorization lu = LuFactorize(view);

    // --- the buffer wake ring --------------------------------------------------
    // Fixed geometry (recreated at the TE each step), carrying gamma(t-1).
    const Vec3 bufStep = Vinf * dtPhys;
    struct Buffer {
        Vec3 C[4];
        Vec3 SVec;
    };
    std::vector<Buffer> buffers(n);
    for (std::size_t i = 0; i < n; ++i) {
        Buffer& b = buffers[i];
        b.C[0] = rings[i].C[3];           // TE_A: near segment runs A -> B,
        b.C[1] = rings[i].C[2];           // TE_B: cancelling the ring's closer.
        b.C[2] = rings[i].C[2] + bufStep;
        b.C[3] = rings[i].C[3] + bufStep;
        b.SVec = Vec3(0, 0, 0);
        for (int k = 0; k < 4; ++k) b.SVec = b.SVec + Cross(b.C[k], b.C[(k + 1) % 4]) * 0.5;
    }
    const auto bufferVelocityUnit = [&](const Vec3& at, const Buffer& b) {
        Vec3 v(0, 0, 0);
        for (int k = 0; k < 4; ++k)
            v = v + Detail::SegmentVelocityUnit(at, b.C[k], b.C[(k + 1) % 4], 1e-12);
        return v;
    };

    // --- state ---------------------------------------------------------------
    std::vector<Detail::WakeParticle> particles;
    particles.reserve(4096);
    std::vector<double> gamma(n, 0.0), bufGamma(n, 0.0);
    std::vector<int> nearIdx(n, -1); ///< Each strip's transient near-particle (see conversion).
    // Spanwise adjacency from the geometry, not the index: strip i abuts
    // strip i-1 only if their bound endpoints actually meet.
    std::vector<bool> adjacentToLeft(n, false);
    for (std::size_t i = 1; i < n; ++i)
        adjacentToLeft[i] =
            (panels[i].A - panels[i - 1].B).Norm() < 0.05 * std::max(strips[i].Width, Math::Tiny);
    Vec3 pPrev(0, 0, 0);
    Vec3 retiredAlpha(0, 0, 0);
    bool havePrev = false;

    const auto particleVelocityAt = [&](const Vec3& at) {
        Vec3 v(0, 0, 0);
        for (const Detail::WakeParticle& p : particles)
            v = v + Detail::ParticleVelocity(at, p, p.Core2);
        for (std::size_t j = 0; j < n; ++j)
            if (std::fabs(bufGamma[j]) > Math::Tiny)
                v = v + bufferVelocityUnit(at, buffers[j]) * bufGamma[j];
        return v;
    };
    const double nuEff = options.TurbulentViscosityCoeff * fc.Vinf * ref.Chord;

    const double avgStart = options.Duration * (1.0 - options.AverageFraction);
    const int batchSteps =
        std::max(1, static_cast<int>(PwBatchConvectiveTimes / options.TimeStep));
    std::vector<double> batchMeans;
    double batchSum = 0.0;
    int batchFill = 0;
    const double qS = 0.5 * fc.rho * fc.Vinf * fc.Vinf *
                      ((ref.Area > 0.0) ? ref.Area : 1.0);
    const Vec3 dragDir = Vinf.Normalized();
    const Vec3 liftDir = Cross(dragDir, Vec3(0, 1, 0)).Normalized();
    const Vec3 sideDir = Cross(liftDir, dragDir).Normalized();
    double sCL = 0, sCD = 0, sCY = 0, sCN = 0, sCL2 = 0, sCN2 = 0, sCircCL = 0;
    int avgCount = 0;

    for (int step = 0; step < steps; ++step) {
        const double t = (step + 1) * options.TimeStep;

        // --- bound solve ------------------------------------------------------
        std::vector<double> rhs(n);
        for (std::size_t i = 0; i < n; ++i)
            rhs[i] = -Dot(Vinf + particleVelocityAt(panels[i].ControlPoint), panels[i].Normal);
        gamma = LuSolve(lu, rhs);

        // --- convert the old buffer (merged), refill it, shed separated edges --
        // The naive four-segment conversion filled the wake with +/-
        // full-strength spanwise pairs whose impulse-difference noise
        // drowned the loads. Merged across steps instead: each strip keeps
        // ONE transient near-particle (+Gamma at the TE); the next step's
        // far-emission (-Gamma at TE+d) is ABSORBED into it -- by then it
        // has convected to that neighbourhood -- leaving a small particle
        // carrying exactly the SHED vorticity. Side legs emit per EDGE as
        // circulation differences, small by the same logic. The wake then
        // carries one bounded full-strength transient per strip and
        // small-strength particles everywhere else, and the impulse loads
        // are usable again.
        res.SheddingStrips = 0;
        const auto emit = [&](const Vec3& at, const Vec3& alpha) {
            Detail::WakeParticle p;
            p.X = at;
            p.Alpha = alpha;
            p.Birth = alpha.Norm();
            p.Core2 = core2;
            if (p.Birth > Math::Tiny) {
                particles.push_back(p);
                return static_cast<int>(particles.size()) - 1;
            }
            return -1;
        };
        for (std::size_t i = 0; i < n; ++i) {
            const Buffer& b = buffers[i];
            const Vec3 spanW = rings[i].Span.Normalized() * strips[i].Width;
            // Absorb this strip's far-emission into its transient near
            // particle from the previous step; fresh transient after.
            if (std::fabs(bufGamma[i]) > Math::Tiny || nearIdx[i] >= 0) {
                if (nearIdx[i] >= 0) {
                    Detail::WakeParticle& nearP = particles[static_cast<std::size_t>(nearIdx[i])];
                    nearP.Alpha = nearP.Alpha - spanW * bufGamma[i];
                    nearP.Birth = std::max(nearP.Alpha.Norm(), nearP.Birth * 0.25);
                } else {
                    emit((b.C[2] + b.C[3]) * 0.5, spanW * (-bufGamma[i]));
                }
                nearIdx[i] = emit((b.C[0] + b.C[1]) * 0.5, spanW * bufGamma[i]);
            }
            // Trailing vorticity, merged per edge: this strip's LEFT edge
            // against its left neighbour -- where a neighbour EXISTS. A
            // trimmed wing's strips are not spanwise-contiguous (the row
            // stops at the body flank and resumes on the far side), and
            // differencing across that gap injected spurious mid-span
            // trailing vorticity strong enough to diverge the normal-plate
            // run on the real geometry. Every segment boundary without a
            // geometric neighbour closes on zero, exactly like a tip.
            const double gLeft = adjacentToLeft[i] ? bufGamma[i - 1] : 0.0;
            emit((b.C[0] + b.C[3]) * 0.5, bufStep * (gLeft - bufGamma[i]));
            const bool rightOpen = (i + 1 == n) || !adjacentToLeft[i + 1];
            if (rightOpen) emit((b.C[1] + b.C[2]) * 0.5, bufStep * bufGamma[i]);

            // Local incidence for the separation decision.
            Vec3 vLoc = Vinf + particleVelocityAt(rings[i].Mid);
            const double alphaEff =
                Math::RadToDeg(std::atan2(Dot(vLoc, strips[i].LiftDir),
                                          Dot(vLoc, strips[i].ChordDir))) -
                strips[i].Alpha0Deg;
            bool separated;
            if (options.SeparationPoint)
                separated = options.SeparationPoint(strips[i].Eta, std::fabs(alphaEff)) <
                            PwLeSheddingF;
            else
                separated = std::fabs(alphaEff) > PwLeSheddingIncidenceDeg;

            // The separated strip's counter-rotating pair: boundary-layer
            // flux out the leading edge, its exact negative folded out the
            // trailing edge -- total circulation untouched.
            if (separated) {
                ++res.SheddingStrips;
                const double vt = std::min(vLoc.Norm(), PwLeFluxSpeedCap * fc.Vinf);
                const double leFlux = std::copysign(0.5 * vt * vt * dtPhys, alphaEff);
                const Vec3 spanDir = rings[i].Span.Normalized();
                const double offset = 0.3 * fc.Vinf * dtPhys;
                emit(rings[i].Mid - strips[i].ChordDir * (0.25 * strips[i].Chord) +
                         strips[i].LiftDir * offset,
                     spanDir * (leFlux * strips[i].Width));
                emit(rings[i].Te + dragDir * offset, spanDir * (-leFlux * strips[i].Width));
            }
        }
        bufGamma = gamma; // the refilled buffer carries the current circulation
        res.MaxParticles = std::max(res.MaxParticles, static_cast<int>(particles.size()));
        if (particles.size() > PwMaxParticles) return res; // Valid stays false

        // --- impulse loads (exact at any incidence; see the header note) --------
        Vec3 impulse(0, 0, 0);
        for (std::size_t i = 0; i < n; ++i)
            impulse = impulse + rings[i].SVec * gamma[i] + buffers[i].SVec * bufGamma[i];
        for (const Detail::WakeParticle& p : particles)
            impulse = impulse + Cross(p.X, p.Alpha) * 0.5;
        if (havePrev) {
            Vec3 force = (pPrev - impulse) * (fc.rho / dtPhys); // F = -rho dP/dt
            force = force - Cross(Vinf, retiredAlpha) * (0.5 * fc.rho);
            const double CL = Dot(force, liftDir) / qS;
            const double CD = Dot(force, dragDir) / qS;
            const double CY = Dot(force, sideDir) / qS;
            const double CN = force.z / qS;
            if (options.KeepHistory)
                res.History.push_back({t, CL, CD, CY, CN, static_cast<int>(particles.size())});
            if (t >= avgStart) {
                sCL += CL; sCD += CD; sCY += CY; sCN += CN;
                sCL2 += CL * CL; sCN2 += CN * CN;
                double circ = 0.0;
                for (std::size_t i = 0; i < n; ++i) circ += gamma[i] * strips[i].Width;
                sCircCL += 2.0 * circ / (fc.Vinf * ((ref.Area > 0.0) ? ref.Area : 1.0));
                ++avgCount;
                batchSum += CN;
                if (++batchFill >= batchSteps) {
                    batchMeans.push_back(batchSum / batchFill);
                    batchSum = 0.0;
                    batchFill = 0;
                }
            }
        }
        pPrev = impulse;
        havePrev = true;

        // --- convect + stretch + diffuse ---------------------------------------
        // Second-order (midpoint) convection: the velocity is evaluated at
        // the current positions, the particles take a half step, and the
        // full step uses the midpoint field. The gradient for stretching
        // and relaxation is taken at the midpoint. Pair cores are
        // symmetrized, (core_i^2 + core_j^2)/2, so momentum exchange stays
        // pairwise consistent as cores spread.
        const std::size_t np = particles.size();
        std::vector<Vec3> u1(np), xm(np), um(np);
        std::vector<Vec3> gm(3 * np);
        const auto filamentVelocity = [&](const Vec3& at) {
            Vec3 v(0, 0, 0);
            for (std::size_t j = 0; j < n; ++j) {
                v = v + ringVelocityUnit(at, rings[j]) * gamma[j];
                if (std::fabs(bufGamma[j]) > Math::Tiny)
                    v = v + bufferVelocityUnit(at, buffers[j]) * bufGamma[j];
            }
            return v;
        };
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
        for (long long ips = 0; ips < static_cast<long long>(np); ++ips) {
            const std::size_t ip = static_cast<std::size_t>(ips);
            Vec3 u = Vinf + filamentVelocity(particles[ip].X);
            for (std::size_t kq = 0; kq < np; ++kq) {
                if (kq == ip) continue;
                const double pair2 = 0.5 * (particles[ip].Core2 + particles[kq].Core2);
                u = u + Detail::ParticleVelocity(particles[ip].X, particles[kq], pair2);
            }
            u1[ip] = u;
            xm[ip] = particles[ip].X + u * (0.5 * dtPhys);
        }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
        for (long long ips = 0; ips < static_cast<long long>(np); ++ips) {
            const std::size_t ip = static_cast<std::size_t>(ips);
            Vec3 u = Vinf + filamentVelocity(xm[ip]);
            Vec3 grad[3] = {Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0)};
            for (std::size_t kq = 0; kq < np; ++kq) {
                if (kq == ip) continue;
                const double pair2 = 0.5 * (particles[ip].Core2 + particles[kq].Core2);
                u = u + Detail::ParticleVelocity(xm[ip], particles[kq], pair2);
                Vec3 g[3];
                Detail::ParticleVelocityGradient(xm[ip], particles[kq], pair2, g);
                for (int l = 0; l < 3; ++l) grad[l] = grad[l] + g[l];
            }
            um[ip] = u;
            for (int l = 0; l < 3; ++l) gm[3 * ip + l] = grad[l];
        }
        const bool relaxNow = (step % PwRelaxEvery) == PwRelaxEvery - 1;
        // The transient near-particles are one-step FILAMENT STAND-INS whose
        // strength the next absorb cancels by subtracting an exactly
        // spanwise vector. Stretching or relaxing them rotates that vector,
        // and the absorb then leaves a full-strength misaligned residue
        // injected at the trailing edge every event -- measured as the
        // deep-incidence blow-up that survived the relaxation fixes (quiet
        // flows barely rotate them, which is why the attached case stayed
        // clean). They convect and spread, nothing else.
        std::vector<char> isTransient(np, 0);
        for (std::size_t i = 0; i < n; ++i)
            if (nearIdx[i] >= 0) isTransient[static_cast<std::size_t>(nearIdx[i])] = 1;
        std::vector<Detail::WakeParticle> next = particles;
        for (std::size_t ip = 0; ip < np; ++ip) {
            next[ip].X = particles[ip].X + um[ip] * dtPhys;
            if (isTransient[ip]) {
                next[ip].Core2 = particles[ip].Core2 + 4.0 * nuEff * dtPhys;
                continue;
            }
            // Stretching in the TRANSPOSE scheme, d(alpha)/dt = (grad u)^T
            // alpha: unlike the classical (alpha . grad) u it conserves the
            // TOTAL vector strength exactly, and that is not a nicety here
            // -- the classical scheme's sum(d alpha) drift, levered by the
            // particles' positions, was measured burying the impulse loads
            // (RMS 3.8 on an O(1) mean at the normal plate). Capped.
            const Vec3& al = particles[ip].Alpha;
            const Vec3* grad = &gm[3 * ip];
            const Vec3 stretch(Dot(grad[0], al), Dot(grad[1], al), Dot(grad[2], al));
            const Vec3 physAlpha = al + stretch * dtPhys;
            Vec3 newAlpha = physAlpha;
            // Pedrizzetti relaxation: periodically blend the strength
            // toward the LOCAL vorticity direction. The local vorticity
            // MUST include the particle's own contribution,
            // omega_self = alpha / (2 pi sigma^3) for this kernel -- it
            // dominates at the particle and is parallel to alpha, so the
            // blend gently confirms a consistent particle and corrects an
            // inconsistent one. Excluding it (the first attempt) relaxed
            // every particle toward its neighbours' noise and was measured
            // INJECTING energy: plate CN 6.8 with RMS 277.
            if (relaxNow) {
                const double sigma2 = std::max(particles[ip].Core2, Math::Tiny);
                const double selfFactor =
                    1.0 / (2.0 * std::numbers::pi * sigma2 * std::sqrt(sigma2));
                const Vec3 curl(grad[1].z - grad[2].y + selfFactor * physAlpha.x,
                                grad[2].x - grad[0].z + selfFactor * physAlpha.y,
                                grad[0].y - grad[1].x + selfFactor * physAlpha.z);
                const double cn2 = curl.Norm();
                const double an = physAlpha.Norm();
                if (cn2 > Math::Tiny && an > Math::Tiny)
                    newAlpha = physAlpha * (1.0 - PwRelaxFraction) +
                               curl * (PwRelaxFraction * an / cn2);
            }
            const double cap = PwStretchCap * particles[ip].Birth;
            const double mag = newAlpha.Norm();
            if (mag > cap && mag > Math::Tiny) newAlpha = newAlpha * (cap / mag);
            next[ip].Alpha = newAlpha;
            next[ip].Core2 = particles[ip].Core2 + 4.0 * nuEff * dtPhys; // core spreading
            // Relaxation and the cap are bookkeeping edits, not physics:
            // their strength change must leave the impulse baseline with
            // them, or the differencer reads each edit as a force spike
            // (the retirement/merging lesson, applied to strengths).
            pPrev = pPrev + Cross(next[ip].X, newAlpha - physAlpha) * 0.5;
        }
        particles = std::move(next);

        // --- retire -------------------------------------------------------------
        std::vector<Detail::WakeParticle> kept;
        kept.reserve(particles.size());
        std::vector<int> remap(particles.size(), -1);
        for (std::size_t ip = 0; ip < particles.size(); ++ip) {
            const Detail::WakeParticle& p = particles[ip];
            if (Dot(p.X - rings[0].Mid, dragDir) > cutoff) {
                retiredAlpha = retiredAlpha + p.Alpha;
                pPrev = pPrev - Cross(p.X, p.Alpha) * 0.5; // the 2-D lesson: leave the
                continue;                                  // baseline as you leave the sum
            }
            remap[ip] = static_cast<int>(kept.size());
            kept.push_back(p);
        }
        particles = std::move(kept);
        for (std::size_t i = 0; i < n; ++i)
            if (nearIdx[i] >= 0) nearIdx[i] = remap[static_cast<std::size_t>(nearIdx[i])];

        // --- far-wake merging (diffusion's bookkeeping half) --------------------
        // Spreading cores overlap; overlapping same-scale particles are one
        // particle's worth of information. Greedy pairwise merge in the far
        // zone, conserving total strength, the strength-weighted centroid,
        // and the second moment (into the core). The merged pair's impulse
        // change is removed from the baseline -- the retirement lesson
        // applies to every event that edits the sum being differenced.
        if (nuEff > 0.0 && (step % PwMergeEvery) == PwMergeEvery - 1 && particles.size() > 1) {
            const double zone = PwMergeZoneChords * ref.Chord;
            for (std::size_t ia = 0; ia < particles.size(); ++ia) {
                Detail::WakeParticle& pa = particles[ia];
                if (pa.Birth < 0.0 || Dot(pa.X - rings[0].Te, dragDir) < zone) continue;
                for (std::size_t ib = ia + 1; ib < particles.size(); ++ib) {
                    Detail::WakeParticle& pb = particles[ib];
                    if (pb.Birth < 0.0 || Dot(pb.X - rings[0].Te, dragDir) < zone) continue;
                    const double pair2 = 0.5 * (pa.Core2 + pb.Core2);
                    const double d2 = Dot(pa.X - pb.X, pa.X - pb.X);
                    if (d2 > PwMergeDistanceFactor * PwMergeDistanceFactor * pair2) continue;
                    const double wa = pa.Alpha.Norm(), wb = pb.Alpha.Norm();
                    if (!(wa + wb > Math::Tiny)) continue;
                    const Vec3 oldImpulse =
                        Cross(pa.X, pa.Alpha) * 0.5 + Cross(pb.X, pb.Alpha) * 0.5;
                    const Vec3 xNew = (pa.X * wa + pb.X * wb) * (1.0 / (wa + wb));
                    const double spread =
                        (wa * Dot(pa.X - xNew, pa.X - xNew) + wb * Dot(pb.X - xNew, pb.X - xNew)) /
                        (wa + wb);
                    pa.Core2 = (wa * pa.Core2 + wb * pb.Core2) / (wa + wb) + spread;
                    pa.Alpha = pa.Alpha + pb.Alpha;
                    pa.X = xNew;
                    pa.Birth = std::max(pa.Birth, pb.Birth);
                    pPrev = pPrev + Cross(pa.X, pa.Alpha) * 0.5 - oldImpulse;
                    pb.Birth = -1.0; // absorbed
                    ++res.Merged;
                    break;
                }
            }
            std::vector<Detail::WakeParticle> alive;
            alive.reserve(particles.size());
            std::vector<int> remap2(particles.size(), -1);
            for (std::size_t ip = 0; ip < particles.size(); ++ip) {
                if (particles[ip].Birth < 0.0) continue;
                remap2[ip] = static_cast<int>(alive.size());
                alive.push_back(particles[ip]);
            }
            particles = std::move(alive);
            for (std::size_t i = 0; i < n; ++i)
                if (nearIdx[i] >= 0) nearIdx[i] = remap2[static_cast<std::size_t>(nearIdx[i])];
        }
    }

    if (avgCount == 0) return res;
    res.Valid = true;
    res.MeanCL = sCL / avgCount;
    res.MeanCD = sCD / avgCount;
    res.MeanCY = sCY / avgCount;
    res.MeanCN = sCN / avgCount;
    res.RmsCL = std::sqrt(std::max(sCL2 / avgCount - res.MeanCL * res.MeanCL, 0.0));
    res.RmsCN = std::sqrt(std::max(sCN2 / avgCount - res.MeanCN * res.MeanCN, 0.0));
    res.MeanCirculationCL = sCircCL / avgCount;
    res.Batches = static_cast<int>(batchMeans.size());
    if (batchMeans.size() >= 2) {
        double bm = 0.0;
        for (double b : batchMeans) bm += b;
        bm /= batchMeans.size();
        double bv = 0.0;
        for (double b : batchMeans) bv += (b - bm) * (b - bm);
        bv /= (batchMeans.size() - 1);
        res.MeanCN_CI = 2.0 * std::sqrt(bv / batchMeans.size());
    }
    return res;
}

} // namespace Aeolion::Solver
