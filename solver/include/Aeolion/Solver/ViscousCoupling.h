// Solver/ViscousCoupling.h
//
// Level-2 viscous-inviscid coupling: sectional lift feedback. The lattice
// stops being the authority on how much lift a section produces and becomes
// the authority on what the section SEES -- the induced flow. Per spanwise
// strip, per iteration:
//
//   1. the lattice provides alpha_eff (local velocity, kinematic + induced,
//      projected in the strip's section plane against its chord),
//   2. a 2-D section model provides cl_sect = f(alpha_eff, Re, Ma),
//   3. the target circulation follows from L' = rho * Vrel * Gamma
//      = 1/2 * rho * Vrel^2 * c * cl  =>  |Gamma_target| = 1/2 * Vrel * c * cl,
//   4. the circulation relaxes toward it,
//      Gamma^(k+1) = (1 - omega) * Gamma^k + omega * Gamma_target,
//   5. induced velocities are recomputed and the loop repeats until the
//      residual  R_i = cl_VLM,i(Gamma) - cl_sect,i(alpha_eff,i(Gamma))
//      is driven to zero.
//
// This is the nonlinear-lifting-line idea posed on the existing lattice
// (relaxed fixed point today; the residual form is written so Newton /
// quasi-Newton / Anderson acceleration can replace the update rule without
// touching anything else). The inviscid linear solve provides the starting
// circulation. Because the section model owns the lift curve, the lattice
// geometry's own camber boundary condition is superseded -- which is why
// StripSection carries the CST-derived zero-lift angle (see
// Geometry::SectionZeroLiftAngleDeg): the section model must represent the
// same camber the geometry would otherwise have supplied, or camber would
// be lost (or double-counted).
//
// The section model is the seam for the future boundary-layer solver: the
// AnalyticSectionModel below (finite lift slope, smooth stall saturation,
// Reynolds-scaled parabolic drag polar) is the placeholder occupying the
// exact interface -- SectionCoefficients f(strip, alphaEffDeg, Re, Ma) --
// that a 2-D integral-BL section solve will implement. Ma rides along in
// the signature for that future model; nothing here uses it yet.
//
// What Level 2 adds to the loads, beyond reshaping the lift: the section
// model's cd acts at each strip along its LOCAL relative wind, so the
// result finally carries PROFILE drag -- for a propeller, the dominant
// missing torque term of the inviscid solve (potential flow is
// d'Alembert-clean and cannot produce it).

#pragma once

#include "Aeolion/Lattice/Panel.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Solver/Solver.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <numbers>
#include <vector>

namespace Aeolion::Solver {

// --- physical constants ----------------------------------------------------
inline constexpr double SeaLevelKinematicViscosity = 1.46e-5; // [m^2/s], ISA
inline constexpr double SeaLevelSpeedOfSound = 340.3;         // [m/s], ISA

// --- coupling defaults -------------------------------------------------------
inline constexpr double DefaultCouplingRelaxation = 0.15; // omega in [0.05, 0.3]
inline constexpr int    DefaultCouplingMaxIterations = 200;
inline constexpr double DefaultCouplingTolerance = 1e-4;  // on max |cl residual|

// The target circulation Gamma = 1/2 * Vrel * c * cl uses the LOCAL speed,
// which itself grows with circulation -- so an errant strip can enter a
// runaway where Gamma inflates Vrel inflates Gamma. No physical section
// carries more than about this cl on its kinematic (rotation + inflow)
// dynamic pressure, so the target is capped against the KINEMATIC speed,
// which the circulation cannot touch. Breaks the feedback loop without
// affecting converged states (see theory.rst).
inline constexpr double MaxTargetSectionCl = 2.0;

/**
 * The section-plane frame and section data of one spanwise strip -- what a
 * 2-D section model needs to know about the geometry it is a section OF.
 * Aligned one-to-one with the lattice's panels (one Weissinger row per
 * strip).
 */
struct StripSection {
    Vec3 ChordDir;          ///< Unit, leading edge -> trailing edge.
    Vec3 LiftDir;           ///< Unit, suction side (positive-lift direction).
    double Chord = 0.0;     ///< [m]
    double Width = 0.0;     ///< Spanwise/radial width [m].
    double Eta = 0.0;       ///< Span/radius fraction keying the section shape (e.g. r/R).
    double Alpha0Deg = 0.0; ///< Thin-airfoil zero-lift angle of the section's camber line.
};

/** What a section model answers with, at one (alpha_eff, Re, Ma) state. */
struct SectionCoefficients {
    double cl = 0.0;
    double cd = 0.0;
};

/**
 * A 2-D viscous section solve: coefficients at this strip's effective angle
 * of attack (measured from the chord line, degrees), Reynolds and Mach
 * number. The boundary-layer section solver implements this signature.
 */
using SectionModel = std::function<SectionCoefficients(const StripSection&, double alphaEffDeg,
                                                       double Re, double Ma)>;

// ---------------------------------------------------------- analytic model
// The placeholder section model: thin-airfoil lift slope about the CST
// camber line's zero-lift angle, saturated smoothly toward +-ClMax (the
// same tanh form the retired BEMT solver's polar used -- it keeps the
// fixed point well behaved near stall instead of handing it a kink), and a
// parabolic drag polar whose cd0 scales with Reynolds number like
// turbulent skin friction, cd0 ~ Re^(-1/5). Mach is accepted and ignored
// (the lattice is incompressible; Prandtl-Glauert belongs in the future
// BL-coupled model).
//
// DEEP stall blends to flat-plate normal-force behavior: past the blend
// band the section is just a bluff plate, cn ~ PlateNormal * sin(alpha),
// so cl ~ cn cos(alpha) DECAYS toward high incidence and cd ~ cn
// sin(alpha) grows. Holding the saturated ClMax to 90 degrees instead is
// not conservative, it is unreachable: a swirl-dominated strip (a vane
// tip in the jet edge sits near 85 degrees) can never drive the lattice
// to a saturated cl at near-perpendicular flow, and the coupling residual
// wedges open (found empirically).
/** Analytic viscous section model -- the stand-in for the future BL solver. */
struct AnalyticSectionModel {
    double ClAlphaPerRad = Math::Two * std::numbers::pi;
    double ClMax = 1.2;
    double Cd0 = 0.012;              ///< At ReferenceReynolds.
    double KCd = 0.015;              ///< cd = cd0 + KCd * cl^2
    double ReferenceReynolds = 2e5;  ///< Small-propeller regime.
    double ReynoldsExponent = 0.2;   ///< cd0 * (RefRe/Re)^this
    double DeepStallStartDeg = 25.0; ///< Attached-polar validity edge (from zero lift).
    double DeepStallEndDeg = 40.0;   ///< Pure flat-plate behavior beyond.
    double PlateNormal = 1.8;        ///< Bluff-plate normal-force coefficient scale.

    [[nodiscard]] SectionCoefficients operator()(const StripSection& strip, double alphaEffDeg,
                                                 double Re, double /*Ma*/) const {
        const double alphaRad = Math::DegToRad(alphaEffDeg - strip.Alpha0Deg);
        const double reScale = (Re > 0.0) ? std::pow(ReferenceReynolds / Re, ReynoldsExponent) : 1.0;

        SectionCoefficients attached;
        attached.cl = ClMax * std::tanh(ClAlphaPerRad * alphaRad / ClMax);
        attached.cd = Cd0 * reScale + KCd * attached.cl * attached.cl;

        const double fromZeroLift = std::fabs(alphaEffDeg - strip.Alpha0Deg);
        if (fromZeroLift <= DeepStallStartDeg) return attached;

        const double cn = PlateNormal * std::sin(alphaRad);
        SectionCoefficients plate;
        plate.cl = cn * std::cos(alphaRad);
        plate.cd = Cd0 * reScale + cn * std::sin(alphaRad);

        const double weight =
            std::clamp((DeepStallEndDeg - fromZeroLift) / (DeepStallEndDeg - DeepStallStartDeg), 0.0, 1.0);
        return {weight * attached.cl + (1.0 - weight) * plate.cl,
                weight * attached.cd + (1.0 - weight) * plate.cd};
    }
};

// -------------------------------------------------------------- the solve
struct ViscousCouplingOptions {
    double Relaxation = DefaultCouplingRelaxation;
    int MaxIterations = DefaultCouplingMaxIterations;
    double Tolerance = DefaultCouplingTolerance;

    /**
     * Anderson-acceleration depth (number of residual differences kept).
     * Zero falls back to plain relaxed fixed-point iteration. The
     * acceleration is what makes strongly-coupled neighbor strips (their
     * induced velocities share trailing legs) converge where fixed-point
     * relaxation stalls -- see theory.rst.
     */
    int AndersonDepth = 4;

    /**
     * Strips whose |alpha_eff| exceeds this are excluded from the
     * CONVERGENCE measure (their circulation is still capped and
     * relaxed). Strip theory carries no cl-matching contract in reversed
     * flow -- a vane tip behind the jet's shear edge can sit near
     * -100 degrees -- and such a strip both wedges the residual open AND,
     * if it runs away, inflates the reference dynamic pressure until the
     * weighted residual of the raw state reads deceptively small.
     */
    double ResidualIncidenceLimitDeg = 75.0;

    /**
     * Optional warm start: when its size matches the panel count, the
     * iteration starts from this circulation instead of the inviscid
     * solve. This is CONTINUATION, and it matters where the section
     * response is multivalued (post-stall strips): a solve near a branch
     * boundary re-run from scratch can land on either branch under tiny
     * input changes, which shows up in an outer driver as a
     * never-decaying pass-to-pass oscillation; warm-started from the
     * previous pass it follows one branch smoothly.
     */
    std::vector<double> InitialGamma;
};

/** Converged per-strip state, for reporting, radial plots, and slipstream
 * reconstruction. */
struct StripState {
    double alphaEffDeg = 0.0; ///< From the chord line, induced flow included.
    double cl = 0.0;          ///< The section model's converged lift coefficient.
    double cd = 0.0;
    double Re = 0.0;
    double Ma = 0.0;
    double Vrel = 0.0;        ///< Local relative speed [m/s].
    double Residual = 0.0;    ///< cl_VLM - cl_sect at exit.
    Vec3 Mid{0, 0, 0};        ///< Bound-vortex midpoint (the strip's load point).
    Vec3 Force{0, 0, 0};      ///< Circulatory + profile force on the strip [N].
};

/**
 * Full result: the familiar SolveResult shape (circulation, stations,
 * dimensional force/moment projections -- so every existing consumer,
 * renderer included, reads it unchanged) plus the coupling's own state.
 * Moments split so a propeller consumer can report induced and profile
 * torque separately.
 */
struct ViscousCoupledResult {
    SolveResult Base;
    std::vector<StripState> Strips;
    Vec3 InducedMoment{0, 0, 0}; ///< From the circulation forces, about RefPoint [N*m].
    Vec3 ProfileMoment{0, 0, 0}; ///< From the section-drag forces, about RefPoint [N*m].
    Vec3 SourceForce{0, 0, 0};   ///< Pressure force on the source panels (e.g. a duct shroud) [N].
    Vec3 SourceMoment{0, 0, 0};  ///< Its moment about RefPoint [N*m].
    std::vector<double> sigma;   ///< Converged source strengths, aligned with `sources`.
    int Iterations = 0;
    bool Converged = false;
    double MaxResidual = 0.0;
};

/**
 * Sectional lift-feedback solve over a single-row lattice. `strips` must
 * align one-to-one with `panels`. The inviscid linear solve seeds the
 * circulation; the accelerated fixed point then drives
 * R_i = cl_VLM,i - cl_sect,i to zero.
 *
 * `sources` (optional) is a source-panel body sharing the solve -- for a
 * propeller, the duct shroud. The interaction is two-way: every coupling
 * iteration re-solves the source strengths against the CURRENT blade
 * circulation (their boundary condition sees the propwash), and the strip
 * velocities include the sources' induced flow (the blades see the duct).
 * The sources' own loads are the same pressure integral the airframe
 * solve uses, reported in SourceForce/SourceMoment and folded into the
 * Base totals. In the rotating frame this is consistent for bodies of
 * revolution about +x: the frame's rotation moves such a surface only
 * tangentially to itself, so its no-through-flow condition is unchanged
 * (up to faceting).
 */
[[nodiscard]] inline ViscousCoupledResult SolveViscousCoupled(
    const std::vector<Panel>& panels, const std::vector<StripSection>& strips,
    const FreestreamConditions& fc, const ReferenceGeometry& ref, double trail,
    const SectionModel& model, const ViscousCouplingOptions& options = {},
    const std::vector<SourcePanel>& sources = {},
    const std::function<Vec3(const Vec3&)>& externalField = nullptr) {
    const std::size_t n = panels.size();
    const std::size_t nb = sources.size();
    ViscousCoupledResult res;
    if (n == 0 || strips.size() != n) return res;

    // Seed with the inviscid solution (blades and sources coupled when
    // sources are present) -- the lattice's own boundary condition is the
    // right starting point for the fixed point -- unless the caller
    // supplied a warm start (continuation across an outer iteration; see
    // ViscousCouplingOptions::InitialGamma).
    const PanelSystem system{panels, sources};
    std::vector<double> gamma;
    if (options.InitialGamma.size() == panels.size()) {
        gamma = options.InitialGamma;
    } else if (nb > 0) {
        gamma = SolveWithSystem(Prepare(system, trail), fc, ref, externalField).gamma;
    } else {
        gamma = Solve(panels, fc, ref, trail, externalField).gamma;
    }

    // Source-only influence block, factored once: geometry-only, reused
    // every iteration to re-solve sigma under the current circulation.
    LUFactorization sourceLu;
    if (nb > 0) {
        std::vector<double> block(nb * nb);
        DenseMatrixView view(block.data(), static_cast<int>(nb), static_cast<int>(nb));
        const int offset = system.VortexCount();
        for (std::size_t i = 0; i < nb; ++i)
            for (std::size_t j = 0; j < nb; ++j)
                view[static_cast<int>(i), static_cast<int>(j)] =
                    system.NormalInfluence(offset + static_cast<int>(i), offset + static_cast<int>(j), trail);
        sourceLu = LuFactorize(view);
    }
    std::vector<double>& sigma = res.sigma;
    sigma.assign(nb, 0.0);

    const double alpha = Math::DegToRad(fc.alphaDeg);
    const double beta = Math::DegToRad(fc.betaDeg);
    const Vec3 Vinf(fc.Vinf * std::cos(alpha) * std::cos(beta), fc.Vinf * std::sin(beta),
                    fc.Vinf * std::sin(alpha) * std::cos(beta));
    const Vec3 omegaRate(fc.p, fc.q, fc.r);
    const auto kinematicVelocity = [&](const Vec3& point) {
        Vec3 v = Vinf - Cross(omegaRate, point - fc.RefPoint);
        if (externalField) v = v + externalField(point);
        return v;
    };

    // Velocity induced by the source body at an arbitrary point, under the
    // current sigma. Empty sources collapse it to zero.
    const auto sourceVelocity = [&](const Vec3& point) {
        Vec3 v(0, 0, 0);
        for (std::size_t k = 0; k < nb; ++k)
            v = v + SourcePanelVelocity(point, sources[k]) * sigma[k];
        return v;
    };

    // Re-solve the source strengths against the current circulation: the
    // body's boundary condition sees the blades' full induced flow.
    const auto updateSources = [&]() {
        if (nb == 0) return;
        std::vector<double> rhs(nb);
        for (std::size_t i = 0; i < nb; ++i) {
            Vec3 v = kinematicVelocity(sources[i].ControlPoint);
            for (std::size_t j = 0; j < n; ++j)
                v = v + HorseshoeVelocity(sources[i].ControlPoint, panels[j], gamma[j], trail);
            rhs[i] = sources[i].PrescribedNormalVelocity - Dot(v, sources[i].Normal);
        }
        sigma = LuSolve(sourceLu, rhs);
    };

    // Local velocity at a strip's bound-vortex midpoint under the CURRENT
    // circulation and source strengths -- own bound vortex excluded
    // (singular there), exactly the near-field force evaluation's
    // convention.
    const auto localVelocity = [&](std::size_t i) {
        const Vec3 mid = (panels[i].A + panels[i].B) * Math::Half;
        Vec3 v = kinematicVelocity(mid) + sourceVelocity(mid);
        for (std::size_t j = 0; j < n; ++j) {
            v = v + ((j == i) ? HorseshoeVelocityNoBound(mid, panels[j], gamma[j], trail)
                              : HorseshoeVelocity(mid, panels[j], gamma[j], trail));
        }
        return v;
    };

    res.Strips.resize(n);

    // Kinematic-speed circulation caps, one per strip (see MaxTargetSectionCl).
    std::vector<double> gammaCap(n);
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 mid = (panels[i].A + panels[i].B) * Math::Half;
        gammaCap[i] = Math::Half * kinematicVelocity(mid).Norm() * strips[i].Chord * MaxTargetSectionCl;
    }

    // --- Anderson-accelerated fixed point on the circulation ----------------
    // The map is G(Gamma) = Gamma_target(Gamma); its residual is
    // f = G(Gamma) - Gamma. Type-II Anderson mixing over the last few
    // residual differences (Walker & Ni) replaces the plain relaxed
    // update: neighboring strips couple strongly through their shared
    // trailing legs, and the resulting stiff modes stall a damped fixed
    // point that Anderson converges. Safeguard: if the residual norm
    // doubles, the history is discarded and a damped plain step taken.
    // Full write-up in theory.rst.
    double omegaCoupling = options.Relaxation;
    double bestResidual = 1e30;
    double fNormPrev = 1e30;
    std::vector<std::vector<double>> historyX, historyF; // iterates and map residuals
    // Tail mean of the circulation: when the iteration exhausts its budget
    // in a limit cycle (deep post-stall strips do this), the returned state
    // is the CYCLE MEAN -- the arithmetic mean over the second half of the
    // budget (several full cycles, no memory of the transient), evaluated
    // consistently in one final sweep. The mean is reproducible across
    // repeated solves, which is what lets an outer driver (the rotor-vane
    // coupling) converge on stable quasi-steady loads instead of sampling
    // a random phase of the cycle. Converged stays false and MaxResidual
    // keeps the cycle's mismatch: the mean is reported, not declared
    // converged.
    std::vector<double> gammaMeanSum(n, 0.0);
    int gammaMeanCount = 0;
    bool finalSweep = false;
    for (res.Iterations = 1; res.Iterations <= options.MaxIterations; ++res.Iterations) {
        updateSources();
        res.MaxResidual = 0.0;
        std::vector<double> target(n, 0.0);
        std::vector<double> rawResidual(n, 0.0);
        double vrelMaxSq = Math::Tiny;

        for (std::size_t i = 0; i < n; ++i) {
            const StripSection& strip = strips[i];
            const Vec3 v = localVelocity(i);
            const double vrel = v.Norm();
            const double alphaEffDeg =
                Math::RadToDeg(std::atan2(Dot(v, strip.LiftDir), Dot(v, strip.ChordDir)));
            const double Re = vrel * strip.Chord / SeaLevelKinematicViscosity;
            const double Ma = vrel / SeaLevelSpeedOfSound;

            const SectionCoefficients sect = model(strip, alphaEffDeg, Re, Ma);

            // cl_VLM is linear in this strip's gamma through Kutta-Joukowski
            // projected on the lift direction:
            //   L'_i * w = rho * gamma_i * (v x dl) . LiftDir
            //   cl_VLM,i = gamma_i * k_i,  k_i = 2 (v x dl) . LiftDir / (Vrel^2 c w).
            // Inverting k for the target keeps every orientation convention
            // (dl sense, rotating frame, camber tilt) out of the update --
            // no hand-derived sign anywhere.
            const Vec3 dl = panels[i].B - panels[i].A;
            const double denom = vrel * vrel * strip.Chord * strip.Width;
            const double k = (denom > Math::Tiny) ? 2.0 * Dot(Cross(v, dl), strip.LiftDir) / denom : 0.0;
            const double clVlm = gamma[i] * k;

            const double residual = clVlm - sect.cl;
            rawResidual[i] = residual;
            vrelMaxSq = std::max(vrelMaxSq, vrel * vrel);
            target[i] = (std::fabs(k) > Math::Tiny)
                            ? std::clamp(sect.cl / k, -gammaCap[i], gammaCap[i])
                            : gamma[i];

            StripState& state = res.Strips[i];
            state.alphaEffDeg = alphaEffDeg;
            state.cl = sect.cl;
            state.cd = sect.cd;
            state.Re = Re;
            state.Ma = Ma;
            state.Vrel = vrel;
            state.Residual = residual;
        }

        // The convergence measure is the DYNAMIC-PRESSURE-WEIGHTED residual:
        // a cl mismatch on a strip carrying a fraction of the reference
        // dynamic pressure is the same fraction of a force mismatch --
        // and it is taken only over strips inside the incidence limit
        // (see ViscousCouplingOptions::ResidualIncidenceLimitDeg):
        // reversed-flow strips have no cl-matching contract to converge.
        const auto eligible = [&](std::size_t i) {
            return std::fabs(res.Strips[i].alphaEffDeg) <= options.ResidualIncidenceLimitDeg;
        };
        double eligibleVrelMaxSq = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            if (eligible(i))
                eligibleVrelMaxSq =
                    std::max(eligibleVrelMaxSq, res.Strips[i].Vrel * res.Strips[i].Vrel);
        const bool anyEligible = eligibleVrelMaxSq > Math::Tiny;
        if (!anyEligible) eligibleVrelMaxSq = vrelMaxSq;
        for (std::size_t i = 0; i < n; ++i) {
            if (anyEligible && !eligible(i)) continue;
            res.MaxResidual = std::max(
                res.MaxResidual, std::fabs(rawResidual[i]) *
                                     std::min(res.Strips[i].Vrel * res.Strips[i].Vrel /
                                                  eligibleVrelMaxSq,
                                              1.0));
        }

        if (finalSweep) break; // strips/residual now consistent with the cycle mean
        if (res.MaxResidual < options.Tolerance) {
            res.Converged = true;
            break;
        }
        if (res.MaxResidual > bestResidual)
            omegaCoupling = std::max(Math::Half * omegaCoupling, 0.02);
        bestResidual = std::min(bestResidual, res.MaxResidual);

        // Map residual f = G(x) - x and its norm.
        std::vector<double> f(n);
        double fNorm = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            f[i] = target[i] - gamma[i];
            fNorm += f[i] * f[i];
        }
        fNorm = std::sqrt(fNorm);

        const bool diverging = fNorm > 2.0 * fNormPrev;
        fNormPrev = fNorm;
        if (diverging) {
            historyX.clear();
            historyF.clear();
        }

        const std::size_t depth = historyX.size();
        if (options.AndersonDepth > 0 && depth >= 1 && !diverging) {
            // Least squares min || f_k - dF theta || over the history of
            // residual DIFFERENCES, via the (tiny) normal equations.
            const std::size_t m = depth;
            std::vector<std::vector<double>> dF(m, std::vector<double>(n));
            std::vector<std::vector<double>> dX(m, std::vector<double>(n));
            for (std::size_t j = 0; j < m; ++j)
                for (std::size_t i = 0; i < n; ++i) {
                    dF[j][i] = f[i] - historyF[j][i];
                    dX[j][i] = gamma[i] - historyX[j][i];
                }
            std::vector<double> normal(m * m, 0.0), rhsLs(m, 0.0);
            for (std::size_t a = 0; a < m; ++a) {
                for (std::size_t b = 0; b < m; ++b)
                    for (std::size_t i = 0; i < n; ++i) normal[a * m + b] += dF[a][i] * dF[b][i];
                normal[a * m + a] += 1e-12; // Tikhonov guard for collinear history
                for (std::size_t i = 0; i < n; ++i) rhsLs[a] += dF[a][i] * f[i];
            }
            // Gaussian elimination on the m x m system (m <= AndersonDepth).
            for (std::size_t col = 0; col < m; ++col) {
                std::size_t best = col;
                for (std::size_t row = col + 1; row < m; ++row)
                    if (std::fabs(normal[row * m + col]) > std::fabs(normal[best * m + col])) best = row;
                for (std::size_t k = 0; k < m; ++k) std::swap(normal[col * m + k], normal[best * m + k]);
                std::swap(rhsLs[col], rhsLs[best]);
                const double diag = normal[col * m + col];
                if (std::fabs(diag) < 1e-30) continue;
                for (std::size_t row = col + 1; row < m; ++row) {
                    const double factor = normal[row * m + col] / diag;
                    for (std::size_t k = col; k < m; ++k) normal[row * m + k] -= factor * normal[col * m + k];
                    rhsLs[row] -= factor * rhsLs[col];
                }
            }
            std::vector<double> theta(m, 0.0);
            for (std::size_t row = m; row-- > 0;) {
                double sum = rhsLs[row];
                for (std::size_t k = row + 1; k < m; ++k) sum -= normal[row * m + k] * theta[k];
                const double diag = normal[row * m + row];
                theta[row] = (std::fabs(diag) > 1e-30) ? sum / diag : 0.0;
            }

            // Anderson update, then the physical circulation cap.
            historyX.push_back(gamma);
            historyF.push_back(f);
            for (std::size_t i = 0; i < n; ++i) {
                double next = gamma[i] + f[i];
                for (std::size_t j = 0; j < m; ++j) next -= theta[j] * (dX[j][i] + dF[j][i]);
                gamma[i] = std::clamp(next, -gammaCap[i], gammaCap[i]);
            }
        } else {
            historyX.push_back(gamma);
            historyF.push_back(f);
            for (std::size_t i = 0; i < n; ++i)
                gamma[i] = std::clamp(gamma[i] + omegaCoupling * f[i], -gammaCap[i], gammaCap[i]);
        }
        while (historyX.size() > static_cast<std::size_t>(std::max(options.AndersonDepth, 1))) {
            historyX.erase(historyX.begin());
            historyF.erase(historyF.begin());
        }

        if (2 * res.Iterations > options.MaxIterations) {
            for (std::size_t i = 0; i < n; ++i) gammaMeanSum[i] += gamma[i];
            ++gammaMeanCount;
        }
        if (res.Iterations == options.MaxIterations - 1 && gammaMeanCount > 0) {
            // Within the caps: a convex mix of capped iterates.
            for (std::size_t i = 0; i < n; ++i)
                gamma[i] = gammaMeanSum[i] / static_cast<double>(gammaMeanCount);
            finalSweep = true;
        }
    }

    // --- loads under the converged circulation ------------------------------
    // Kutta-Joukowski with the full local velocity (circulation forces:
    // lift + induced drag), plus the section model's profile drag along the
    // local relative wind. Wind-axis projections match SolveResult's so
    // consumers read this like any solve.
    const Vec3 dragDir = Vinf.Normalized();
    const Vec3 liftDir = Cross(dragDir, Vec3(0, 1, 0)).Normalized();
    const Vec3 sideDir = Cross(liftDir, dragDir).Normalized();
    updateSources(); // sigma current for the final circulation

    Vec3 totalForce(0, 0, 0);
    res.Base.gamma = gamma;
    res.Base.Stations.reserve(n);
    double areaSum = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        const StripSection& strip = strips[i];
        StripState& state = res.Strips[i];
        const Vec3 mid = (panels[i].A + panels[i].B) * Math::Half;
        const Vec3 v = localVelocity(i);

        const Vec3 circulatory = Cross(v, panels[i].B - panels[i].A) * (fc.rho * gamma[i]);
        const double q = Math::Half * fc.rho * state.Vrel * state.Vrel;
        const Vec3 profile = v.Normalized() * (q * strip.Chord * strip.Width * state.cd);
        state.Mid = mid;
        state.Force = circulatory + profile;

        totalForce = totalForce + circulatory + profile;
        res.InducedMoment = res.InducedMoment + Cross(mid - fc.RefPoint, circulatory);
        res.ProfileMoment = res.ProfileMoment + Cross(mid - fc.RefPoint, profile);

        StationResult sr;
        sr.y = mid.y;
        sr.Surface = panels[i].Surface;
        sr.PanelIndex = static_cast<int>(i);
        sr.gamma = gamma[i];
        sr.LiftPerSpan = (strip.Width > Math::Tiny) ? Dot(circulatory, strip.LiftDir) / strip.Width : 0.0;
        sr.cl_local = state.cl;
        res.Base.Stations.push_back(sr);

        res.Base.LiftBySurface[panels[i].Surface] += Dot(circulatory + profile, liftDir);
        res.Base.DragBySurface[panels[i].Surface] += Dot(circulatory + profile, dragDir);
        res.Base.AreaBySurface[panels[i].Surface] += panels[i].PlanformArea;
        areaSum += panels[i].PlanformArea;
    }

    // Source-body pressure loads under the converged state -- identical in
    // form to the airframe solve's body-pressure integral, with the local
    // velocity carrying the blades' full induced flow (this is where the
    // duct's lip-suction thrust appears).
    const double qInf = Math::Half * fc.rho * fc.Vinf * fc.Vinf;
    for (std::size_t k = 0; k < nb; ++k) {
        const SourcePanel& body = sources[k];
        if (body.Permeable) continue; // a transpiring face has no wall to press on
        Vec3 v = kinematicVelocity(body.ControlPoint);
        for (std::size_t j = 0; j < n; ++j)
            v = v + HorseshoeVelocity(body.ControlPoint, panels[j], gamma[j], trail);
        for (std::size_t other = 0; other < nb; ++other) {
            v = v + ((other == k)
                         ? body.Normal * (sigma[other] * SourceSelfInfluence)
                         : SourcePanelVelocity(body.ControlPoint, sources[other]) * sigma[other]);
        }
        // Gauge pressure p - p_inf = qInf - 1/2 rho |V|^2, acting inward:
        // the outward force is (1/2 rho |V|^2 - qInf) * A * n. Identical to
        // the airframe solve's cp form, and well-behaved at the hover floor
        // speed where qInf is negligible.
        const Vec3 F = body.Normal * ((Math::Half * fc.rho * Dot(v, v) - qInf) * body.Area);
        res.SourceForce = res.SourceForce + F;
        res.SourceMoment = res.SourceMoment + Cross(body.ControlPoint - fc.RefPoint, F);
        res.Base.LiftBySurface[body.Surface] += Dot(F, liftDir);
        res.Base.DragBySurface[body.Surface] += Dot(F, dragDir);
    }
    totalForce = totalForce + res.SourceForce;

    const Vec3 totalMoment = res.InducedMoment + res.ProfileMoment + res.SourceMoment;
    res.Base.L = Dot(totalForce, liftDir);
    res.Base.Di = Dot(totalForce, dragDir);
    res.Base.Y = Dot(totalForce, sideDir);
    res.Base.Mx = totalMoment.x;
    res.Base.My = totalMoment.y;
    res.Base.Mz = totalMoment.z;
    res.Base.ReferenceArea = (ref.Area > 0.0) ? ref.Area : areaSum;
    res.Base.ReferenceChord = ref.Chord;
    res.Base.ReferenceSpan = ref.Span;

    return res;
}

// ------------------------------------------------------- slipstream field
// The time-mean propwash a STATIC surface downstream of the rotor sees,
// reconstructed from the CONVERGED radial load distribution by annular
// momentum theory -- not by averaging the lattice's own induced field. The
// distinction matters: the quasi-steady solve's prescribed straight
// trailing legs never wrap into the downstream helical wake cylinder, so
// their Biot-Savart field carries almost no axial jet a few chord lengths
// aft of the disk (verified empirically: vanes placed there saw swirl but
// no jet). Momentum theory applied to the loads the lattice actually
// solved is consistent by construction:
//
//   per annulus:  dT = 4 pi r rho (Vax + vi) vi dr
//                     =>  vi = ( -Vax + sqrt(Vax^2 + dT/(pi r rho dr)) ) / 2
//                 dQ = 4 pi r^3 rho (Vax + vi) wi dr   =>  wi
//
// with dT, dQ summed over the blades from the per-strip forces. Downstream
// the axial component develops from its at-disk value toward the classical
// far-wake doubling over about a diameter; swirl is carried unchanged (the
// same engineering shape the retired BEMT project's SlipstreamField used).
// The duct's contribution rides along exactly (an axisymmetric static
// source body needs no averaging). The returned field is the induced
// PERTURBATION only -- the shape of the externalField hook -- and the
// coupling is one-way: the vanes' influence back on the rotor is not
// represented.
inline constexpr double SlipstreamFarWakeDiameters = 2.0; // development length = this * radius

/** One annular band of the collapsed rotor loading, with its
 * momentum-theory induced velocities. */
struct SlipstreamBand {
    double r = 0.0, width = 0.0, dT = 0.0, dQ = 0.0, vi = 0.0, wi = 0.0;
};

/** The banded rotor loading (sorted by radius) and its tip extent. */
struct SlipstreamBands {
    std::vector<SlipstreamBand> Bands;
    double TipRadius = 0.0;
};

/**
 * Collapse a converged rotor solve's strips (all blades) into radial
 * bands of dT and dQ, and solve each band's annular momentum balance for
 * the induced velocities vi (axial) and wi (swirl). This is the shared
 * core of the slipstream reconstruction AND the self-consistent wake
 * pitch (theory.rst).
 */
[[nodiscard]] inline SlipstreamBands ComputeSlipstreamBands(const ViscousCoupledResult& rotor,
                                                            double axialSpeed, double rho) {
    SlipstreamBands result;
    std::vector<SlipstreamBand>& bands = result.Bands;
    for (const StripState& strip : rotor.Strips) {
        const double r = std::hypot(strip.Mid.y, strip.Mid.z);
        const double thrust = -strip.Force.x;
        const double torque = -(strip.Mid.y * strip.Force.z - strip.Mid.z * strip.Force.y);
        bool merged = false;
        for (SlipstreamBand& band : bands) {
            if (std::fabs(band.r - r) < 1e-9) {
                band.dT += thrust;
                band.dQ += torque;
                merged = true;
                break;
            }
        }
        if (!merged) bands.push_back({r, 0.0, thrust, torque, 0.0, 0.0});
    }
    std::sort(bands.begin(), bands.end(),
              [](const SlipstreamBand& a, const SlipstreamBand& b) { return a.r < b.r; });
    for (std::size_t i = 0; i < bands.size(); ++i) {
        SlipstreamBand& band = bands[i];
        const double below = (i > 0) ? bands[i - 1].r : band.r;
        const double above = (i + 1 < bands.size()) ? bands[i + 1].r : band.r;
        band.width = std::max(Math::Half * (above - below) + ((i == 0 || i + 1 == bands.size())
                                                                  ? Math::Half * (above - below)
                                                                  : 0.0),
                              1e-6);
        const double denomT = std::numbers::pi * band.r * rho * band.width;
        const double c = std::max(band.dT, 0.0) / std::max(denomT, Math::Tiny);
        band.vi = Math::Half * (-axialSpeed + std::sqrt(axialSpeed * axialSpeed + c));
        const double denomQ = 4.0 * std::numbers::pi * band.r * band.r * band.r * rho *
                              std::max(axialSpeed + band.vi, 1e-6) * band.width;
        band.wi = band.dQ / std::max(denomQ, Math::Tiny);
        result.TipRadius = std::max(result.TipRadius, band.r + Math::Half * band.width);
    }
    return result;
}

/**
 * The axial induced-velocity distribution vi(r) of banded rotor loading,
 * as a callable (linear interpolation between band centers, zero outside
 * the disk) -- the wake-pitch input BuildPropellerLattice consumes for
 * the self-consistent wake.
 */
[[nodiscard]] inline std::function<double(double)> AxialInflowFromBands(SlipstreamBands bands) {
    return [bands = std::move(bands)](double r) {
        const std::vector<SlipstreamBand>& b = bands.Bands;
        if (b.empty() || r > bands.TipRadius || r < 0.0) return 0.0;
        std::size_t hi = 0;
        while (hi + 1 < b.size() && b[hi + 1].r < r) ++hi;
        const std::size_t lo = hi;
        if (hi + 1 < b.size()) ++hi;
        const double span = b[hi].r - b[lo].r;
        const double t = (span > Math::Tiny) ? std::clamp((r - b[lo].r) / span, 0.0, 1.0) : 0.0;
        return b[lo].vi + t * (b[hi].vi - b[lo].vi);
    };
}

/** Time-mean slipstream of a converged rotor(+duct) solve, for downstream
 * static surfaces. axialSpeed is the freestream through the disk.
 * swirlFactor scales the tangential component -- the two-way coupling's
 * momentum-budget handle (see theory.rst, "Two-way rotor-vane
 * coupling"). */
[[nodiscard]] inline std::function<Vec3(const Vec3&)> SlipstreamField(
    const ViscousCoupledResult& rotor, double axialSpeed, double rho,
    std::vector<SourcePanel> sources = {}, double swirlFactor = 1.0) {
    SlipstreamBands banded = ComputeSlipstreamBands(rotor, axialSpeed, rho);
    const double tipRadius = banded.TipRadius;

    return [bands = std::move(banded.Bands), sources = std::move(sources), sigma = rotor.sigma,
            tipRadius, swirlFactor](const Vec3& point) {
        Vec3 v(0, 0, 0);
        // The duct's induced flow, exact (static, axisymmetric body).
        for (std::size_t m = 0; m < sources.size() && m < sigma.size(); ++m)
            v = v + SourcePanelVelocity(point, sources[m]) * sigma[m];

        const double r = std::hypot(point.y, point.z);
        if (bands.empty() || point.x < 0.0 || r > tipRadius || r < 1e-9) return v;

        // Interpolate vi, wi radially between band centers.
        std::size_t hi = 0;
        while (hi + 1 < bands.size() && bands[hi + 1].r < r) ++hi;
        const std::size_t lo = hi;
        if (hi + 1 < bands.size()) ++hi;
        const double span = bands[hi].r - bands[lo].r;
        const double t = (span > Math::Tiny) ? std::clamp((r - bands[lo].r) / span, 0.0, 1.0) : 0.0;
        const double vi = bands[lo].vi + t * (bands[hi].vi - bands[lo].vi);
        const double wi = bands[lo].wi + t * (bands[hi].wi - bands[lo].wi);

        // Axial development toward the far-wake doubling; swirl carried.
        const double devLength = SlipstreamFarWakeDiameters * tipRadius;
        const double development = std::clamp(point.x / devLength, 0.0, 1.0);
        const Vec3 rhat(0.0, point.y / r, point.z / r);
        const Vec3 that = Cross(Vec3(1.0, 0.0, 0.0), rhat); // +Omega sense about +x
        return v + Vec3(1.0, 0.0, 0.0) * (vi * (1.0 + development)) + that * (wi * swirlFactor);
    };
}

// The azimuthal-mean induced field of a STATIC panel system, as seen by
// the rotating frame: a blade sweeps past the static system once per
// revolution, and its quasi-steady boundary condition sees the time mean
// of that encounter,
//
//     v(P) = (1/K) sum_j R_x(-theta_j) . v( R_x(theta_j) . P ).
//
// This is the vane-to-rotor half of the two-way coupling (theory.rst),
// and it is legitimate in THIS direction because the static system's
// near field is fully represented -- its wakes are explicit legs, unlike
// the rotor's prescribed wake, which is why the rotor-to-vane direction
// uses the momentum reconstruction above instead.
inline constexpr int DefaultAzimuthSamples = 16;

/** Azimuthal-mean induced field of static panels (e.g. the vane system). */
[[nodiscard]] inline std::function<Vec3(const Vec3&)> MeanInducedField(
    std::vector<Panel> panels, std::vector<double> gamma, double trail,
    int azimuthSamples = DefaultAzimuthSamples) {
    return [panels = std::move(panels), gamma = std::move(gamma), trail,
            azimuthSamples](const Vec3& point) {
        Vec3 mean(0, 0, 0);
        const int samples = std::max(azimuthSamples, 1);
        for (int j = 0; j < samples; ++j) {
            const double theta = Math::Two * std::numbers::pi * j / samples;
            const Vec3 rotated = RotateAboutX(point, theta);
            Vec3 v(0, 0, 0);
            for (std::size_t i = 0; i < panels.size() && i < gamma.size(); ++i)
                v = v + HorseshoeVelocity(rotated, panels[i], gamma[i], trail);
            mean = mean + RotateAboutX(v, -theta);
        }
        return mean * (1.0 / samples);
    };
}

} // namespace Aeolion::Solver
