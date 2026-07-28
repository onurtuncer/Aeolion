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
/** Analytic viscous section model -- the stand-in for the future BL solver. */
struct AnalyticSectionModel {
    double ClAlphaPerRad = Math::Two * std::numbers::pi;
    double ClMax = 1.2;
    double Cd0 = 0.012;              ///< At ReferenceReynolds.
    double KCd = 0.015;              ///< cd = cd0 + KCd * cl^2
    double ReferenceReynolds = 2e5;  ///< Small-propeller regime.
    double ReynoldsExponent = 0.2;   ///< cd0 * (RefRe/Re)^this

    [[nodiscard]] SectionCoefficients operator()(const StripSection& strip, double alphaEffDeg,
                                                 double Re, double /*Ma*/) const {
        const double alphaRad = Math::DegToRad(alphaEffDeg - strip.Alpha0Deg);
        SectionCoefficients c;
        c.cl = ClMax * std::tanh(ClAlphaPerRad * alphaRad / ClMax);
        const double reScale = (Re > 0.0) ? std::pow(ReferenceReynolds / Re, ReynoldsExponent) : 1.0;
        c.cd = Cd0 * reScale + KCd * c.cl * c.cl;
        return c;
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
};

/** Converged per-strip state, for reporting and radial plots. */
struct StripState {
    double alphaEffDeg = 0.0; ///< From the chord line, induced flow included.
    double cl = 0.0;          ///< The section model's converged lift coefficient.
    double cd = 0.0;
    double Re = 0.0;
    double Ma = 0.0;
    double Vrel = 0.0;        ///< Local relative speed [m/s].
    double Residual = 0.0;    ///< cl_VLM - cl_sect at exit.
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
    const std::vector<SourcePanel>& sources = {}) {
    const std::size_t n = panels.size();
    const std::size_t nb = sources.size();
    ViscousCoupledResult res;
    if (n == 0 || strips.size() != n) return res;

    // Seed with the inviscid solution (blades and sources coupled when
    // sources are present) -- the lattice's own boundary condition is the
    // right starting point for the fixed point.
    const PanelSystem system{panels, sources};
    std::vector<double> gamma;
    if (nb > 0) {
        gamma = SolveWithSystem(Prepare(system, trail), fc, ref).gamma;
    } else {
        gamma = Solve(panels, fc, ref, trail).gamma;
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
        return Vinf - Cross(omegaRate, point - fc.RefPoint);
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
    for (res.Iterations = 1; res.Iterations <= options.MaxIterations; ++res.Iterations) {
        updateSources();
        res.MaxResidual = 0.0;
        std::vector<double> target(n, 0.0);

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
            res.MaxResidual = std::max(res.MaxResidual, std::fabs(residual));
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
        const StripState& state = res.Strips[i];
        const Vec3 mid = (panels[i].A + panels[i].B) * Math::Half;
        const Vec3 v = localVelocity(i);

        const Vec3 circulatory = Cross(v, panels[i].B - panels[i].A) * (fc.rho * gamma[i]);
        const double q = Math::Half * fc.rho * state.Vrel * state.Vrel;
        const Vec3 profile = v.Normalized() * (q * strip.Chord * strip.Width * state.cd);

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

} // namespace Aeolion::Solver
