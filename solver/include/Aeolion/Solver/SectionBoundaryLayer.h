// Solver/SectionBoundaryLayer.h
//
// The 2-D viscous section solver behind the Level-2 coupling's SectionModel
// seam (ViscousCoupling.h): a transpiration-coupled integral boundary layer
// over a lumped-vortex model of the section's camber line. This is Level-3
// viscous-inviscid interaction at the SECTION level -- the boundary layer
// changes the section's own inviscid problem, not just its polar -- while
// the 3-D lattice keeps talking to it through the same
// f(alpha_eff, Re, Ma) -> {cl, cd} interface as the analytic model.
//
// The full numerical method -- discretization, boundary-layer correlations,
// the transpiration boundary condition and its collapse onto the camber
// line, the iteration strategy, and every limitation -- is documented in
// doc/theory.rst ("Level-3: transpiration-coupled boundary layer"). The
// comments here only anchor code to that writeup.
//
// In brief:
//   - The section is a 2-D lumped-vortex method on the camber line
//     (cosine-spaced panels, point vortex at each panel quarter point,
//     flow tangency at each 3/4 point) -- the exact 2-D analog of the
//     3-D lattice, so the two levels agree by construction in the
//     inviscid, thin limit.
//   - The displacement effect enters as TRANSPIRATION: the equivalent
//     inviscid flow satisfies v_n = d(U_e * delta*)/ds on each true
//     surface; collapsed onto the camber line, the vortex sheet's
//     tangency condition receives the antisymmetric part
//        v_t(s) = 1/2 * d/ds [ U_e^u delta*_u  -  U_e^l delta*_l ].
//     The geometry never moves -- no recambering -- the boundary
//     condition carries the whole effect.
//   - The boundary layer is a classical integral march per surface:
//     Thwaites (laminar) -> Michel's criterion or laminar separation
//     (transition) -> Head's entrainment method with Ludwieg-Tillmann
//     skin friction (turbulent), Squire-Young for the drag.
//   - The section iterates: solve vortices -> edge velocities -> march
//     boundary layers -> new transpiration -> relax -> repeat.
//
// Everything is nondimensionalized by the chord and the relative speed;
// the Reynolds number carries the physics in. Mach is accepted and unused
// (incompressible), as everywhere else in the toolkit.

#pragma once

#include "Aeolion/Math/Constants.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numbers>
#include <vector>

namespace Aeolion::Solver {

// --- discretization ---------------------------------------------------------
inline constexpr int DefaultSectionStations = 40;   // camber-line panels
inline constexpr int CamberQuadratureSubsamples = 2; // Simpson half-intervals per panel for z(psi)

// --- transpiration iteration -------------------------------------------------
inline constexpr int    DefaultSectionMaxIterations = 100;
inline constexpr double DefaultSectionTolerance = 1e-3;   // on |delta cl|
inline constexpr double DefaultSectionRelaxation = 0.3;   // on the transpiration distribution
inline constexpr int    SectionOutputAverageWindow = 48;  // phase-averaging window (see theory.rst)
inline constexpr int    SectionConvergenceStreak = 3;     // consecutive small deltas required: a cycle
                                                          // crossing itself fakes single-delta convergence

// --- boundary-layer correlations (see theory.rst for provenance) ------------
inline constexpr double ThwaitesConstant = 0.45;
inline constexpr double LaminarSeparationLambda = -0.09;
inline constexpr double MichelFactor = 1.174;
inline constexpr double MichelOffset = 22400.0;
inline constexpr double MichelExponent = 0.46;
inline constexpr double TurbulentInitialH = 1.35;
inline constexpr double TurbulentSeparationH = 2.4; // cf reaches zero here...
inline constexpr double SeparationRampH = 2.1;      // ...ramping smoothly from here (no latch: see theory.rst)
inline constexpr double MaxShapeFactor = 3.0;       // cap via the H1 floor; separated wake plateau
inline constexpr double LudwiegTillmannA = 0.246;
inline constexpr double LudwiegTillmannB = 0.678; // 10^(-B*H)
inline constexpr double LudwiegTillmannC = 0.268; // Re_theta^(-C)
inline constexpr double HeadEntrainmentA = 0.0306;
inline constexpr double HeadEntrainmentB = 0.6169; // (H1-3)^(-B)

// --- guards ------------------------------------------------------------------
inline constexpr double MinEdgeVelocityFraction = 0.05; // floor on U_e / V_rel in the BL march
inline constexpr double PostStallDragFactor = 2.0;      // cd growth past the envelope, ~flat-plate normal force
inline constexpr double MaxTranspiration = 0.15;        // |v_t|/V_rel clamp (see theory.rst, "stabilization")
inline constexpr double MinSectionRelaxation = 0.02;    // adaptive-damping floor

// The validity envelope, in incidence from the zero-lift line: the
// transpiration-coupled solve is trusted up to BlendStartDeg, blends
// linearly into the saturated analytic polar by BlendEndDeg, and beyond
// that the analytic polar alone continues (deep stall is outside any
// integral method's physics; a smooth hand-off is worth more than a
// heroic extrapolation -- see theory.rst, "validity").
inline constexpr double BlendStartDeg = 9.0;
inline constexpr double BlendEndDeg = 13.0;

// Continuation in alpha: past this incidence the section problem grows
// multiple stable solutions (early-bubble vs late-transition states), and
// a cold start lands on either unpredictably -- so the solve walks up
// from an attached, single-solution incidence in steps, carrying the
// transpiration state, and follows ONE branch smoothly (see theory.rst).
inline constexpr double ContinuationStartDeg = 6.0;
inline constexpr double ContinuationStepDeg = 1.0;
inline constexpr int    ContinuationStepIterations = 30;

namespace Detail {

// Small dense LU with partial pivoting: the section systems are ~40x40 and
// self-contained, so the solver's LAPACK plumbing would be overkill here.
struct SmallLu {
    int N = 0;
    std::vector<double> A;   // row-major, factored in place
    std::vector<int> Pivot;

    void Factor(std::vector<double> matrix, int n) {
        N = n;
        A = std::move(matrix);
        Pivot.resize(static_cast<std::size_t>(n));
        for (int col = 0; col < n; ++col) {
            int best = col;
            for (int row = col + 1; row < n; ++row)
                if (std::fabs(A[static_cast<std::size_t>(row) * n + col]) >
                    std::fabs(A[static_cast<std::size_t>(best) * n + col]))
                    best = row;
            Pivot[static_cast<std::size_t>(col)] = best;
            if (best != col)
                for (int k = 0; k < n; ++k)
                    std::swap(A[static_cast<std::size_t>(col) * n + k],
                              A[static_cast<std::size_t>(best) * n + k]);
            const double diag = A[static_cast<std::size_t>(col) * n + col];
            if (std::fabs(diag) < 1e-30) continue; // singular column: leave zero row, solve degrades gracefully
            for (int row = col + 1; row < n; ++row) {
                const double factor = A[static_cast<std::size_t>(row) * n + col] / diag;
                A[static_cast<std::size_t>(row) * n + col] = factor;
                for (int k = col + 1; k < n; ++k)
                    A[static_cast<std::size_t>(row) * n + k] -=
                        factor * A[static_cast<std::size_t>(col) * n + k];
            }
        }
    }

    [[nodiscard]] std::vector<double> Solve(std::vector<double> rhs) const {
        for (int col = 0; col < N; ++col) {
            std::swap(rhs[static_cast<std::size_t>(col)], rhs[static_cast<std::size_t>(Pivot[static_cast<std::size_t>(col)])]);
            for (int row = col + 1; row < N; ++row)
                rhs[static_cast<std::size_t>(row)] -= A[static_cast<std::size_t>(row) * N + col] * rhs[static_cast<std::size_t>(col)];
        }
        for (int row = N - 1; row >= 0; --row) {
            for (int k = row + 1; k < N; ++k)
                rhs[static_cast<std::size_t>(row)] -= A[static_cast<std::size_t>(row) * N + k] * rhs[static_cast<std::size_t>(k)];
            const double diag = A[static_cast<std::size_t>(row) * N + row];
            rhs[static_cast<std::size_t>(row)] = (std::fabs(diag) > 1e-30) ? rhs[static_cast<std::size_t>(row)] / diag : 0.0;
        }
        return rhs;
    }
};

// Thwaites shape-factor and shear correlations, Cebeci-Bradshaw piecewise
// fits (see theory.rst).
[[nodiscard]] inline double ThwaitesH(double lambda) {
    const double l = std::clamp(lambda, -0.1, 0.1);
    if (l >= 0.0) return 2.61 - 3.75 * l + 5.24 * l * l;
    return 2.088 + 0.0731 / (l + 0.14);
}

[[nodiscard]] inline double ThwaitesShear(double lambda) {
    const double l = std::clamp(lambda, -0.1, 0.1);
    if (l >= 0.0) return 0.22 + 1.57 * l - 1.8 * l * l;
    return 0.22 + 1.402 * l + 0.018 * l / (l + 0.107);
}

// Head's H1(H) and its inverse direction handled by marching H1 and
// recovering H (see theory.rst).
[[nodiscard]] inline double HeadH1(double H) {
    if (H <= 1.6) return 3.3 + 0.8234 * std::pow(std::max(H - 1.1, 1e-3), -1.287);
    return 3.3 + 1.5501 * std::pow(H - 0.6778, -3.064);
}

[[nodiscard]] inline double HeadHFromH1(double H1) {
    const double h1 = std::max(H1, 3.32);
    if (h1 >= HeadH1(1.6)) // low-H branch
        return 1.1 + std::pow(0.8234 / (h1 - 3.3), 1.0 / 1.287);
    return 0.6778 + std::pow(1.5501 / (h1 - 3.3), 1.0 / 3.064);
}

} // namespace Detail

// -------------------------------------------------------------- the model
/**
 * Transpiration-coupled boundary-layer section solver, implementing the
 * SectionModel interface. CamberSlope(eta, psi) is d(camber)/d(psi) of the
 * mean line at span/radius fraction eta (e.g. Geometry::CamberSlopeAt
 * bound over a contract's blade sections); an empty function means flat
 * sections. Full method in doc/theory.rst.
 */
struct BoundaryLayerSectionModel {
    std::function<double(double eta, double psi)> CamberSlope; ///< d(z/c)/d(psi); empty = flat.
    AnalyticSectionModel Fallback; ///< The saturated polar the envelope blends into.
    int Stations = DefaultSectionStations;
    int MaxIterations = DefaultSectionMaxIterations;
    double Tolerance = DefaultSectionTolerance;
    double Relaxation = DefaultSectionRelaxation;

    [[nodiscard]] SectionCoefficients operator()(const StripSection& strip, double alphaEffDeg,
                                                 double Re, double Ma) const {
        // Blend weight from the incidence relative to the zero-lift line
        // (camber-fair): 1 inside the trusted envelope, 0 past it.
        const double fromZeroLift = std::fabs(alphaEffDeg - strip.Alpha0Deg);
        const double weight =
            std::clamp((BlendEndDeg - fromZeroLift) / (BlendEndDeg - BlendStartDeg), 0.0, 1.0);

        SectionCoefficients analytic = Fallback(strip, alphaEffDeg, Re, Ma);
        if (fromZeroLift > BlendEndDeg) {
            // Deep stall: analytic continuation plus a bluff post-stall
            // drag rise the attached-flow polar has no business knowing.
            const double excess = Math::DegToRad(fromZeroLift - BlendEndDeg);
            analytic.cd += PostStallDragFactor * std::sin(excess) * std::sin(excess);
        }
        if (weight <= 0.0) return analytic;

        const SectionCoefficients viscous = SolveSection(strip, alphaEffDeg, std::max(Re, 1e3));
        return {weight * viscous.cl + (1.0 - weight) * analytic.cl,
                weight * viscous.cd + (1.0 - weight) * analytic.cd};
    }

private:
    [[nodiscard]] SectionCoefficients SolveSection(const StripSection& strip, double alphaDeg,
                                                   double Re) const {
        const int m = std::max(Stations, 8);
        const std::size_t mm = static_cast<std::size_t>(m);

        // --- camber-line geometry (unit chord) ------------------------------
        // Cosine-spaced nodes; z from Simpson integration of the slope.
        std::vector<double> node(mm + 1), zNode(mm + 1);
        for (std::size_t k = 0; k <= mm; ++k)
            node[k] = Math::Half * (1.0 - std::cos(std::numbers::pi * static_cast<double>(k) / m));
        const auto slopeAt = [&](double psi) {
            return CamberSlope ? CamberSlope(strip.Eta, psi) : 0.0;
        };
        zNode[0] = 0.0;
        for (std::size_t k = 0; k < mm; ++k) {
            const double a = node[k], b = node[k + 1];
            const double mid = Math::Half * (a + b);
            zNode[k + 1] = zNode[k] + (b - a) / 6.0 * (slopeAt(a) + 4.0 * slopeAt(mid) + slopeAt(b));
        }
        const auto zAt = [&](std::size_t panel, double psi) {
            // Linear within a panel is enough at this resolution.
            const double t = (psi - node[panel]) / std::max(node[panel + 1] - node[panel], 1e-12);
            return zNode[panel] + t * (zNode[panel + 1] - zNode[panel]);
        };

        struct Station {
            double xv, zv;   // point vortex (quarter point)
            double xc, zc;   // control point (3/4 point)
            double nx, nz;   // unit normal at the control point (+ suction side)
            double tx, tz;   // unit tangent, LE -> TE
            double ds;       // panel arc length
            double sMid;     // arc position of the panel midpoint
        };
        std::vector<Station> st(mm);
        double sAccum = 0.0;
        for (std::size_t k = 0; k < mm; ++k) {
            const double d = node[k + 1] - node[k];
            Station& s = st[k];
            s.xv = node[k] + Math::QuarterChord * d;
            s.zv = zAt(k, s.xv);
            s.xc = node[k] + Math::ThreeQuarterChord * d;
            s.zc = zAt(k, s.xc);
            const double slope = slopeAt(s.xc);
            const double norm = std::sqrt(1.0 + slope * slope);
            s.tx = 1.0 / norm;
            s.tz = slope / norm;
            s.nx = -slope / norm;
            s.nz = 1.0 / norm;
            s.ds = d * norm;
            s.sMid = sAccum + Math::Half * s.ds;
            sAccum += s.ds;
        }

        // --- lumped-vortex influence matrix, factored once ------------------
        // Unit clockwise vortex at (xv,zv): u = dz/(2 pi r^2), w = -dx/(2 pi r^2).
        std::vector<double> aic(mm * mm);
        for (std::size_t i = 0; i < mm; ++i) {
            for (std::size_t j = 0; j < mm; ++j) {
                const double dx = st[i].xc - st[j].xv;
                const double dz = st[i].zc - st[j].zv;
                const double r2 = std::max(dx * dx + dz * dz, 1e-12);
                const double u = dz / (2.0 * std::numbers::pi * r2);
                const double w = -dx / (2.0 * std::numbers::pi * r2);
                aic[i * mm + j] = u * st[i].nx + w * st[i].nz;
            }
        }
        Detail::SmallLu lu;
        lu.Factor(std::move(aic), m);

        // --- transpiration fixed point, with continuation in alpha ----------
        // Adaptively damped: at low Reynolds number the displacement
        // thickness -- and with it the feedback gain -- grows, so a fixed
        // omega that converges at Re ~ 5e5 can cycle at 1e5. Direct
        // viscous-inviscid iteration is furthermore non-convergent once
        // separation is present (the classical Goldstein-singularity
        // behavior that semi-inverse and simultaneous couplings exist to
        // cure), and the low-Re section problem grows MULTIPLE stable
        // states past moderate incidence -- so the solve (a) walks alpha up
        // from the attached regime carrying the transpiration state, to
        // follow one solution branch smoothly, and (b) when the iterate
        // still cycles at the target, returns the ARITHMETIC MEAN over the
        // final iteration window: a phase average of the cycle's attractor,
        // smooth in alpha where the exit iterate is noise. The outer
        // coupling's fixed point needs cl(alpha) smooth far more than it
        // needs any single interior iterate. Full discussion in
        // theory.rst, "stabilization".
        std::vector<double> transp(mm, 0.0); // v_t at control points, per unit V_rel
        std::vector<double> gamma(mm, 0.0);

        const auto solveAtAlpha = [&](double stepAlphaDeg, int budget) -> SectionCoefficients {
            const double a = Math::DegToRad(stepAlphaDeg);
            const double ux = std::cos(a), uz = std::sin(a); // unit freestream
            double cl = 0.0, cd = 0.0, clPrev = 1e30, deltaPrev = 1e30;
            double omega = Relaxation;
            double clWindowSum = 0.0, cdWindowSum = 0.0;
            int windowCount = 0;
            int convergedStreak = 0;

            for (int iter = 0; iter < budget; ++iter) {
                // Inviscid pass under the current transpiration.
                std::vector<double> rhs(mm);
                for (std::size_t i = 0; i < mm; ++i)
                    rhs[i] = transp[i] - (ux * st[i].nx + uz * st[i].nz);
                gamma = lu.Solve(rhs);

                // Edge velocities at the control points: tangential mean plus
                // half the local sheet jump on each side. The own vortex's
                // tangential contribution at its own control point is ~zero
                // for mild camber (the offset is chordwise), so the full sum
                // serves.
                std::vector<double> ueUp(mm), ueLo(mm), ut(mm);
                for (std::size_t i = 0; i < mm; ++i) {
                    double u = ux, w = uz;
                    for (std::size_t j = 0; j < mm; ++j) {
                        const double dx = st[i].xc - st[j].xv;
                        const double dz = st[i].zc - st[j].zv;
                        const double r2 = std::max(dx * dx + dz * dz, 1e-12);
                        u += gamma[j] * dz / (2.0 * std::numbers::pi * r2);
                        w += gamma[j] * (-dx) / (2.0 * std::numbers::pi * r2);
                    }
                    ut[i] = u * st[i].tx + w * st[i].tz;
                    const double sheet = gamma[i] / st[i].ds;
                    ueUp[i] = ut[i] + Math::Half * sheet;
                    ueLo[i] = ut[i] - Math::Half * sheet;
                }

                // Integral boundary layers, LE -> TE on each surface.
                const BoundaryLayerState upper = MarchBoundaryLayer(st, ueUp, Re);
                const BoundaryLayerState lower = MarchBoundaryLayer(st, ueLo, Re);

                // New transpiration: the antisymmetric displacement-flux
                // derivative, central differences over the arc coordinate.
                std::vector<double> mass(mm);
                for (std::size_t i = 0; i < mm; ++i)
                    mass[i] =
                        Math::Half * (upper.ue[i] * upper.deltaStar[i] - lower.ue[i] * lower.deltaStar[i]);
                for (std::size_t i = 0; i < mm; ++i) {
                    const std::size_t lo = (i == 0) ? 0 : i - 1;
                    const std::size_t hi = (i + 1 < mm) ? i + 1 : i;
                    const double dsSpan = std::max(st[hi].sMid - st[lo].sMid, 1e-9);
                    const double fresh =
                        std::clamp((mass[hi] - mass[lo]) / dsSpan, -MaxTranspiration, MaxTranspiration);
                    transp[i] = (1.0 - omega) * transp[i] + omega * fresh;
                }

                cl = 0.0;
                for (std::size_t i = 0; i < mm; ++i) cl += 2.0 * gamma[i];
                cd = upper.SquireYoungDrag() + lower.SquireYoungDrag();

                // Plain mean over the final window (phase average of the cycle).
                if (iter >= budget - SectionOutputAverageWindow) {
                    clWindowSum += cl;
                    cdWindowSum += cd;
                    ++windowCount;
                }

                const double delta = std::fabs(cl - clPrev);
                convergedStreak = (delta < Tolerance) ? convergedStreak + 1 : 0;
                if (convergedStreak >= SectionConvergenceStreak) return {cl, cd}; // genuinely converged
                // Adapt the damping only BEFORE the averaging window: the
                // phase average needs a stationary cycle, and a
                // still-shrinking omega keeps deforming it.
                if (iter < budget - SectionOutputAverageWindow && delta > deltaPrev)
                    omega = std::max(Math::Half * omega, MinSectionRelaxation);
                deltaPrev = delta;
                clPrev = cl;
            }

            return {clWindowSum / std::max(windowCount, 1), cdWindowSum / std::max(windowCount, 1)};
        };

        // Continuation: warm-start the transpiration along an alpha ramp from
        // the attached regime, then solve the target with the full budget.
        const double magnitude = std::fabs(alphaDeg);
        const double direction = (alphaDeg >= 0.0) ? 1.0 : -1.0;
        for (double a = ContinuationStartDeg; a < magnitude; a += ContinuationStepDeg)
            (void)solveAtAlpha(direction * a, ContinuationStepIterations);
        return solveAtAlpha(alphaDeg, MaxIterations);
    }

    struct BoundaryLayerState {
        std::vector<double> ue;        // floored edge velocity used by the march
        std::vector<double> deltaStar; // H * theta
        double thetaTe = 0.0;
        double hTe = 1.4;
        double ueTe = 1.0;

        [[nodiscard]] double SquireYoungDrag() const {
            // cd = 2 theta_TE * Ue_TE^((H_TE+5)/2), unit chord / unit V_rel.
            return 2.0 * thetaTe * std::pow(std::max(ueTe, 1e-3), Math::Half * (hTe + 5.0));
        }
    };

    template <typename Stations>
    [[nodiscard]] BoundaryLayerState MarchBoundaryLayer(const Stations& st,
                                                        const std::vector<double>& ueRaw,
                                                        double Re) const {
        const std::size_t mm = ueRaw.size();
        BoundaryLayerState bl;
        bl.ue.resize(mm);
        bl.deltaStar.assign(mm, 0.0);
        for (std::size_t i = 0; i < mm; ++i) bl.ue[i] = std::max(ueRaw[i], MinEdgeVelocityFraction);

        bool turbulent = false;
        double theta = 0.0, H = 2.61, H1 = Detail::HeadH1(TurbulentInitialH);
        double thwaitesIntegral = 0.0;
        double sPrev = 0.0;
        double thetaLamPrev = 0.0;
        double michelExcessPrev = -1e30;
        double lambdaPrev = 0.0;

        // One Euler step of Head's method over dsLen at edge state
        // (ue, dueds); mutates theta/H1/H. Separation is handled WITHOUT a
        // latch: the skin friction ramps smoothly to zero over
        // [SeparationRampH, TurbulentSeparationH] and the shape factor is
        // capped through the H1 floor, so the marched state is a
        // memoryless, continuous function of the edge-velocity
        // distribution. A latched "separated" flag would make it
        // hysteretic, and hysteresis inside the section is what feeds
        // limit cycles in BOTH coupling loops (see theory.rst,
        // "stabilization").
        const auto advanceTurbulent = [&](double dsLen, double ue, double dueds) {
            if (dsLen <= 0.0) return;
            const double reTheta = std::max(Re * ue * theta, 20.0);
            const double ramp = std::clamp((TurbulentSeparationH - H) /
                                               (TurbulentSeparationH - SeparationRampH),
                                           0.0, 1.0);
            const double cf = ramp * LudwiegTillmannA * std::pow(10.0, -LudwiegTillmannB * H) *
                              std::pow(reTheta, -LudwiegTillmannC);
            const double dThetaDs = Math::Half * cf - (H + 2.0) * theta / ue * dueds;
            const double entrainment =
                HeadEntrainmentA * std::pow(std::max(H1 - 3.0, 1e-3), -HeadEntrainmentB);
            const double dH1Ds = entrainment / std::max(theta, 1e-9) -
                                 H1 * (dThetaDs / std::max(theta, 1e-9) + dueds / ue);
            theta = std::max(theta + dThetaDs * dsLen, 1e-9);
            H1 = std::max(H1 + dH1Ds * dsLen, 3.32);
            H = std::min(Detail::HeadHFromH1(H1), MaxShapeFactor);
        };

        for (std::size_t i = 0; i < mm; ++i) {
            const double s = st[i].sMid;
            const double dsStep = std::max(s - sPrev, 1e-9);
            const double ue = bl.ue[i];
            const double uePrev = (i == 0) ? ue : bl.ue[i - 1];
            const double dueds = (ue - uePrev) / dsStep;

            if (!turbulent) {
                // Thwaites: theta^2 = 0.45/(Re Ue^6) * integral of Ue^5 ds.
                thwaitesIntegral += Math::Half * (std::pow(uePrev, 5) + std::pow(ue, 5)) * dsStep;
                const double theta2 = ThwaitesConstant / (Re * std::pow(ue, 6)) * thwaitesIntegral;
                theta = std::sqrt(std::max(theta2, 1e-16));
                const double lambda = Re * theta2 * dueds;
                H = Detail::ThwaitesH(lambda);

                const double reTheta = Re * ue * theta;
                const double reS = Re * ue * s;
                // Both transition triggers -- Michel's criterion and laminar
                // separation (bubble treated as transition) -- are located
                // WITHIN the step by interpolating their signed indicators:
                // a transition point quantized to whole stations makes
                // cl(alpha) discontinuous, and the coupling's fixed points
                // cannot iterate below that jump (see theory.rst,
                // "stabilization"). Whichever trigger crosses earlier in the
                // step wins.
                const double excess = (reS > 1.0)
                                          ? reTheta - MichelFactor * (1.0 + MichelOffset / reS) *
                                                          std::pow(reS, MichelExponent)
                                          : -1e30;
                bool trigger = false;
                double frac = 1.0;
                if (excess >= 0.0) {
                    trigger = true;
                    frac = (michelExcessPrev < 0.0 && excess > michelExcessPrev)
                               ? -michelExcessPrev / (excess - michelExcessPrev)
                               : 0.0;
                }
                if (lambda < LaminarSeparationLambda) {
                    const double fracSep = (lambdaPrev > LaminarSeparationLambda && lambda < lambdaPrev)
                                               ? (lambdaPrev - LaminarSeparationLambda) / (lambdaPrev - lambda)
                                               : 0.0;
                    frac = trigger ? std::min(frac, fracSep) : fracSep;
                    trigger = true;
                }

                if (trigger) {
                    // Laminar theta at the crossing, then turbulent growth
                    // over the remainder of this step.
                    theta = std::max(thetaLamPrev + frac * (theta - thetaLamPrev), 1e-9);
                    turbulent = true;
                    H = TurbulentInitialH;
                    H1 = Detail::HeadH1(H);
                    advanceTurbulent((1.0 - frac) * dsStep, ue, dueds);
                }
                michelExcessPrev = excess;
                lambdaPrev = lambda;
                thetaLamPrev = theta;
            } else {
                advanceTurbulent(dsStep, ue, dueds);
            }

            bl.deltaStar[i] = H * theta;
            bl.thetaTe = theta;
            bl.hTe = H;
            bl.ueTe = ue;
            sPrev = s;
        }
        return bl;
    }
};

} // namespace Aeolion::Solver
