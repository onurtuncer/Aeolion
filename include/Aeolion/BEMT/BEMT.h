// BEMT/BEMT.h
//
// Blade Element Momentum Theory (BEMT) for a propeller, plus a slipstream
// velocity field meant to be fed into Solver.h's Solve() externalField hook
// so downstream control vanes see the actual propwash (axial acceleration
// + swirl) instead of a uniform freestream.
//
// HOVER-SAFE FORMULATION: classical BEMT is usually written in terms of
// induction FACTORS a = vi/Vinf and a' = wi/(Omega*r). That blows up as
// Vinf -> 0 -- exactly the hover condition that matters most for a
// tail-sitter. This implementation instead solves directly for the
// induced velocities vi (axial) and wi (swirl) themselves, so Vinf=0
// (pure hover) is just an ordinary, well-posed case.
//
// Method, per radial station:
//   Uax = Vinf + vi                  (axial velocity through the disk)
//   Ut  = Omega*r - wi               (tangential velocity through the disk)
//   phi = atan2(Uax, Ut)             (inflow angle)
//   alpha = twist(r) - phi           (blade element angle of attack)
//   cl, cd from the airfoil polar (analytic default, or a lookup table --
//     see Polar -- fill TableAlphaDeg/TableCl/TableCd with real section data)
//   Blade-element thrust/torque (per unit radius, summed over N blades)
//     dT/dr = 0.5*rho*Urel^2*N*chord*Cn
//     dQ/dr = 0.5*rho*Urel^2*N*chord*Ct*r
//   set equal to momentum-theory (annular streamtube) thrust/torque,
//   written directly in terms of vi, wi (this is what keeps it hover-safe):
//     dT/dr = 4*pi*r*rho*Uax*vi*F
//     dQ/dr = 4*pi*r^3*rho*Uax*wi*F
//   iterate (relaxed fixed point) to convergence, with Prandtl tip+hub
//   loss factor F.
//
// This is a mid-fidelity, well-established method (same family as used in
// most propeller/rotor design codes) -- not a substitute for measured
// prop data, but a solid physically-grounded estimate, especially useful
// here because it also gives you the RADIAL DISTRIBUTION of slipstream
// velocity that a downstream vane actually sees, not just total thrust.
//
// No external dependencies beyond the standard library and Aeolion/Math
// (Vec3/Cross, used by the slipstream field function).

#pragma once
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <numbers>
#include "Aeolion/BEMT/PropGeometry.h"
#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Math/Constants.h"

namespace Aeolion::BEMT {

using Solver::Vec3;
using Solver::Cross;
using Solver::Dot;
using Math::Half;
using Math::Two;
using Math::Tiny;
using Math::DegToRad;
using Math::RadToDeg;
using Math::SecondsPerMinute;
using Math::UnitClampLo;
using Math::UnitClampHi;

// --- polar defaults -------------------------------------------------------
inline constexpr double ThinAirfoilClAlpha = Two * std::numbers::pi; // per rad
inline constexpr double DefaultClMax = 1.2;
inline constexpr double DefaultCd0   = 0.02;
inline constexpr double DefaultKCd   = 0.02;   // cd = cd0 + kCd*cl^2

// --- solver defaults ------------------------------------------------------
inline constexpr double SeaLevelDensity   = 1.225; // kg/m^3, ISA
inline constexpr int    DefaultMaxIter    = 3000;
inline constexpr double DefaultTolerance  = 1e-7;
inline constexpr double DefaultRelaxation = 0.25;

// --- momentum-theory / annular-streamtube constants -----------------------
inline constexpr double MomentumDiskFactor = 4.0; // 4*pi*r*rho*Uax*vi*F (thrust); 4*pi*r^3... (torque)

// --- Prandtl tip/hub loss -------------------------------------------------
inline constexpr double EndpointClampFraction = 1e-4; // clamp radius this fraction of span inboard of hub/tip
inline constexpr double MinSinPhi             = 1e-4; // floor on |sin(phi)| in the loss denominator
inline constexpr double MinLossFactor         = 1e-3; // floor so the momentum solve never divides by ~0

// A station this close to the hub or the tip (as a fraction of blade span)
// is treated as a BOUNDARY station carrying no load, rather than iterated.
//
// This is Prandtl's own statement, not a numerical dodge: the loss factor
// goes to zero at both ends because the tip and hub vortices relieve the
// loading there, so a blade element at the boundary produces nothing. The
// momentum balance is correspondingly degenerate -- it divides by that
// vanishing factor -- so asking the fixed point to satisfy it means asking
// for an infinite induced velocity, which the clamps turn into a limit
// cycle rather than an answer.
//
// This is not a rare case. A geometry contract states its blade from the
// hub cutout to the tip (Geometry::PropulsionSpec), so the FIRST station of
// every contract-derived propeller lands exactly on the hub.
inline constexpr double BoundaryStationFraction = 1e-3;

// --- fixed-point iteration guards & seeds ---------------------------------
inline constexpr double SmallVelocity        = 1e-4; // floor on axial velocity in the momentum denominator
inline constexpr double SmallDenominator     = 1e-9; // floor before dividing by the momentum denominator
inline constexpr double AxialInductionSeed   = 0.05; // seed vi ~ this * |omega|*R + floor
inline constexpr double AxialInductionFloor  = 0.1;
inline constexpr double SwirlInductionSeed   = 0.01; // seed wi ~ this * |omega|*r
inline constexpr double MinInductionFactor   = -0.5; // clamp vi/wi lower bound (times |omega|*length)
inline constexpr double MaxAxialInductionFactor = 3.0; // clamp vi upper bound (times |omega|*R)
inline constexpr double MaxSwirlInductionFactor = 0.5; // clamp |wi| (times |omega|*r)

// --- slipstream field -----------------------------------------------------
inline constexpr double SlipstreamRadialEps  = 1e-6; // on-axis guard
inline constexpr double FarWakeDiametersAuto = Two;  // auto development length = this * radius

// ---------------------------------------------------------- Airfoil polar
// Analytic default (smooth-saturating lift curve + parabolic drag polar)
// or an optional table lookup (e.g. from measured or CFD section data) --
// set UseTable=true and fill TableAlphaDeg/TableCl/TableCd (sorted by
// alpha) to use real section data instead.
/**
 * Section lift/drag polar: analytic default (smooth-saturating lift curve +
 * parabolic drag polar), or a table lookup when UseTable is set (fill
 * TableAlphaDeg/TableCl/TableCd, sorted by alpha, with real section data).
 */
struct Polar {
    double cl_alpha = ThinAirfoilClAlpha; ///< Per rad, thin-airfoil default.
    double alpha0Deg = 0.0;     ///< Zero-lift angle.
    double clMax = DefaultClMax;
    double cd0 = DefaultCd0;
    double kCd = DefaultKCd;

    bool UseTable = false;
    std::vector<double> TableAlphaDeg, TableCl, TableCd;
};

/** Evaluate a Polar at alphaDeg, writing the resulting cl/cd. */
inline void EvalPolar(const Polar& p, double alphaDeg, double& cl, double& cd) {
    if (p.UseTable && p.TableAlphaDeg.size() >= 2) {
        const auto& A = p.TableAlphaDeg;
        double a = std::clamp(alphaDeg, A.front(), A.back());
        size_t i = 0;
        while (i + 1 < A.size() && A[i + 1] < a) ++i;
        double t = (A[i + 1] - A[i] > Tiny) ? (a - A[i]) / (A[i + 1] - A[i]) : 0.0;
        cl = p.TableCl[i] + t * (p.TableCl[i + 1] - p.TableCl[i]);
        cd = p.TableCd[i] + t * (p.TableCd[i + 1] - p.TableCd[i]);
        return;
    }
    double alpha = DegToRad(alphaDeg - p.alpha0Deg);
    // smooth saturation toward +-clMax instead of a hard clip -- keeps the
    // BEMT fixed-point iteration well behaved near stall rather than
    // handing it a kink.
    cl = p.clMax * std::tanh(p.cl_alpha * alpha / p.clMax);
    cd = p.cd0 + p.kCd * cl * cl;
}

// Blade geometry lives in BEMT/PropGeometry.h so that describing a
// propeller does not require including this solver.

// -------------------------------------------------------------- Solving
/** BEMT solution at one radial blade station. */
struct StationResult {
    double r = 0, vi = 0, wi = 0, phiDeg = 0, alphaDeg = 0, cl = 0, cd = 0;
    double dT_dr = 0, dQ_dr = 0;

    /**
     * Per-station convergence. Result::Converged is the AND of these, which
     * on its own cannot say whether one awkward station failed or the whole
     * blade did -- and those call for very different responses.
     */
    bool Converged = true;
    int Iterations = 0;
    double Residual = 0.0; ///< Last |change| in induced velocity [m/s].
};

/** Full BEMT solve result: per-station breakdown plus integrated thrust/torque/power. */
struct Result {
    std::vector<StationResult> Stations;
    double Thrust = 0, Torque = 0, Power = 0;
    double rpm = 0, omega = 0;
    PropGeometry Geom;
    bool Converged = true;
};

/** Prandtl tip/hub loss factor F at radius rIn and inflow angle phi. */
[[nodiscard]] inline double TipHubLoss(const PropGeometry& g, double rIn, double phi) {
    // The Prandtl loss factor is *designed* to go to exactly zero right at
    // the tip and hub -- correct for the physical loading, but a singular
    // denominator for the momentum-theory induced-velocity solve, which
    // divides by F. Stations placed exactly at r=HubRadius or r=Radius
    // (an easy mistake when generating a station distribution) will blow
    // up the iteration. Clamp the radius used here slightly inboard of
    // the true endpoints so F stays small-but-finite instead of exactly
    // zero; the blade geometry itself (chord/twist/actual r) is untouched.
    double span = std::max(g.Radius - g.HubRadius, Tiny);
    double r = std::clamp(rIn, g.HubRadius + EndpointClampFraction * span, g.Radius - EndpointClampFraction * span);

    double sinPhi = std::max(std::fabs(std::sin(phi)), MinSinPhi);
    double fTip = (g.NBlades * Half) * (g.Radius - r) / (r * sinPhi);
    double fHub = (g.NBlades * Half) * (r - g.HubRadius) / (r * sinPhi);
    fTip = std::max(fTip, 0.0); fHub = std::max(fHub, 0.0);
    double Ftip = (Two / std::numbers::pi) * std::acos(std::clamp(std::exp(-fTip), UnitClampLo, UnitClampHi));
    double Fhub = (Two / std::numbers::pi) * std::acos(std::clamp(std::exp(-fHub), UnitClampLo, UnitClampHi));
    return std::clamp(Ftip * Fhub, MinLossFactor, UnitClampHi);
}

/**
 * Solve the hover-safe BEMT fixed point for every radial station of geom,
 * at rotor speed rpm and freestream Vinf, then integrate thrust/torque/power.
 */
[[nodiscard]] inline Result Solve(const PropGeometry& geom, const Polar& polar, double rpm, double Vinf,
                     double rho = SeaLevelDensity, int maxIter = DefaultMaxIter,
                     double tol = DefaultTolerance, double relax = DefaultRelaxation) {
    Result res;
    res.Geom = geom;
    res.rpm = rpm;
    res.omega = rpm * Two * std::numbers::pi / SecondsPerMinute;
    double omega = res.omega;

    const double bladeSpan = std::max(geom.Radius - geom.HubRadius, Tiny);

    for (const auto& st : geom.Stations) {
        double r = st.r;

        // Boundary station: unloaded by the tip/hub relief, so there is
        // nothing to solve. Recording it as converged is honest -- the
        // answer is known, not merely unreached.
        if (r <= geom.HubRadius + BoundaryStationFraction * bladeSpan ||
            r >= geom.Radius - BoundaryStationFraction * bladeSpan) {
            StationResult boundary;
            boundary.r = r;
            boundary.phiDeg = RadToDeg(std::atan2(Vinf, omega * r));
            boundary.alphaDeg = st.TwistDeg - boundary.phiDeg;
            res.Stations.push_back(boundary);
            continue;
        }

        double vi = AxialInductionSeed * std::fabs(omega) * geom.Radius + AxialInductionFloor; // seed, refined by iteration
        double wi = SwirlInductionSeed * std::fabs(omega) * r;
        bool ok = false;
        int used = 0;
        double residual = 0.0;

        for (int iter = 0; iter < maxIter; ++iter) {
            used = iter + 1;
            double Uax = Vinf + vi;
            double Ut = omega * r - wi;
            double phi = std::atan2(Uax, Ut);
            double alphaDeg = st.TwistDeg - RadToDeg(phi);

            double cl, cd;
            EvalPolar(polar, alphaDeg, cl, cd);

            double Cn = cl * std::cos(phi) - cd * std::sin(phi);
            double Ct = cl * std::sin(phi) + cd * std::cos(phi);
            double Urel2 = Uax * Uax + Ut * Ut;

            double F = TipHubLoss(geom, r, phi);

            // Blade-element dT/dr, dQ/dr:
            double dTdr_be = Half * rho * Urel2 * geom.NBlades * st.Chord * Cn;
            double dQdr_be = Half * rho * Urel2 * geom.NBlades * st.Chord * Ct * r;

            // Momentum-theory dT/dr, dQ/dr, solved for vi_new/wi_new by
            // equating to the blade-element values above:
            double denomT = MomentumDiskFactor * std::numbers::pi * r * rho * std::max(Uax, SmallVelocity) * F;
            double viNew = std::clamp(dTdr_be / std::max(denomT, SmallDenominator),
                                      MinInductionFactor * std::fabs(omega) * geom.Radius,
                                      MaxAxialInductionFactor * std::fabs(omega) * geom.Radius);
            double denomQ = MomentumDiskFactor * std::numbers::pi * r * r * rho * std::max(Uax, SmallVelocity) * F;
            double wiNew = std::clamp(dQdr_be / std::max(denomQ, SmallDenominator),
                                      -MaxSwirlInductionFactor * std::fabs(omega) * r,
                                      MaxSwirlInductionFactor * std::fabs(omega) * r);

            double dv = viNew - vi, dw = wiNew - wi;
            residual = std::max(std::fabs(dv), std::fabs(dw));
            vi += relax * dv;
            wi += relax * dw;

            if (std::fabs(dv) < tol * std::max(1.0, std::fabs(vi)) &&
                std::fabs(dw) < tol * std::max(1.0, std::fabs(wi))) {
                ok = true;
                break;
            }
        }
        if (!ok) res.Converged = false;

        double Uax = Vinf + vi, Ut = omega * r - wi;
        double phi = std::atan2(Uax, Ut);
        double alphaDeg = st.TwistDeg - RadToDeg(phi);
        double cl, cd; EvalPolar(polar, alphaDeg, cl, cd);
        double Cn = cl * std::cos(phi) - cd * std::sin(phi);
        double Ct = cl * std::sin(phi) + cd * std::cos(phi);
        double Urel2 = Uax * Uax + Ut * Ut;

        StationResult sr;
        sr.Converged = ok; sr.Iterations = used; sr.Residual = residual;
        sr.r = r; sr.vi = vi; sr.wi = wi; sr.phiDeg = RadToDeg(phi); sr.alphaDeg = alphaDeg;
        sr.cl = cl; sr.cd = cd;
        sr.dT_dr = Half * rho * Urel2 * geom.NBlades * st.Chord * Cn;
        sr.dQ_dr = Half * rho * Urel2 * geom.NBlades * st.Chord * Ct * r;
        res.Stations.push_back(sr);
    }

    // Integrate thrust/torque (trapezoidal over r).
    for (size_t i = 0; i + 1 < res.Stations.size(); ++i) {
        double dr = res.Stations[i + 1].r - res.Stations[i].r;
        res.Thrust += Half * (res.Stations[i].dT_dr + res.Stations[i + 1].dT_dr) * dr;
        res.Torque += Half * (res.Stations[i].dQ_dr + res.Stations[i + 1].dQ_dr) * dr;
    }
    res.Power = res.Torque * std::fabs(omega);
    return res;
}

// ------------------------------------------------------ Performance metrics
// The two hard bounds a propeller result must respect. Both are ratios of
// an ideal power to the power actually absorbed, so both are <= 1 for any
// physically sensible answer -- they are thermodynamic constraints, not
// tuning targets. TestBEMT checks exactly this, and it is how two real
// sign/scale bugs were caught during development.

/**
 * Figure of merit: hover efficiency, ideal momentum-theory induced power
 * over shaft power. Meaningful only at (or near) zero forward speed.
 */
[[nodiscard]] inline double FigureOfMerit(const Result& result, double rho = SeaLevelDensity) {
    if (result.Power <= Tiny || result.Thrust <= 0.0) return 0.0;
    const double diskArea = std::numbers::pi * result.Geom.Radius * result.Geom.Radius;
    const double idealPower = result.Thrust * std::sqrt(result.Thrust / (Two * rho * diskArea));
    return idealPower / result.Power;
}

/**
 * Propulsive efficiency in forward flight: useful power out (thrust times
 * airspeed) over shaft power in. Zero in hover by definition, since a
 * stationary propeller does no useful work however much thrust it makes.
 */
[[nodiscard]] inline double PropulsiveEfficiency(const Result& result, double Vinf) {
    if (result.Power <= Tiny || result.Thrust <= 0.0 || Vinf <= 0.0) return 0.0;
    return result.Thrust * Vinf / result.Power;
}

/**
 * Disk loading [N/m^2] -- thrust per unit disk area, the quantity that sets
 * induced velocity and therefore hover efficiency.
 */
[[nodiscard]] inline double DiskLoading(const Result& result) {
    const double diskArea = std::numbers::pi * result.Geom.Radius * result.Geom.Radius;
    return (diskArea > Tiny) ? result.Thrust / diskArea : 0.0;
}

// --------------------------------------------------------- Slipstream field
// Meant to be passed as Solver::Solve()'s externalField callback: returns the
// INDUCED velocity perturbation (axial + swirl) at a global point P, to be
// added on top of the vane surface's own freestream. Assumes the prop axis
// is aligned with AxisDir (unit vector) through HubCenter, thrust directed
// along +AxisDir, and rotation sense given by Geom.RotationSign about that
// axis (right-hand rule).
//
// DevelopmentLength: distance downstream (along AxisDir) over which the
// axial slipstream velocity develops from its at-disk value toward the
// classical far-wake value (~2x the disk value, from simple momentum
// theory / mass continuity). A common engineering default is about one
// prop diameter; tune shorter if your vanes sit very close behind the
// disk and you want a more "at-disk" estimate. Swirl is modeled as
// roughly constant with downstream distance (no development ramp) --
// viscous mixing will erode it further downstream than this model
// accounts for, so don't trust this far aft of the disk.
/**
 * Callable slipstream velocity field, meant to be passed as
 * Solver::Solve()'s externalField callback: returns the induced velocity
 * perturbation (axial + swirl) at a global point P, to be added on top of
 * the vane surface's own freestream.
 */
struct SlipstreamField {
    Result BEMTResult;
    Vec3 HubCenter{0, 0, 0};
    Vec3 AxisDir{1, 0, 0}; ///< Unit vector, thrust direction.
    double DevelopmentLength = -1.0; ///< <0 => auto (1 prop diameter).

    /** Induced velocity perturbation at global point P. */
    [[nodiscard]] Vec3 operator()(const Vec3& P) const {
        Vec3 rel = P - HubCenter;
        double x = Dot(rel, AxisDir);
        Vec3 radial = rel - AxisDir * x;
        double r = radial.Norm();
        if (x < 0.0 || r < SlipstreamRadialEps || BEMTResult.Stations.empty()) return Vec3(0, 0, 0);

        double R = BEMTResult.Geom.Radius, rHub = BEMTResult.Geom.HubRadius;
        if (r > R) return Vec3(0, 0, 0); // no slipstream contraction modeled -- outside disk tip = no effect

        // interpolate vi(r), wi(r) from BEMT stations
        const auto& st = BEMTResult.Stations;
        double rc = std::clamp(r, rHub, R);
        size_t i = 0;
        while (i + 1 < st.size() && st[i + 1].r < rc) ++i;
        double vi, wi;
        if (i + 1 < st.size()) {
            double t = (st[i + 1].r - st[i].r > Tiny) ? (rc - st[i].r) / (st[i + 1].r - st[i].r) : 0.0;
            vi = st[i].vi + t * (st[i + 1].vi - st[i].vi);
            wi = st[i].wi + t * (st[i + 1].wi - st[i].wi);
        } else {
            vi = st.back().vi; wi = st.back().wi;
        }

        double devLen = (DevelopmentLength > 0.0) ? DevelopmentLength : FarWakeDiametersAuto * R;
        double dev = std::clamp(x / devLen, 0.0, 1.0); // 0 at disk -> 1 by devLen downstream
        double viEff = vi * (1.0 + dev); // ramps disk value -> ~2x far-wake value
        double wiEff = wi;               // swirl modeled as roughly constant downstream (see note above)

        Vec3 rhat = radial * (1.0 / r);
        Vec3 that = Cross(AxisDir, rhat) * (double)BEMTResult.Geom.RotationSign;

        return AxisDir * viEff + that * wiEff;
    }
};

} // namespace Aeolion::BEMT
