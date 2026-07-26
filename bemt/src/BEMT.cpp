// BEMT/BEMT.cpp
//
// Implementation of the hover-safe BEMT fixed point and the slipstream
// field it feeds downstream vanes. See BEMT.h for the method and the
// rationale behind solving for induced velocities directly rather than
// induction factors.

#include "Aeolion/BEMT/BEMT.h"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace Aeolion::BEMT {

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

namespace {

// --- momentum-theory / annular-streamtube constants -----------------------
constexpr double MomentumDiskFactor = 4.0; // 4*pi*r*rho*Uax*vi*F (thrust); 4*pi*r^3... (torque)

// --- Prandtl tip/hub loss -------------------------------------------------
constexpr double EndpointClampFraction = 1e-4; // clamp radius this fraction of span inboard of hub/tip
constexpr double MinSinPhi             = 1e-4; // floor on |sin(phi)| in the loss denominator
constexpr double MinLossFactor         = 1e-3; // floor so the momentum solve never divides by ~0

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
constexpr double BoundaryStationFraction = 1e-3;

// --- fixed-point iteration guards & seeds ---------------------------------
constexpr double SmallVelocity        = 1e-4; // floor on axial velocity in the momentum denominator
constexpr double SmallDenominator     = 1e-9; // floor before dividing by the momentum denominator
constexpr double AxialInductionSeed   = 0.05; // seed vi ~ this * |omega|*R + floor
constexpr double AxialInductionFloor  = 0.1;
constexpr double SwirlInductionSeed   = 0.01; // seed wi ~ this * |omega|*r
constexpr double MinInductionFactor   = -0.5; // clamp vi/wi lower bound (times |omega|*length)
constexpr double MaxAxialInductionFactor = 3.0; // clamp vi upper bound (times |omega|*R)
constexpr double MaxSwirlInductionFactor = 0.5; // clamp |wi| (times |omega|*r)

// --- slipstream field -----------------------------------------------------
constexpr double SlipstreamRadialEps  = 1e-6; // on-axis guard

} // namespace

void EvalPolar(const Polar& p, double alphaDeg, double& cl, double& cd) {
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

double TipHubLoss(const PropGeometry& g, double rIn, double phi) {
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

Result Solve(const PropGeometry& geom, const Polar& polar, double rpm, double Vinf, double rho, int maxIter,
            double tol, double relax) {
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

double FigureOfMerit(const Result& result, double rho) {
    if (result.Power <= Tiny || result.Thrust <= 0.0) return 0.0;
    const double diskArea = std::numbers::pi * result.Geom.Radius * result.Geom.Radius;
    const double idealPower = result.Thrust * std::sqrt(result.Thrust / (Two * rho * diskArea));
    return idealPower / result.Power;
}

double PropulsiveEfficiency(const Result& result, double Vinf) {
    if (result.Power <= Tiny || result.Thrust <= 0.0 || Vinf <= 0.0) return 0.0;
    return result.Thrust * Vinf / result.Power;
}

double DiskLoading(const Result& result) {
    const double diskArea = std::numbers::pi * result.Geom.Radius * result.Geom.Radius;
    return (diskArea > Tiny) ? result.Thrust / diskArea : 0.0;
}

Vec3 SlipstreamField::operator()(const Vec3& P) const {
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

} // namespace Aeolion::BEMT
