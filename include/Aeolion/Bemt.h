// Bemt.h
//
// Blade Element Momentum Theory (BEMT) for a propeller, plus a slipstream
// velocity field meant to be fed into VLM.h's Solve() externalField hook
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
#include "Aeolion/Math/Vec3.h"

namespace Aeolion::Bemt {

using VLM::Vec3;
using VLM::Cross;
using VLM::Dot;

// ---------------------------------------------------------- Airfoil polar
// Analytic default (smooth-saturating lift curve + parabolic drag polar)
// or an optional table lookup (e.g. from measured or CFD section data) --
// set UseTable=true and fill TableAlphaDeg/TableCl/TableCd (sorted by
// alpha) to use real section data instead.
struct Polar {
    double cl_alpha = 2.0 * std::numbers::pi; // per rad, thin-airfoil default
    double alpha0Deg = 0.0;     // zero-lift angle
    double clMax = 1.2;
    double cd0 = 0.02;
    double kCd = 0.02;          // cd = cd0 + kCd*cl^2

    bool UseTable = false;
    std::vector<double> TableAlphaDeg, TableCl, TableCd;
};

inline void EvalPolar(const Polar& p, double alphaDeg, double& cl, double& cd) {
    if (p.UseTable && p.TableAlphaDeg.size() >= 2) {
        const auto& A = p.TableAlphaDeg;
        double a = std::clamp(alphaDeg, A.front(), A.back());
        size_t i = 0;
        while (i + 1 < A.size() && A[i + 1] < a) ++i;
        double t = (A[i + 1] - A[i] > 1e-9) ? (a - A[i]) / (A[i + 1] - A[i]) : 0.0;
        cl = p.TableCl[i] + t * (p.TableCl[i + 1] - p.TableCl[i]);
        cd = p.TableCd[i] + t * (p.TableCd[i + 1] - p.TableCd[i]);
        return;
    }
    double alpha = (alphaDeg - p.alpha0Deg) * std::numbers::pi / 180.0;
    // smooth saturation toward +-clMax instead of a hard clip -- keeps the
    // BEMT fixed-point iteration well behaved near stall rather than
    // handing it a kink.
    cl = p.clMax * std::tanh(p.cl_alpha * alpha / p.clMax);
    cd = p.cd0 + p.kCd * cl * cl;
}

// -------------------------------------------------------- Blade geometry
struct BladeStation {
    double r;         // radial position [m]
    double Chord;     // local chord [m]
    double TwistDeg;  // local geometric pitch angle [deg] (blade chordline
                       // angle relative to the rotor plane)
};

struct PropGeometry {
    int NBlades = 2;
    double Radius = 0.15;     // tip radius [m]
    double HubRadius = 0.02;  // hub radius [m]
    std::vector<BladeStation> Stations; // sorted by increasing r, spans [HubRadius, Radius]
    int RotationSign = 1;     // +1 or -1: sense of blade rotation about the prop axis (right-hand rule about +thrust axis)
};

// -------------------------------------------------------------- Solving
struct StationResult {
    double r = 0, vi = 0, wi = 0, phiDeg = 0, alphaDeg = 0, cl = 0, cd = 0;
    double dT_dr = 0, dQ_dr = 0;
};

struct Result {
    std::vector<StationResult> Stations;
    double Thrust = 0, Torque = 0, Power = 0;
    double rpm = 0, omega = 0;
    PropGeometry Geom;
    bool Converged = true;
};

inline double TipHubLoss(const PropGeometry& g, double rIn, double phi) {
    // The Prandtl loss factor is *designed* to go to exactly zero right at
    // the tip and hub -- correct for the physical loading, but a singular
    // denominator for the momentum-theory induced-velocity solve, which
    // divides by F. Stations placed exactly at r=HubRadius or r=Radius
    // (an easy mistake when generating a station distribution) will blow
    // up the iteration. Clamp the radius used here slightly inboard of
    // the true endpoints so F stays small-but-finite instead of exactly
    // zero; the blade geometry itself (chord/twist/actual r) is untouched.
    double span = std::max(g.Radius - g.HubRadius, 1e-9);
    double r = std::clamp(rIn, g.HubRadius + 1e-4 * span, g.Radius - 1e-4 * span);

    double sinPhi = std::max(std::fabs(std::sin(phi)), 1e-4);
    double fTip = (g.NBlades / 2.0) * (g.Radius - r) / (r * sinPhi);
    double fHub = (g.NBlades / 2.0) * (r - g.HubRadius) / (r * sinPhi);
    fTip = std::max(fTip, 0.0); fHub = std::max(fHub, 0.0);
    double Ftip = (2.0 / std::numbers::pi) * std::acos(std::clamp(std::exp(-fTip), -1.0, 1.0));
    double Fhub = (2.0 / std::numbers::pi) * std::acos(std::clamp(std::exp(-fHub), -1.0, 1.0));
    return std::clamp(Ftip * Fhub, 1e-3, 1.0);
}

inline Result Solve(const PropGeometry& geom, const Polar& polar, double rpm, double Vinf,
                     double rho = 1.225, int maxIter = 3000, double tol = 1e-7, double relax = 0.25) {
    Result res;
    res.Geom = geom;
    res.rpm = rpm;
    res.omega = rpm * 2.0 * std::numbers::pi / 60.0;
    double omega = res.omega;

    for (const auto& st : geom.Stations) {
        double r = st.r;
        double vi = 0.05 * std::fabs(omega) * geom.Radius + 0.1; // seed guesses, refined by iteration
        double wi = 0.01 * std::fabs(omega) * r;
        bool ok = false;

        for (int iter = 0; iter < maxIter; ++iter) {
            double Uax = Vinf + vi;
            double Ut = omega * r - wi;
            double phi = std::atan2(Uax, Ut);
            double alphaDeg = st.TwistDeg - phi * 180.0 / std::numbers::pi;

            double cl, cd;
            EvalPolar(polar, alphaDeg, cl, cd);

            double Cn = cl * std::cos(phi) - cd * std::sin(phi);
            double Ct = cl * std::sin(phi) + cd * std::cos(phi);
            double Urel2 = Uax * Uax + Ut * Ut;

            double F = TipHubLoss(geom, r, phi);

            // Blade-element dT/dr, dQ/dr:
            double dTdr_be = 0.5 * rho * Urel2 * geom.NBlades * st.Chord * Cn;
            double dQdr_be = 0.5 * rho * Urel2 * geom.NBlades * st.Chord * Ct * r;

            // Momentum-theory dT/dr, dQ/dr, solved for vi_new/wi_new by
            // equating to the blade-element values above:
            double denomT = 4.0 * std::numbers::pi * r * rho * std::max(Uax, 1e-4) * F;
            double viNew = std::clamp(dTdr_be / std::max(denomT, 1e-9), -0.5 * std::fabs(omega) * geom.Radius, 3.0 * std::fabs(omega) * geom.Radius);
            double denomQ = 4.0 * std::numbers::pi * r * r * rho * std::max(Uax, 1e-4) * F;
            double wiNew = std::clamp(dQdr_be / std::max(denomQ, 1e-9), -0.5 * std::fabs(omega) * r, 0.5 * std::fabs(omega) * r);

            double dv = viNew - vi, dw = wiNew - wi;
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
        double alphaDeg = st.TwistDeg - phi * 180.0 / std::numbers::pi;
        double cl, cd; EvalPolar(polar, alphaDeg, cl, cd);
        double Cn = cl * std::cos(phi) - cd * std::sin(phi);
        double Ct = cl * std::sin(phi) + cd * std::cos(phi);
        double Urel2 = Uax * Uax + Ut * Ut;

        StationResult sr;
        sr.r = r; sr.vi = vi; sr.wi = wi; sr.phiDeg = phi * 180.0 / std::numbers::pi; sr.alphaDeg = alphaDeg;
        sr.cl = cl; sr.cd = cd;
        sr.dT_dr = 0.5 * rho * Urel2 * geom.NBlades * st.Chord * Cn;
        sr.dQ_dr = 0.5 * rho * Urel2 * geom.NBlades * st.Chord * Ct * r;
        res.Stations.push_back(sr);
    }

    // Integrate thrust/torque (trapezoidal over r).
    for (size_t i = 0; i + 1 < res.Stations.size(); ++i) {
        double dr = res.Stations[i + 1].r - res.Stations[i].r;
        res.Thrust += 0.5 * (res.Stations[i].dT_dr + res.Stations[i + 1].dT_dr) * dr;
        res.Torque += 0.5 * (res.Stations[i].dQ_dr + res.Stations[i + 1].dQ_dr) * dr;
    }
    res.Power = res.Torque * std::fabs(omega);
    return res;
}

// --------------------------------------------------------- Slipstream field
// Meant to be passed as VLM::Solve()'s externalField callback: returns the
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
struct SlipstreamField {
    Result BemtResult;
    Vec3 HubCenter{0, 0, 0};
    Vec3 AxisDir{1, 0, 0}; // unit vector, thrust direction
    double DevelopmentLength = -1.0; // <0 => auto (1 prop diameter)

    Vec3 operator()(const Vec3& P) const {
        Vec3 rel = P - HubCenter;
        double x = Dot(rel, AxisDir);
        Vec3 radial = rel - AxisDir * x;
        double r = radial.Norm();
        if (x < 0.0 || r < 1e-6 || BemtResult.Stations.empty()) return Vec3(0, 0, 0);

        double R = BemtResult.Geom.Radius, rHub = BemtResult.Geom.HubRadius;
        if (r > R) return Vec3(0, 0, 0); // no slipstream contraction modeled -- outside disk tip = no effect

        // interpolate vi(r), wi(r) from BEMT stations
        const auto& st = BemtResult.Stations;
        double rc = std::clamp(r, rHub, R);
        size_t i = 0;
        while (i + 1 < st.size() && st[i + 1].r < rc) ++i;
        double vi, wi;
        if (i + 1 < st.size()) {
            double t = (st[i + 1].r - st[i].r > 1e-9) ? (rc - st[i].r) / (st[i + 1].r - st[i].r) : 0.0;
            vi = st[i].vi + t * (st[i + 1].vi - st[i].vi);
            wi = st[i].wi + t * (st[i + 1].wi - st[i].wi);
        } else {
            vi = st.back().vi; wi = st.back().wi;
        }

        double devLen = (DevelopmentLength > 0.0) ? DevelopmentLength : 2.0 * R;
        double dev = std::clamp(x / devLen, 0.0, 1.0); // 0 at disk -> 1 by devLen downstream
        double viEff = vi * (1.0 + dev); // ramps disk value -> ~2x far-wake value
        double wiEff = wi;               // swirl modeled as roughly constant downstream (see note above)

        Vec3 rhat = radial * (1.0 / r);
        Vec3 that = Cross(AxisDir, rhat) * (double)BemtResult.Geom.RotationSign;

        return AxisDir * viEff + that * wiEff;
    }
};

} // namespace Aeolion::Bemt
