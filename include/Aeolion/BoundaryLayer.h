// BoundaryLayer.h
//
// 2D integral boundary-layer method, marched along a known inviscid
// surface velocity distribution Ue(s) (e.g. from AirfoilPanel.h) from
// the stagnation point to the trailing edge, on each surface (upper/
// lower) separately. This is the same family of method XFOIL uses for its
// boundary-layer module (minus the coupled viscous-inviscid iteration --
// this is a one-shot analysis on a fixed inviscid Ue(s), which is the
// standard "uncoupled" simplification and is adequate as long as the
// boundary layer stays thin and attached).
//
// Method:
//   1. LAMINAR region: Thwaites' method (closed-form momentum-thickness
//      integral + the standard two-branch H(lambda)/l(lambda) closure).
//   2. TRANSITION: Michel's criterion (Re_theta vs Re_x correlation).
//   3. TURBULENT region: Head's entrainment method (von Karman momentum
//      integral + entrainment equation, Ludwieg-Tillmann skin friction,
//      Head's H1(H) and CE(H1) correlations).
//   4. Section profile drag: Squire-Young formula, applied at the
//      trailing edge of each surface and summed.
//
// This is a mid-fidelity method (same accuracy class as XFOIL's viscous
// module on an uncoupled/fixed inviscid solution): good for attached,
// mildly-separated-at-most flow. It will flag (not correctly compute)
// significant separation -- see SurfaceBLResult::Separated.
//
// No external dependencies beyond the standard library.

#pragma once
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace Aeolion { namespace BoundaryLayer {

// ---- Thwaites closure (standard two-branch fit, Thwaites/Curle) --------
inline double ThwaitesL(double lambda) {
    if (lambda >= 0.0) return 0.22 + 1.57 * lambda - 1.8 * lambda * lambda;
    return 0.22 + 1.402 * lambda + 0.018 * lambda / (lambda + 0.107);
}
inline double ThwaitesH(double lambda) {
    if (lambda >= 0.0) return 2.61 - 3.75 * lambda + 5.24 * lambda * lambda;
    return 2.088 + 0.0731 / (lambda + 0.14);
}
constexpr double ThwaitesSeparationLambda = -0.09;

// ---- Michel transition criterion ----------------------------------------
inline double MichelReThetaCrit(double Rex) {
    return 1.174 * (1.0 + 22400.0 / std::max(Rex, 1.0)) * std::pow(Rex, 0.46);
}

// ---- Head's method closures ---------------------------------------------
inline double HeadH1OfH(double H) {
    if (H <= 1.6) return 3.3 + 0.8234 * std::pow(std::max(H - 1.1, 1e-6), -1.287);
    return 3.732 + 0.04 * std::pow(std::max(H - 1.5, 1e-6), -2.5);
}
// Inverse of HeadH1OfH via bisection (H1 is monotonically decreasing in H
// over the physically relevant range).
inline double HeadHOfH1(double H1) {
    double lo = 1.05, hi = 6.0;
    // ensure bracket: HeadH1OfH is decreasing, so H1(lo) should be > target > H1(hi)
    for (int iter = 0; iter < 60; ++iter) {
        double mid = 0.5 * (lo + hi);
        double h1mid = HeadH1OfH(mid);
        if (h1mid > H1) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}
inline double HeadCE(double H1) {
    return 0.0306 * std::pow(std::max(H1 - 3.0, 1e-6), -0.6169);
}
inline double LudwiegTillmannCf(double Re_theta, double H) {
    return 0.246 * std::pow(10.0, -0.678 * H) * std::pow(std::max(Re_theta, 1.0), -0.268);
}
constexpr double TurbulentSeparationH = 2.4; // approximate onset-of-separation flag threshold

// ---- one surface (stagnation point -> trailing edge) --------------------
struct SurfaceBLResult {
    double theta_te = 0.0;      // momentum thickness at TE
    double H_te = 0.0;          // shape factor at TE
    double Ue_te = 0.0;         // edge velocity at TE (same units as input Ue)
    double cd = 0.0;            // Squire-Young section drag contribution from this surface
    double sTransition = -1.0;  // arc length at transition (-1 if never transitioned / stayed laminar)
    bool Separated = false;     // true if a separation criterion was crossed anywhere along the run
    double sSeparation = -1.0;
};

// s: cumulative arc length from stagnation point (s[0]=0, increasing toward TE)
// Ue: edge speed magnitude at each s (same length as s), assumed >0 (flow moving away from stagnation)
// Vinf: freestream reference speed, nu: kinematic viscosity [m^2/s], chord: for Squire-Young normalization
inline SurfaceBLResult AnalyzeSurface(const std::vector<double>& s, const std::vector<double>& Ue,
                                       double Vinf, double nu, double chord,
                                       bool forceTurbulentFromStart = false) {
    if (s.size() < 3 || s.size() != Ue.size()) throw std::invalid_argument("AnalyzeSurface: bad input arrays");
    size_t n = s.size();

    SurfaceBLResult res;

    // ---- Laminar region: Thwaites (closed-form integral) ----
    // theta^2(x) = 0.45*nu/Ue(x)^6 * integral_0^x Ue(xi)^5 dxi
    std::vector<double> intUe5(n, 0.0);
    for (size_t i = 1; i < n; ++i) {
        double ds = s[i] - s[i - 1];
        double avg = 0.5 * (std::pow(Ue[i - 1], 5) + std::pow(Ue[i], 5));
        intUe5[i] = intUe5[i - 1] + avg * ds;
    }

    size_t transIdx = n - 1; // default: never transitions (stays laminar to TE)
    double theta = 0.0, H = 2.2, lambda = 0.0;
    bool transitioned = forceTurbulentFromStart;

    if (!forceTurbulentFromStart) {
        double sEnd = s.back();
        for (size_t i = 1; i < n; ++i) {
            double Uei = std::max(Ue[i], 1e-9);
            double theta2 = 0.45 * nu / std::pow(Uei, 6) * intUe5[i];
            theta = std::sqrt(std::max(theta2, 0.0));

            double dUedx = (i + 1 < n) ? (Ue[i + 1] - Ue[i - 1]) / (s[i + 1] - s[i - 1])
                                        : (Ue[i] - Ue[i - 1]) / std::max(s[i] - s[i - 1], 1e-12);
            lambda = (theta * theta / nu) * dUedx;
            lambda = std::clamp(lambda, -0.25, 0.25); // keep correlations in their valid range

            // Skip the separation check in the last ~2% of chord: constant-
            // strength panel methods have known Ue noise right at the
            // trailing edge (Kutta-condition discretization artifact) that
            // can produce a transient dUedx spike unrelated to real flow
            // behavior. A momentary flag there is not trustworthy.
            if (lambda < ThwaitesSeparationLambda && s[i] < 0.98 * sEnd) {
                res.Separated = true;
                res.sSeparation = s[i];
            }

            double Re_theta = Uei * theta / nu;
            double Re_x = Uei * s[i] / nu;
            double Re_theta_crit = MichelReThetaCrit(std::max(Re_x, 1.0));

            if (Re_theta > Re_theta_crit && s[i] > 1e-6) {
                transIdx = i;
                transitioned = true;
                H = ThwaitesH(lambda);
                res.sTransition = s[i];
                break;
            }
        }
        if (!transitioned) {
            // Stayed laminar all the way to the trailing edge.
            H = ThwaitesH(std::clamp(lambda, -0.25, 0.25));
            res.theta_te = theta;
            res.H_te = H;
            res.Ue_te = Ue.back();
            double UeRatio = res.Ue_te / std::max(Vinf, 1e-9);
            res.cd = 2.0 * res.theta_te / chord * std::pow(UeRatio, (res.H_te + 5.0) / 2.0);
            return res;
        }
    } else {
        // Seed with a small nonzero initial momentum thickness (via the
        // same Thwaites closed form) rather than exactly theta=0: the
        // turbulent Cf correlation is singular at Re_theta=0, so even a
        // "fully turbulent from the leading edge" idealization needs a
        // numerically safe, physically tiny starting theta -- this is
        // standard practice in BL codes, not just for this implementation.
        size_t seedIdx = std::min<size_t>(2, n - 1);
        double Uei = std::max(Ue[seedIdx], 1e-9);
        double theta2 = 0.45 * nu / std::pow(Uei, 6) * intUe5[seedIdx];
        theta = std::sqrt(std::max(theta2, 1e-16));
        transIdx = seedIdx;
        H = 1.4;
    }

    // ---- Turbulent region: Head's method (RK2 marching) ----
    // State: theta, H1 (Head's entrainment shape factor). H derived from H1
    // each step via HeadHOfH1.
    double th = theta;
    double H1 = HeadH1OfH(forceTurbulentFromStart ? 1.4 : std::clamp(H, 1.1, 3.0));
    double Hcur = HeadHOfH1(H1);

    auto derivatives = [&](double s_, double th_, double H1_, double Ue_, double dUedx_) {
        double H_ = HeadHOfH1(H1_);
        double Re_theta_ = std::max(Ue_, 1e-9) * std::max(th_, 1e-9) / nu;
        double Cf_ = LudwiegTillmannCf(Re_theta_, H_);
        double dthdx_ = Cf_ / 2.0 - (H_ + 2.0) * (th_ / std::max(Ue_, 1e-9)) * dUedx_;
        double CE_ = HeadCE(H1_);
        // d(Ue*th*H1)/dx = Ue*CE  =>  dH1/dx = (Ue*CE - H1*d(Ue*th)/dx) / (Ue*th)
        // d(Ue*th)/dx = Ue*dthdx + th*dUedx
        double dUeThdx_ = Ue_ * dthdx_ + th_ * dUedx_;
        double dH1dx_ = (Ue_ * CE_ - H1_ * dUeThdx_) / std::max(Ue_ * th_, 1e-9);
        (void)s_;
        return std::make_pair(dthdx_, dH1dx_);
    };

    for (size_t i = transIdx; i + 1 < n; ++i) {
        double ds = s[i + 1] - s[i];
        if (ds < 1e-12) continue;
        double Uei = std::max(Ue[i], 1e-9);
        double Uei1 = std::max(Ue[i + 1], 1e-9);
        double dUedx = (Uei1 - Uei) / ds;

        auto [k1th, k1H1] = derivatives(s[i], th, H1, Uei, dUedx);
        double th_mid = th + 0.5 * ds * k1th;
        double H1_mid = H1 + 0.5 * ds * k1H1;
        double Ue_mid = 0.5 * (Uei + Uei1);
        auto [k2th, k2H1] = derivatives(s[i] + 0.5 * ds, th_mid, H1_mid, Ue_mid, dUedx);

        th += ds * k2th;
        H1 += ds * k2H1;
        th = std::max(th, 1e-9);
        H1 = std::max(H1, 3.01); // H1 must stay > 3 for CE correlation to be valid

        Hcur = HeadHOfH1(H1);
        if (Hcur > TurbulentSeparationH && !res.Separated && s[i + 1] < 0.98 * s.back()) {
            res.Separated = true;
            res.sSeparation = s[i + 1];
        }
    }

    res.theta_te = th;
    res.H_te = Hcur;
    res.Ue_te = Ue.back();
    double UeRatio = res.Ue_te / std::max(Vinf, 1e-9);
    res.cd = 2.0 * res.theta_te / chord * std::pow(UeRatio, (res.H_te + 5.0) / 2.0);
    return res;
}

struct AirfoilBLResult {
    SurfaceBLResult Upper, Lower;
    double cd = 0.0; // total section profile drag coefficient (upper + lower)
};

inline AirfoilBLResult AnalyzeAirfoil(const std::vector<double>& sUpper, const std::vector<double>& UeUpper,
                                       const std::vector<double>& sLower, const std::vector<double>& UeLower,
                                       double Vinf, double nu, double chord) {
    AirfoilBLResult res;
    res.Upper = AnalyzeSurface(sUpper, UeUpper, Vinf, nu, chord);
    res.Lower = AnalyzeSurface(sLower, UeLower, Vinf, nu, chord);
    res.cd = res.Upper.cd + res.Lower.cd;
    return res;
}

}} // namespace Aeolion::BoundaryLayer
