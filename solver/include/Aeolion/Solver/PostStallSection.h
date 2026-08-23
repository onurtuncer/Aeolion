// Solver/PostStallSection.h
//
// Phase 1 of the post-separation study: a post-stall section model whose
// every ingredient is either classical theory with published validation or
// the solver's own computed separation state -- replacing the four
// hand-tuned constants of AnalyticSectionModel's deep-stall blend, which is
// exactly what TODO.md 3b's caveat demanded before any post-stall number
// could be quoted.
//
// The construction, from attached flow to 90 degrees:
//
//   1. ATTACHED / PARTIALLY SEPARATED (Kirchhoff).  Free-streamline theory
//      gives the lift attenuation of a section whose suction-side flow
//      separates at chord fraction f (1 = trailing edge, 0 = leading edge)
//      as K(f) = ((1 + sqrt(f)) / 2)^2 -- the factor the Beddoes-Leishman
//      family builds its static backbone on (leishman2006rotor). Here it
//      attenuates the exact thin-airfoil normal/chordwise pair,
//
//          cn = cl_alpha * sin(a) * cos(a) * K(f),
//          cc = cl_alpha * sin^2(a) * sqrt(f),
//
//      chosen so BOTH classical limits are exact: f = 1 recovers
//      cl = cl_alpha * sin(a) with identically zero pressure drag
//      (d'Alembert), f = 0 leaves the force perpendicular to the chord.
//      The separation point f is NOT modelled here: it is supplied by the
//      caller, computed by Solver/AttachmentBoundaryLayer.h's march from
//      the real attachment point. Stall is therefore an OUTPUT -- the
//      incidence where the branch's lift peaks -- not a ClMax constant.
//
//   2. DEEP STALL (Viterna).  Past the computed stall angle the Viterna-
//      Corrigan extension carries the polar to 90 degrees anchored on the
//      finite-wing bluff-plate drag
//
//          CdMax = 1.11 + 0.018 * AR      (viternaCorrigan1982),
//
//      continuous with the Kirchhoff branch at the junction by
//      construction. Hoerner's 2-D post-stall normal force
//      cn = 1/(0.222 + 0.283/sin(alpha)) (hoerner1965fluiddynamicdrag) is
//      carried as the infinite-AR cross-check: CdMax(AR -> 50) = 2.01 vs
//      Hoerner's 1.98 at 90 degrees.
//
//   3. CENTRE OF PRESSURE (Rayleigh).  Kirchhoff flow past an inclined
//      plate has the closed-form centre of pressure
//
//          x_cp / c = 1/2 - (3/4) cos(alpha) / (4 + pi sin(alpha)),
//
//      5/16 at small incidence, exactly mid-chord at 90. The section cm
//      (about the quarter chord) interpolates between the attached
//      quarter-chord (cm = 0) and the Rayleigh point with the separated
//      fraction (1 - f) -- this is what lets the 3-D solve's centre of
//      pressure finally walk aft of c/4, which the Phase-0 map proved
//      unrepresentable without a section cm.
//
// Known limitations, stated rather than hidden: the f-table is queried at
// |alpha from zero lift|, so camber asymmetry of separation between upper
// and lower surface is not represented; attached-camber cm0 is not carried
// (the lattice geometry supplies camber through Alpha0Deg only); Kirchhoff
// theory's missing base suction is exactly why the deep end is anchored on
// Viterna rather than on Kirchhoff's own CD(90) = 0.88.

#pragma once

#include "Aeolion/Math/Constants.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <numbers>

namespace Aeolion::Solver {

// --- Viterna-Corrigan --------------------------------------------------------
/** The AR fit's stated validity edge; beyond it the 2-D plate value is held. */
inline constexpr double ViternaMaxAspectRatio = 50.0;

/** Finite-wing bluff-plate drag ceiling, CdMax = 1.11 + 0.018 AR. */
[[nodiscard]] inline double ViternaCdMax(double aspectRatio) {
    return 1.11 + 0.018 * std::clamp(aspectRatio, 0.0, ViternaMaxAspectRatio);
}

/** The Viterna extension's two matching constants, fixed by continuity at the stall junction. */
struct ViternaConstants {
    double CdMax = 2.01;
    double KL = 0.0; ///< cl = CdMax/2 sin(2a) + KL cos^2(a)/sin(a)
    double KD = 0.0; ///< cd = CdMax sin^2(a) + KD cos(a)
};

/**
 * Match the Viterna extension to the polar at the stall junction
 * (alphaStallRad, clStall, cdStall) so both curves are continuous there.
 */
[[nodiscard]] inline ViternaConstants MatchViterna(double aspectRatio, double alphaStallRad,
                                                   double clStall, double cdStall) {
    ViternaConstants v;
    v.CdMax = ViternaCdMax(aspectRatio);
    const double sa = std::sin(alphaStallRad), ca = std::cos(alphaStallRad);
    v.KL = (clStall - v.CdMax * sa * ca) * sa / std::max(ca * ca, Math::Tiny);
    v.KD = (cdStall - v.CdMax * sa * sa) / std::max(ca, Math::Tiny);
    return v;
}

[[nodiscard]] inline double ViternaCl(const ViternaConstants& v, double alphaRad) {
    const double sa = std::sin(alphaRad), ca = std::cos(alphaRad);
    return Math::Half * v.CdMax * std::sin(2.0 * alphaRad) +
           v.KL * ca * ca / std::max(sa, 1e-6);
}

[[nodiscard]] inline double ViternaCd(const ViternaConstants& v, double alphaRad) {
    const double sa = std::sin(alphaRad);
    return v.CdMax * sa * sa + v.KD * std::cos(alphaRad);
}

// --- classical anchors -------------------------------------------------------
/** Hoerner's 2-D post-stall normal force, cn = 1/(0.222 + 0.283/sin a): 1.98 at 90 deg. */
[[nodiscard]] inline double HoernerPostStallCn(double alphaRad) {
    const double sa = std::max(std::fabs(std::sin(alphaRad)), 1e-6);
    return 1.0 / (0.222 + 0.283 / sa);
}

/** Kirchhoff lift attenuation ((1 + sqrt(f))/2)^2: 1 attached, 1/4 fully separated. */
[[nodiscard]] inline double KirchhoffAttenuation(double f) {
    const double root = Math::Half * (1.0 + std::sqrt(std::clamp(f, 0.0, 1.0)));
    return root * root;
}

/** Rayleigh's centre of pressure for Kirchhoff flow past an inclined plate. */
[[nodiscard]] inline double RayleighCenterOfPressure(double alphaRad) {
    const double sa = std::fabs(std::sin(alphaRad));
    return Math::Half - 0.75 * std::fabs(std::cos(alphaRad)) / (4.0 + std::numbers::pi * sa);
}

// ------------------------------------------------------------------ the model
/**
 * How far the suction-side flow stays attached, f = psi_sep in [0, 1], at
 * span fraction eta and incidence |alpha - alpha_0| in degrees. Supplied by
 * the caller from AttachmentBoundaryLayer's computed separation survey
 * (interpolated per strip); an empty function means fully attached until
 * the fallback stall scan finds no peak, i.e. no stall below the scan edge.
 */
using SeparationPointFunction = std::function<double(double eta, double alphaFromZeroLiftDeg)>;

// The stall-junction scan: the Kirchhoff branch's lift peaks where
// d/da [sin(a) K(f(a))] = 0, located by a coarse-to-exact scan over this
// range. Past the scan edge without a peak, the edge is the junction.
inline constexpr double StallScanStartDeg = 4.0;
inline constexpr double StallScanEndDeg = 45.0;
inline constexpr double StallScanStepDeg = 0.5;

/**
 * Anchored post-stall section model (SectionModel signature): Kirchhoff
 * attenuation from the COMPUTED separation point up to the emergent stall,
 * Viterna-Corrigan with AR-aware CdMax beyond, Rayleigh centre of pressure
 * for the section cm. See the header comment for the construction and its
 * stated limitations.
 */
struct PostStallSectionModel {
    double ClAlphaPerRad = Math::Two * std::numbers::pi;
    double AspectRatio = 6.0;
    double Cd0 = 0.012;             ///< At ReferenceReynolds (same polar constants as
    double KCd = 0.015;             ///< AnalyticSectionModel, so the attached limit agrees).
    double ReferenceReynolds = 2e5;
    double ReynoldsExponent = 0.2;
    double SuctionRecovery = 1.0;   ///< eta_e on the chordwise force (Beddoes-Leishman's ~0.95-1).
    SeparationPointFunction SeparationPoint;

    [[nodiscard]] SectionCoefficients operator()(const StripSection& strip, double alphaEffDeg,
                                                 double Re, double /*Ma*/) const {
        const double alphaZDeg = alphaEffDeg - strip.EffectiveAlpha0Deg();
        const double sign = (alphaZDeg >= 0.0) ? 1.0 : -1.0;
        const double aDeg = std::fabs(alphaZDeg);
        const double aRad = Math::DegToRad(aDeg);
        const double reScale =
            (Re > 0.0) ? std::pow(ReferenceReynolds / Re, ReynoldsExponent) : 1.0;
        const double cdVisc0 = Cd0 * reScale;

        // The emergent stall junction for THIS strip's f(alpha).
        const double aStallDeg = StallAngleDeg(strip.Eta);

        SectionCoefficients out;
        double cn = 0.0;
        if (aDeg <= aStallDeg) {
            const auto branch = KirchhoffBranch(strip.Eta, aRad);
            cn = branch.cn;
            out.cl = branch.cl;
            out.cd = cdVisc0 + KCd * branch.cl * branch.cl + branch.cdPressure;
        } else {
            // Viterna, matched where the Kirchhoff branch left off (with its
            // viscous drag included so cd is continuous too).
            const double asRad = Math::DegToRad(aStallDeg);
            const auto junction = KirchhoffBranch(strip.Eta, asRad);
            const double cdJunction =
                cdVisc0 + KCd * junction.cl * junction.cl + junction.cdPressure;
            const ViternaConstants v =
                MatchViterna(AspectRatio, asRad, junction.cl, cdJunction);
            out.cl = ViternaCl(v, aRad);
            out.cd = ViternaCd(v, aRad);
            cn = out.cl * std::cos(aRad) + (out.cd - cdVisc0) * std::sin(aRad);
        }

        // Centre of pressure: quarter chord while attached, the Rayleigh
        // point as the separated fraction takes over. cm about c/4.
        const double f = SeparationFraction(strip.Eta, aDeg);
        const double xcp = f * Math::QuarterChord + (1.0 - f) * RayleighCenterOfPressure(aRad);
        out.cm = -(xcp - Math::QuarterChord) * cn;

        out.cl *= sign;
        out.cm *= sign;
        return out;
    }

    /** The emergent stall angle (deg from zero lift) at span fraction eta. */
    [[nodiscard]] double StallAngleDeg(double eta) const {
        double bestDeg = StallScanEndDeg;
        double bestCl = -1.0;
        for (double aDeg = StallScanStartDeg; aDeg <= StallScanEndDeg + 1e-9;
             aDeg += StallScanStepDeg) {
            const double cl = KirchhoffBranch(eta, Math::DegToRad(aDeg)).cl;
            if (cl > bestCl) {
                bestCl = cl;
                bestDeg = aDeg;
            }
        }
        return bestDeg;
    }

private:
    [[nodiscard]] double SeparationFraction(double eta, double alphaFromZeroLiftDeg) const {
        if (!SeparationPoint) return 1.0;
        return std::clamp(SeparationPoint(eta, alphaFromZeroLiftDeg), 0.0, 1.0);
    }

    struct BranchState {
        double cl = 0.0, cn = 0.0, cdPressure = 0.0;
    };

    // The attenuated thin-airfoil pair (see the header comment): with
    // K(1) = 1 the wind-axis rotation gives cl = cl_alpha sin(a) and zero
    // pressure drag EXACTLY, so the attached limit needs no special case.
    [[nodiscard]] BranchState KirchhoffBranch(double eta, double aRad) const {
        const double f = SeparationFraction(eta, Math::RadToDeg(aRad));
        const double sa = std::sin(aRad), ca = std::cos(aRad);
        BranchState b;
        b.cn = ClAlphaPerRad * sa * ca * KirchhoffAttenuation(f);
        const double cc = SuctionRecovery * ClAlphaPerRad * sa * sa * std::sqrt(f);
        b.cl = b.cn * ca + cc * sa;
        b.cdPressure = std::max(b.cn * sa - cc * ca, 0.0);
        return b;
    }
};

} // namespace Aeolion::Solver
