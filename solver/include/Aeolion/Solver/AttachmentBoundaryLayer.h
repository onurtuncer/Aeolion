// Solver/AttachmentBoundaryLayer.h
//
// Where the boundary layer separates, marched from the attachment point
// that Solver/AttachmentLine.h locates rather than from the camber line's
// leading edge.
//
// --- what this is, and what it is not ---------------------------------------
// This is a DIAGNOSTIC march, not a coupling. It consumes a converged
// attachment line and reports where each strip's surface flow separates;
// it feeds nothing back into the lattice and no fixed point iterates on
// it. That distinction is why it may do something
// Solver/SectionBoundaryLayer.h deliberately must not: LATCH the first
// separation crossing and stop.
//
// SectionBoundaryLayer marches inside two nested fixed points, so it ramps
// the skin friction smoothly to zero across [SeparationRampH,
// TurbulentSeparationH] and never latches -- a latched flag would make the
// marched state hysteretic, and hysteresis inside the section is what feeds
// limit cycles in the coupling loops (theory.rst, "stabilization"). Nothing
// iterates on THIS march, so the honest first crossing is the right answer
// and the smoothing would only blur the location being reported.
//
// The correlations themselves are shared with that module rather than
// restated -- Thwaites' shape factor and shear, Michel's transition
// criterion, Head's entrainment with Ludwieg-Tillmann skin friction -- so
// the two marches cannot drift apart in their physics, only in what they
// do at the crossing.
//
// --- why starting at the attachment point matters ----------------------------
// SectionBoundaryLayer starts both surfaces at the camber line's leading
// edge with theta = 0, which assumes the stagnation point sits at
// psi = 0. It does not: a cambered section stagnates below its leading
// edge even at zero incidence, and the offset grows like
// sqrt(r_LE) * alpha_e (Section "Where the stagnation point lies" in the
// paper). Starting the march at the real attachment point costs nothing
// here because SectionSolution::UpperRun() already measures its arc length
// from that point.
//
// Thwaites supplies its own initial condition when the integral starts at
// the stagnation point, which is worth stating because it looks like a
// missing input. Near the attachment point Ue ~ a s, so
//
//     theta^2 = 0.45/(Re Ue^6) * integral_0^s Ue^5 ds
//             -> 0.45 a^5 s^6 / (6 Re a^6 s^6)  =  0.075 / (Re a),
//
// which is exactly the stagnation momentum thickness
// AttachmentStation::MomentumThickness reports independently. The march
// therefore needs no seed, and the agreement between the two is a
// consistency check the test asserts rather than a coincidence.
//
// --- units -------------------------------------------------------------------
// A SurfaceRun states arc length in fractions of the NORMAL chord and edge
// speed in multiples of the NORMAL speed (the leading-edge-normal plane is
// the one the section problem is posed in -- see AttachmentLine.h), so the
// Reynolds number this march takes is likewise the normal-plane one,
// Re_n = U_n c_n / nu. Everything returned is in those same nondimensional
// terms.

#pragma once

#include "Aeolion/Math/Constants.h"
#include "Aeolion/Solver/AttachmentLine.h"
#include "Aeolion/Solver/SectionBoundaryLayer.h"
#include "Aeolion/Solver/SectionPanelMethod.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace Aeolion::Solver {

// --- numerical guards --------------------------------------------------------
/** A run with fewer stations than this cannot resolve a gradient, let alone a crossing. */
inline constexpr std::size_t MinMarchStations = 4;
/** Below this edge speed the station is inside the stagnation region; Ue^6 there is noise. */
inline constexpr double MarchEdgeSpeedFloor = 1e-6;
/** Re_s below this makes Michel's correlation meaningless (it is a fit, not a limit). */
inline constexpr double MinMichelReynolds = 1.0;

// --- what counts as separation ----------------------------------------------
// At the Reynolds number of a small airframe (Re_n ~ 3e5 here) the laminar
// layer ALWAYS reaches Thwaites' lambda = -0.09 before Michel's criterion
// trips, at every incidence including negative ones. Reporting that as
// "the wing has separated" is true of the laminar layer and useless as an
// answer: it says every attitude is separated and draws no boundary.
//
// What happens physically is a laminar separation BUBBLE -- the sheared
// layer transitions just downstream of separating and reattaches turbulent
// a few percent chord later. That is the defining feature of low-Reynolds
// airfoils, not a failure of the flow. So the bubble is treated as the
// transition trigger (the standard engineering treatment, and the one
// Solver/SectionBoundaryLayer.h already uses inside its fixed point), the
// march continues turbulent, and the VERDICT is turbulent separation:
// the trailing-edge separation point moving forward, which is what stalls
// a wing and what an attitude limit should be drawn from.
//
// The laminar separation point is still reported, because it is the part
// of this march that is verifiable against a closed-form answer (Howarth's
// retarded flow, s/L = 0.123) and because a designer wants to know where
// the bubble sits.
//
// STATED LIMITATION: bubble BURSTING is not modelled. A short bubble that
// fails to reattach is what actually ends the lift curve of a thin section
// at Re below ~2e5, and predicting it needs a bubble-length correlation
// (Gaster) or an e^N envelope this method does not carry. The turbulent
// separation boundary reported here is therefore an upper bound on the
// usable incidence, not a stall prediction.

/**
 * How a surface run ended. There is deliberately no "laminar separation"
 * state: laminar separation triggers transition here rather than ending
 * the march, so it is reported through SurfaceMarch::BubbleFormed and its
 * location, not as an outcome.
 */
enum class SeparationMode {
    Attached,            ///< Reached the trailing edge attached (possibly via a bubble).
    TurbulentSeparation, ///< Head's shape factor reached 2.4 -- the separation verdict.
};

/** One surface of one strip, marched from its attachment point. */
struct SurfaceMarch {
    bool Valid = false;
    SeparationMode Mode = SeparationMode::Attached;

    [[nodiscard]] bool Separated() const { return Mode != SeparationMode::Attached; }

    /**
     * Arc length from the attachment point at which the layer separates, in
     * fractions of the normal chord. When the run stays attached this is
     * the run's full length -- so "how far did the layer get" is one field
     * whether or not it separated.
     */
    double SeparationArc = 0.0;
    /**
     * Chordwise station psi = x/c_n of the separation point. This is the
     * number a reader wants: "separates at 0.62 chord". Equal to the run's
     * last psi when attached.
     */
    double SeparationPsi = 1.0;

    bool Transitioned = false;
    double TransitionArc = 0.0;
    double TransitionPsi = 0.0;
    /**
     * Transition was triggered by laminar separation (a bubble) rather than
     * by Michel's criterion -- so TransitionArc/Psi is also the laminar
     * separation point. Worth distinguishing: a bubble sitting near the
     * leading edge at high incidence is the thing that bursts, and this
     * march does not model bursting (see the note above).
     */
    bool BubbleFormed = false;

    /**
     * The momentum thickness the Thwaites integral produces at the first
     * station, theta_0/c_n. Compare against
     * AttachmentStation::MomentumThickness / NormalChord: they are the same
     * quantity reached two different ways.
     */
    double StartingTheta = 0.0;
    double TrailingTheta = 0.0;      ///< theta/c_n where the march ended.
    double TrailingShapeFactor = 0.0;
    double ReynoldsNormal = 0.0;     ///< The Re_n this march ran at.

    // --- the marched profile ---------------------------------------------------
    // Kept because a separation location is only as trustworthy as the
    // distribution behind it, and a caller plotting cf(s) or H(s) should not
    // have to re-run the march to get them. Aligned with the run's stations.
    std::vector<double> S, Ue, Theta, ShapeFactor;
};

namespace Detail {

/**
 * Locate a crossing WITHIN a step by linear interpolation of a signed
 * indicator, rather than quantizing it to the station that first tripped.
 * A separation point rounded to whole stations moves in visible jumps as
 * alpha varies, which reads as physics and is discretization.
 */
[[nodiscard]] inline double CrossingFraction(double previous, double current, double threshold) {
    const double span = current - previous;
    if (!(std::fabs(span) > Math::Tiny)) return 0.0;
    return std::clamp((threshold - previous) / span, 0.0, 1.0);
}

/**
 * Exact integral of a piecewise-LINEAR edge velocity's fifth power over one
 * step, which is what Thwaites' integral wants and what the trapezoid rule
 * gets badly wrong at exactly the station that matters most.
 *
 * On the FIRST step out of an attachment point Ue rises linearly from zero,
 * so the trapezoid rule returns h Ue1^5 / 2 where the true value is
 * h Ue1^5 / 6 -- a factor of three in the integral, hence a factor sqrt(3)
 * in theta_0, the one number a march from a stagnation point exists to get
 * right. (Starting from a leading edge with theta = 0 hides this, since
 * there is no theta_0 to be wrong about.)
 *
 * For a linear u over [0, h],  integral u^5 = h (u1^6 - u0^6) / (6 (u1 - u0)),
 * which tends to h u0^5 as u1 -> u0. Ue between adjacent stations is linear
 * to the accuracy of the panelling, so this uses the known local behaviour
 * rather than fitting one; the near-equal branch avoids the cancellation in
 * the difference of sixth powers.
 */
[[nodiscard]] inline double FifthPowerIntegral(double u0, double u1, double step) {
    const double difference = u1 - u0;
    const double scale = std::max(std::fabs(u0), std::fabs(u1));
    if (!(std::fabs(difference) > 1e-9 * std::max(scale, Math::Tiny)))
        return step * std::pow(Math::Half * (u0 + u1), 5);
    return step * (std::pow(u1, 6) - std::pow(u0, 6)) / (6.0 * difference);
}

} // namespace Detail

/**
 * March one surface run from its attachment point and report where it
 * separates.
 *
 * `reynoldsNormal` is U_n c_n / nu, matching the run's own
 * nondimensionalization. Thwaites runs until either Michel's criterion or
 * laminar separation trips; a laminar separation is reported and the march
 * STOPS, because a laminar separation bubble that reattaches turbulent is a
 * different physical claim than the attached layer this method models, and
 * SectionBoundaryLayer's treatment of it (transition at the bubble) exists
 * to keep a fixed point smooth, not because the bubble is resolved.
 */
[[nodiscard]] inline SurfaceMarch MarchSurfaceRun(const SurfaceRun& run, double reynoldsNormal) {
    SurfaceMarch march;
    const std::size_t n = run.Count();
    if (n < MinMarchStations || !(reynoldsNormal > 0.0)) return march;

    march.Valid = true;
    march.ReynoldsNormal = reynoldsNormal;
    march.S.reserve(n);
    march.Ue.reserve(n);
    march.Theta.reserve(n);
    march.ShapeFactor.reserve(n);

    const double Re = reynoldsNormal;
    bool turbulent = false;
    double theta = 0.0, H = 2.61, H1 = Detail::HeadH1(TurbulentInitialH);
    double thwaitesIntegral = 0.0;

    // The run's arc starts AT the attachment point, where Ue = 0, so the
    // integral's lower limit is a genuine zero rather than an extrapolation.
    double sPrev = 0.0, uePrev = 0.0;
    double lambdaPrev = 0.0, michelExcessPrev = -1e30, thetaLamPrev = 0.0, psiPrev = run.Psi.front();

    const auto advanceTurbulent = [&](double dsLen, double ue, double dueds) {
        if (dsLen <= 0.0) return;
        const double reTheta = std::max(Re * ue * theta, 20.0);
        const double cf = LudwiegTillmannA * std::pow(10.0, -LudwiegTillmannB * H) *
                          std::pow(reTheta, -LudwiegTillmannC);
        const double dThetaDs = Math::Half * cf - (H + 2.0) * theta / ue * dueds;
        const double entrainment =
            HeadEntrainmentA * std::pow(std::max(H1 - 3.0, 1e-3), -HeadEntrainmentB);
        const double dH1Ds =
            entrainment / std::max(theta, 1e-9) - H1 * (dThetaDs / std::max(theta, 1e-9) + dueds / ue);
        theta = std::max(theta + dThetaDs * dsLen, 1e-9);
        H1 = std::max(H1 + dH1Ds * dsLen, 3.32);
        H = Detail::HeadHFromH1(H1);
    };

    for (std::size_t i = 0; i < n; ++i) {
        const double s = run.S[i];
        const double psi = run.Psi[i];
        const double dsStep = std::max(s - sPrev, 1e-12);
        const double ue = std::max(run.Ue[i], MarchEdgeSpeedFloor);
        const double dueds = (ue - uePrev) / dsStep;

        if (!turbulent) {
            thwaitesIntegral += Detail::FifthPowerIntegral(uePrev, ue, dsStep);
            const double theta2 = ThwaitesConstant / (Re * std::pow(ue, 6)) * thwaitesIntegral;
            theta = std::sqrt(std::max(theta2, 1e-16));
            const double lambda = Re * theta2 * dueds;
            H = Detail::ThwaitesH(lambda);
            if (i == 0) march.StartingTheta = theta;

            const double reTheta = Re * ue * theta;
            const double reS = Re * ue * s;
            const double excess = (reS > MinMichelReynolds)
                                      ? reTheta - MichelFactor * (1.0 + MichelOffset / reS) *
                                                      std::pow(reS, MichelExponent)
                                      : -1e30;

            // Either trigger transitions the layer; whichever crosses
            // earlier within the step wins. Laminar separation is a bubble,
            // not the end of the march (see the note above).
            bool trigger = false, bubble = false;
            double frac = 1.0;
            if (excess >= 0.0) {
                trigger = true;
                frac = (michelExcessPrev < -1e29)
                           ? 0.0
                           : Detail::CrossingFraction(michelExcessPrev, excess, 0.0);
            }
            if (lambda < LaminarSeparationLambda) {
                const double fracSep =
                    (i == 0) ? 0.0
                             : Detail::CrossingFraction(lambdaPrev, lambda, LaminarSeparationLambda);
                if (!trigger || fracSep < frac) bubble = true;
                frac = trigger ? std::min(frac, fracSep) : fracSep;
                trigger = true;
            }

            if (trigger) {
                march.Transitioned = true;
                march.BubbleFormed = bubble;
                march.TransitionArc = sPrev + frac * (s - sPrev);
                march.TransitionPsi = psiPrev + frac * (psi - psiPrev);
                theta = std::max(thetaLamPrev + frac * (theta - thetaLamPrev), 1e-9);
                turbulent = true;
                H = TurbulentInitialH;
                H1 = Detail::HeadH1(H);
                advanceTurbulent((1.0 - frac) * dsStep, ue, dueds);
                if (H >= TurbulentSeparationH) {
                    march.Mode = SeparationMode::TurbulentSeparation;
                    march.SeparationArc = s;
                    march.SeparationPsi = psi;
                    march.TrailingTheta = theta;
                    march.TrailingShapeFactor = H;
                    march.S.push_back(s);
                    march.Ue.push_back(ue);
                    march.Theta.push_back(theta);
                    march.ShapeFactor.push_back(H);
                    return march;
                }
            }

            michelExcessPrev = excess;
            lambdaPrev = lambda;
            thetaLamPrev = theta;
        } else {
            const double HPrev = H;
            advanceTurbulent(dsStep, ue, dueds);
            if (H >= TurbulentSeparationH) {
                const double frac = Detail::CrossingFraction(HPrev, H, TurbulentSeparationH);
                march.Mode = SeparationMode::TurbulentSeparation;
                march.SeparationArc = sPrev + frac * (s - sPrev);
                march.SeparationPsi = psiPrev + frac * (psi - psiPrev);
                march.TrailingTheta = theta;
                march.TrailingShapeFactor = H;
                march.S.push_back(s);
                march.Ue.push_back(ue);
                march.Theta.push_back(theta);
                march.ShapeFactor.push_back(H);
                return march;
            }
        }

        march.S.push_back(s);
        march.Ue.push_back(ue);
        march.Theta.push_back(theta);
        march.ShapeFactor.push_back(H);

        sPrev = s;
        uePrev = ue;
        psiPrev = psi;
    }

    // Reached the trailing edge attached.
    march.SeparationArc = run.Length();
    march.SeparationPsi = run.Psi.back();
    march.TrailingTheta = theta;
    march.TrailingShapeFactor = H;
    return march;
}

/** One strip's two surfaces, marched. */
struct StationSeparation {
    int Strip = -1;
    double Eta = 0.0;      ///< Semi-span fraction, unsigned (as StripSection states it).
    bool Resolved = false; ///< The strip had an attachment point to march from.
    SurfaceMarch Upper;    ///< The suction side -- the one that separates first.
    SurfaceMarch Lower;
};

/**
 * The separation picture across a lifting surface at one flight condition.
 *
 * `SeparatedStations` counts UPPER-surface separations only. That is not an
 * omission: at positive incidence the lower surface runs a favourable
 * gradient nearly to the trailing edge and separating there would mean
 * something else entirely had gone wrong, so it is marched and reported but
 * does not drive the verdict.
 */
struct SeparationSurvey {
    std::vector<StationSeparation> Stations;
    int ResolvedStations = 0;
    int SeparatedStations = 0;
    /** Most forward upper-surface separation anywhere on the wing, as psi = x/c_n. */
    double ForwardmostSeparationPsi = 1.0;
    double ForwardmostSeparationEta = 0.0;

    [[nodiscard]] bool AnySeparated() const { return SeparatedStations > 0; }
    /** Fraction of resolved stations whose upper surface separates. */
    [[nodiscard]] double SeparatedFraction() const {
        return (ResolvedStations > 0)
                   ? static_cast<double>(SeparatedStations) / static_cast<double>(ResolvedStations)
                   : 0.0;
    }
};

/**
 * March every station of a converged attachment line and report where the
 * wing separates.
 *
 * Stations whose attachment point was not resolved are carried with
 * `Resolved == false` rather than dropped, so the survey's indices still
 * line up with the attachment line's.
 */
[[nodiscard]] inline SeparationSurvey SurveySeparation(const AttachmentLine& line,
                                                       double kinematicViscosity =
                                                           SeaLevelKinematicViscosity) {
    SeparationSurvey survey;
    survey.Stations.reserve(line.Stations.size());

    for (const AttachmentStation& station : line.Stations) {
        StationSeparation entry;
        entry.Strip = station.Strip;
        entry.Eta = station.Eta;

        if (!station.Found || !station.Section.Valid || !(station.NormalChord > 0.0) ||
            !(station.NormalSpeed > 0.0) || !(kinematicViscosity > 0.0)) {
            survey.Stations.push_back(entry);
            continue;
        }

        const double Re = station.NormalSpeed * station.NormalChord / kinematicViscosity;
        entry.Upper = MarchSurfaceRun(station.Section.UpperRun(), Re);
        entry.Lower = MarchSurfaceRun(station.Section.LowerRun(), Re);
        entry.Resolved = entry.Upper.Valid;

        if (entry.Resolved) {
            ++survey.ResolvedStations;
            if (entry.Upper.Separated()) {
                ++survey.SeparatedStations;
                if (entry.Upper.SeparationPsi < survey.ForwardmostSeparationPsi) {
                    survey.ForwardmostSeparationPsi = entry.Upper.SeparationPsi;
                    survey.ForwardmostSeparationEta = station.Eta;
                }
            }
        }
        survey.Stations.push_back(entry);
    }
    return survey;
}

} // namespace Aeolion::Solver
