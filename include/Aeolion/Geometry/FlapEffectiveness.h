// Geometry/FlapEffectiveness.h
//
// Thin-airfoil-theory effectiveness of a plain trailing-edge flap: how
// much of a whole-section rotation a partial-chord flap is worth.
//
// A flap deflected by delta shifts the section's zero-lift angle by
//
//     delta_alpha0 = -tau * delta
//
// with tau the effectiveness below. That is the entire content, and it is
// why a flap belongs in the SECTION description rather than only in the
// panel geometry: to a strip-theory or lifting-line method, a deflected
// flap IS a camber change, and the natural place to express it is the
// zero-lift angle the section model is posed against.
//
// Writing x_h for the hinge position as a fraction of chord, thin-airfoil
// theory gives
//
//     theta_h = arccos(1 - 2 x_h)
//     dcl/ddelta = 2 (pi - theta_h + sin theta_h)
//     tau = (dcl/ddelta) / (dcl/dalpha) = 1 - theta_h/pi + sin(theta_h)/pi
//
// Both limits are exact and are worth stating because they are what make
// the formula checkable by inspection: a flap hinged at the leading edge
// (x_h = 0) is the whole section rotating, tau = 1; a flap of zero chord
// (x_h = 1) does nothing, tau = 0. In between tau is strongly concave --
// the contract's 12%-chord aileron is worth 0.43, not 0.12, which is the
// single most useful number here and the reason a linear-in-chord guess
// underestimates control power by a factor of three and a half.
//
// LIMITATIONS, stated rather than absorbed. This is inviscid thin-airfoil
// theory: it carries no gap leakage, no flap-edge effects, and no
// viscous decay of effectiveness at large deflection, all of which reduce
// tau in reality (a common engineering allowance is 10-20% at moderate
// deflections, more once the flap's own boundary layer separates). It
// also gives only the LIFT shift; the accompanying section pitching-moment
// increment is not modelled here.
#pragma once

#include "Aeolion/Math/Constants.h"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace Aeolion::Geometry {

/**
 * Thin-airfoil flap effectiveness tau for a plain flap whose hinge sits
 * at `hingeChordFraction` of the chord (0 = leading edge, 1 = trailing
 * edge). Returns a value in [0, 1]; out-of-range inputs clamp.
 *
 * Multiply by the deflection to get the NEGATIVE of the zero-lift-angle
 * shift: alpha0_effective = alpha0_camber - tau * deflection.
 */
[[nodiscard]] inline double FlapEffectiveness(double hingeChordFraction) {
    const double xh = std::clamp(hingeChordFraction, 0.0, 1.0);
    const double thetaH = std::acos(1.0 - 2.0 * xh);
    return std::clamp(
        1.0 - thetaH / std::numbers::pi + std::sin(thetaH) / std::numbers::pi, 0.0, 1.0);
}

/**
 * The zero-lift-angle shift, in the same angular unit as `deflection`,
 * produced by deflecting a plain flap hinged at `hingeChordFraction`.
 * Trailing edge down (positive deflection) lowers the zero-lift angle,
 * hence the sign.
 */
[[nodiscard]] inline double FlapZeroLiftShift(double hingeChordFraction, double deflection) {
    return -FlapEffectiveness(hingeChordFraction) * deflection;
}


/**
 * The section pitching-moment increment about the QUARTER CHORD from
 * deflecting a plain flap, in the same angular unit convention as
 * thin-airfoil theory (deflection in RADIANS):
 *
 *     delta_cm_c/4 = -(delta/2) * sin(theta_h) * (1 - cos(theta_h))
 *
 * This is the companion of FlapEffectiveness and comes from the same
 * place -- the Glauert coefficients of the flap's camber slope. Writing
 * the flap as dz/dx = -delta aft of the hinge gives
 * A1 = (2 delta / pi) sin(theta_h) and A2 = (delta / pi) sin(2 theta_h),
 * and cm_c/4 = (pi/4)(A2 - A1) reduces to the expression above.
 *
 * Sign: trailing edge DOWN adds aft camber, so the moment is NOSE DOWN,
 * which is why the result is negative for positive deflection.
 *
 * Both limits vanish, and both are worth stating because they make the
 * formula checkable by inspection. A flap hinged at the leading edge is
 * the whole section rotating, which adds incidence but no camber, so it
 * produces no quarter-chord moment. A flap of zero chord does nothing at
 * all. The maximum lies in between -- near the hinge position where the
 * added camber is most lopsided about the quarter chord.
 *
 * Omitting this term is not neutral. The flap's load acts aft of the
 * quarter chord, so the increment is a nose-down couple that a lift-only
 * flap model simply loses; on the contract's 12%-chord aileron at 10
 * degrees it is delta_cm = -0.100 against delta_cl = 0.476, i.e. the
 * increment acts about 0.21 chords aft of the quarter chord.
 */
[[nodiscard]] inline double FlapMomentIncrement(double hingeChordFraction,
                                                double deflectionRad) {
    const double xh = std::clamp(hingeChordFraction, 0.0, 1.0);
    const double thetaH = std::acos(1.0 - 2.0 * xh);
    return -Math::Half * deflectionRad * std::sin(thetaH) * (1.0 - std::cos(thetaH));
}

} // namespace Aeolion::Geometry
