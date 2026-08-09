// Solver/FreestreamConditions.h
//
// The flight condition a solve is evaluated at: freestream, attitude, body
// rates, and the moment/rotation reference point.
#pragma once
#include <cmath>

#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"

namespace Aeolion::Solver {

/** Flight condition a solve is evaluated at: freestream, attitude, body rates, and moment reference point. */
struct FreestreamConditions {
    double Vinf = 20.0;      ///< [m/s]
    double alphaDeg = 5.0;   ///< Angle of attack [deg].
    double betaDeg = 0.0;    ///< Sideslip [deg] -- positive beta: relative wind has a +y (right) component.
    double rho = 1.225;      ///< [kg/m^3]
    double p = 0.0;          ///< Body roll  rate about x (aft axis)  [rad/s].
    double q = 0.0;          ///< Body pitch rate about y (right axis) [rad/s].
    double r = 0.0;          ///< Body yaw   rate about z (up axis)    [rad/s].
    Vec3 RefPoint = Vec3(0, 0, 0); ///< Moment reference point AND rotation center (e.g. CG).
};

/**
 * The TRANSLATIONAL freestream vector of a flight condition, in solver axes
 * (x aft, y right, z up).
 *
 *     Vinf = V * ( cos(alpha) cos(beta), sin(beta), sin(alpha) cos(beta) )
 *
 * Body rates are deliberately NOT included: this is the uniform oncoming
 * stream, which is what wind-axis lift/drag/side directions are built from,
 * whereas the rotational part varies from point to point (see
 * FlowField::KinematicVelocity, which adds it).
 *
 * Sign check worth keeping in mind for anything that reasons about which
 * side of a body is windward: at positive alpha the stream carries a +z
 * (upward) component in body axes, so it impinges on the LOWER surface.
 */
[[nodiscard]] inline Vec3 FreestreamVelocity(const FreestreamConditions& fc) {
    const double alpha = Math::DegToRad(fc.alphaDeg);
    const double beta = Math::DegToRad(fc.betaDeg);
    return {fc.Vinf * std::cos(alpha) * std::cos(beta), fc.Vinf * std::sin(beta),
            fc.Vinf * std::sin(alpha) * std::cos(beta)};
}

} // namespace Aeolion::Solver
