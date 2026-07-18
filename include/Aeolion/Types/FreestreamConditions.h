// Types/FreestreamConditions.h
//
// The flight condition a solve is evaluated at: freestream, attitude, body
// rates, and the moment/rotation reference point.
#pragma once
#include "Aeolion/Math/Vec3.h"

namespace Aeolion::VLM {

struct FreestreamConditions {
    double Vinf = 20.0;      // [m/s]
    double alphaDeg = 5.0;   // angle of attack [deg]
    double betaDeg = 0.0;    // sideslip [deg] -- positive beta: relative wind has a +y (right) component
    double rho = 1.225;      // [kg/m^3]
    double p = 0.0;          // body roll  rate about x (aft axis)  [rad/s]
    double q = 0.0;          // body pitch rate about y (right axis) [rad/s]
    double r = 0.0;          // body yaw   rate about z (up axis)    [rad/s]
    Vec3 RefPoint = Vec3(0, 0, 0); // moment reference point AND rotation center (e.g. CG)
};

} // namespace Aeolion::VLM
