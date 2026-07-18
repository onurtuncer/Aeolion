// Types/StabilityDerivatives.h
//
// Central-difference stability & control derivatives about a baseline
// flight condition. Angular derivatives (CL_alpha, CY_beta, ...) are per
// RADIAN. Rate derivatives are provided two ways:
//   - "dimensional": d(coefficient) / d(rate in rad/s)
//   - "_nd"         : the conventional nondimensional stability-derivative
//                      form, e.g. Cm_q_nd = dCm / d(q*cbar/(2*Vinf)), which
//                      is what you'll want for a 6-DOF sim or a DAVE-ML
//                      style derivative table.
// NOTE: CDi derivatives reflect INDUCED drag only (see SolveResult::CDi).
#pragma once

namespace Aeolion::VLM {

struct StabilityDerivatives {
    double CL0 = 0, CDi0 = 0, CY0 = 0, Cm0 = 0, Croll0 = 0, Cn0 = 0;

    double CL_alpha = 0, CDi_alpha = 0, Cm_alpha = 0;      // per rad
    double CY_beta = 0, Croll_beta = 0, Cn_beta = 0;       // per rad

    double CL_q = 0, CDi_q = 0, Cm_q = 0;                  // per rad/s
    double CY_p = 0, Croll_p = 0, Cn_p = 0;                // per rad/s
    double CY_r = 0, Croll_r = 0, Cn_r = 0;                // per rad/s

    double CL_q_nd = 0, Cm_q_nd = 0;                       // x cbar/(2V)
    double Croll_p_nd = 0, Cn_p_nd = 0;                    // x b/(2V)
    double Croll_r_nd = 0, Cn_r_nd = 0;                    // x b/(2V)
};

} // namespace Aeolion::VLM
