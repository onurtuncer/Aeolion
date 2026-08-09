// Solver/StabilityDerivatives.h
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

namespace Aeolion::Solver {

/** Central-difference stability & control derivatives about a baseline flight condition. */
struct StabilityDerivatives {
    double CL0 = 0, CDi0 = 0, CY0 = 0, Cm0 = 0, Croll0 = 0, Cn0 = 0;

    double CL_alpha = 0, CDi_alpha = 0, Cm_alpha = 0;      ///< Per rad.
    double CY_beta = 0, Croll_beta = 0, Cn_beta = 0;       ///< Per rad.

    double CL_q = 0, CDi_q = 0, Cm_q = 0;                  ///< Per rad/s.
    double CY_p = 0, Croll_p = 0, Cn_p = 0;                ///< Per rad/s.
    double CY_r = 0, Croll_r = 0, Cn_r = 0;                ///< Per rad/s.

    // The CONVENTIONAL nondimensional rate derivatives: the response per
    // unit reduced rate, C_x_q = dC_x/d(q cbar/(2V)) and likewise
    // dC_x/d(p b/(2V)), dC_x/d(r b/(2V)). Since the reduced rate is the
    // rate TIMES length/(2V), converting from the per-rad/s derivative
    // above divides by that factor -- i.e. multiplies by 2V/length.
    // Getting this backwards scales roll damping by (2V/b)^2, which at
    // 25 m/s on a 1 m span is a factor of two thousand: it turns a
    // textbook Cl_p = -0.45 into -0.0002 and reads as an aircraft with no
    // roll damping at all.
    double CL_q_nd = 0, Cm_q_nd = 0;                       ///< per q cbar/(2V).
    double Croll_p_nd = 0, Cn_p_nd = 0;                    ///< per p b/(2V).
    double Croll_r_nd = 0, Cn_r_nd = 0;                    ///< per r b/(2V).
};

} // namespace Aeolion::Solver
