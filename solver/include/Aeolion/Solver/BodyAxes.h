// Solver/BodyAxes.h
//
// The conversion from the solver's working frame to the CONTRACT body
// frame, and the body-axis quantities a flight-model tabulation is
// written in. This exists as a header, rather than as a helper inside
// whichever driver needs it, because the conversion has already produced
// two classes of silent defect in this repository and both are the kind a
// casual check waves through.
//
// --- the frames -----------------------------------------------------------
// The contract states its frame once and the consumer converts explicitly
// at ingest (ADR-0016). "aetherion_body_frd" is x FORWARD, y right, z
// DOWN -- the standard aeronautical body axis system of ANSI/AIAA
// R-004-1992, itself an adaptation of ISO 1151-1. The solver works in x
// aft, y right, z up. The two are related by a rotation of 180 degrees
// about y:
//
//     x_frd = -x_vlm,   y_frd = +y_vlm,   z_frd = -z_vlm
//
// That is a PROPER rotation, so handedness survives, a rotation axis
// transforms like an ordinary vector, and the right-hand rule about a
// converted axis still means what it meant before.
//
// --- trap 1: a rate derivative does NOT follow the wrench rule ------------
// For a force or moment the rule is simply "x and z flip, y does not".
// For a DERIVATIVE the flip applies to BOTH the response and the rate,
// and the two cancel whenever the pair shares x/z character:
//
//     Clp, Clr, Cnp, Cnr   moment flips, rate flips   -> INVARIANT
//     Cmq                  both about y               -> INVARIANT
//     CZq                  CZ flips, q invariant      -> flips
//     CYp, CYr             CY invariant, rate flips   -> flips
//
// Applying the wrench rule blindly negates the first four. MEASURED: it
// reported Cl_p = +0.4547 for a wing whose textbook value is -0.45 and
// which TestSolverCore pins at -0.452 -- the right magnitude with the
// wrong sign, i.e. an aircraft with roll ANTI-damping, which any
// magnitude-only check passes happily. TestBodyAxes pins the sign.
//
// --- trap 2: the coupled result's coefficient members are not filled -----
// SolveViscousCoupled reports DIMENSIONAL forces and moments; its
// SolveResult's CL/CDi/Croll/... members are left at zero, and every
// consumer forms its own coefficients from the dimensional fields. Reading
// res.Base.CL therefore yields a clean, converged, entirely zero answer.
// MEASURED: the first run of the DAVE-ML aero map produced a full alpha
// sweep of exactly zero forces at residual 1e-4. BodyAxisFromCoupled is
// the one place that assembly is written down.
#pragma once

#include "Aeolion/Math/Constants.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <vector>

namespace Aeolion::Solver {

/** Solver frame -> contract FRD. An involution: applying it twice is the identity. */
[[nodiscard]] inline Vec3 FrdFromSolver(const Vec3& v) { return Vec3(-v.x, v.y, -v.z); }

/** Body-axis (FRD) force and moment coefficients about the solve's RefPoint. */
struct BodyAxisCoefficients {
    double CX = 0.0, CY = 0.0, CZ = 0.0;
    double Cl = 0.0, Cm = 0.0, Cn = 0.0;
};

/**
 * Assemble body-axis coefficients from a Level-2 coupled result. Total
 * force is the strips' circulatory + profile forces plus the source
 * bodies' pressure integral; total moment adds the sections' own
 * quarter-chord couples. `q` is the free-stream dynamic pressure.
 *
 * Read the dimensional fields, never res.Base.CL -- see trap 2 above.
 */
[[nodiscard]] inline BodyAxisCoefficients BodyAxisFromCoupled(const ViscousCoupledResult& res,
                                                              double q,
                                                              const ReferenceGeometry& ref) {
    Vec3 force(0, 0, 0);
    for (const StripState& strip : res.Strips) force = force + strip.Force;
    force = force + res.SourceForce;
    const Vec3 moment =
        res.InducedMoment + res.ProfileMoment + res.SectionMoment + res.SourceMoment;

    const double qS = q * ref.Area;
    if (!(qS > Math::Tiny)) return {};

    const Vec3 forceFrd = FrdFromSolver(force);
    const Vec3 momentFrd = FrdFromSolver(moment);

    BodyAxisCoefficients c;
    c.CX = forceFrd.x / qS;
    c.CY = forceFrd.y / qS;
    c.CZ = forceFrd.z / qS;
    c.Cl = momentFrd.x / (qS * ref.Span);
    c.Cm = momentFrd.y / (qS * ref.Chord);
    c.Cn = momentFrd.z / (qS * ref.Span);
    return c;
}

/**
 * Reduced-rate stability derivatives in body axes, in the CONVENTIONAL
 * nondimensional form dC/d(rate * length / 2V). The reciprocal of that
 * factor is wrong by (2V/b)^2 -- a factor of two thousand at these
 * numbers -- which is why StabilityDerivatives.h carries its own warning
 * about it.
 */
struct BodyAxisRateDerivatives {
    double CZq = 0.0, Cmq = 0.0;
    double Clp = 0.0, Cnp = 0.0, CYp = 0.0;
    double Clr = 0.0, Cnr = 0.0, CYr = 0.0;
};

inline constexpr double DefaultRateStepRadps = 0.05;

/**
 * Central-difference rate derivatives on an already-factorized system.
 * The factorization depends only on geometry, so all six perturbed solves
 * reuse it. Signs follow the derivative rule of trap 1, NOT the wrench
 * rule.
 */
[[nodiscard]] inline BodyAxisRateDerivatives ComputeBodyAxisRateDerivatives(
    const PreparedSystem& prepared, const FreestreamConditions& base,
    const ReferenceGeometry& ref, double rateStep = DefaultRateStepRadps) {
    BodyAxisRateDerivatives d;
    if (!(base.Vinf > Math::Tiny)) return d;
    const double spanFactor = ref.Span / (2.0 * base.Vinf);
    const double chordFactor = ref.Chord / (2.0 * base.Vinf);

    const auto diff = [&](auto setter, double reduceFactor, auto&& assign) {
        FreestreamConditions p = base, m = base;
        setter(p, +rateStep);
        setter(m, -rateStep);
        const SolveResult rp = SolveWithSystem(prepared, p, ref);
        const SolveResult rm = SolveWithSystem(prepared, m, ref);
        assign(rp, rm, 2.0 * rateStep * reduceFactor);
    };

    diff([](FreestreamConditions& fc, double v) { fc.q = v; }, chordFactor,
         [&](const SolveResult& rp, const SolveResult& rm, double h) {
             d.CZq = -(rp.CL - rm.CL) / h; // CZ = -CL; q invariant
             d.Cmq = (rp.Cm - rm.Cm) / h;  // both about y
         });
    diff([](FreestreamConditions& fc, double v) { fc.p = v; }, spanFactor,
         [&](const SolveResult& rp, const SolveResult& rm, double h) {
             d.Clp = (rp.Croll - rm.Croll) / h; // double flip cancels
             d.Cnp = (rp.Cn - rm.Cn) / h;       // double flip cancels
             d.CYp = -(rp.CY - rm.CY) / h;      // CY invariant, p flips
         });
    diff([](FreestreamConditions& fc, double v) { fc.r = v; }, spanFactor,
         [&](const SolveResult& rp, const SolveResult& rm, double h) {
             d.Clr = (rp.Croll - rm.Croll) / h; // double flip cancels
             d.Cnr = (rp.Cn - rm.Cn) / h;       // double flip cancels
             d.CYr = -(rp.CY - rm.CY) / h;      // CY invariant, r flips
         });
    return d;
}

} // namespace Aeolion::Solver
