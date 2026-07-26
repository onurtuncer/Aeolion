// BEMT/BEMT.h
//
// Blade Element Momentum Theory (BEMT) for a propeller, plus a slipstream
// velocity field meant to be fed into Solver.h's Solve() externalField hook
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
// This header declares the public types and functions; BEMT.cpp holds the
// implementation. Separated from the header-only `aeolion` target -- like
// PanelBuilder -- because it is a self-contained solver with real
// implementation weight, not a thin data-vocabulary header consumers only
// need to construct/read (that's still true of PropGeometry.h, which is
// why it stays in the header-only core: see PropGeometry.h and
// Geometry/HandoffContract.h, which builds a PropGeometry but never calls
// into the solver itself).

#pragma once
#include <numbers>
#include <vector>
#include "Aeolion/BEMT/PropGeometry.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"

namespace Aeolion::BEMT {

using Solver::Vec3;

// --- polar defaults -------------------------------------------------------
inline constexpr double ThinAirfoilClAlpha = Math::Two * std::numbers::pi; // per rad
inline constexpr double DefaultClMax = 1.2;
inline constexpr double DefaultCd0   = 0.02;
inline constexpr double DefaultKCd   = 0.02;   // cd = cd0 + kCd*cl^2

// --- solver defaults ------------------------------------------------------
inline constexpr double SeaLevelDensity   = 1.225; // kg/m^3, ISA
inline constexpr int    DefaultMaxIter    = 3000;
inline constexpr double DefaultTolerance  = 1e-7;
inline constexpr double DefaultRelaxation = 0.25;

// --- slipstream field -------------------------------------------------------
inline constexpr double FarWakeDiametersAuto = 2.0; // auto development length = this * radius

// ---------------------------------------------------------- Airfoil polar
// Analytic default (smooth-saturating lift curve + parabolic drag polar)
// or an optional table lookup (e.g. from measured or CFD section data) --
// set UseTable=true and fill TableAlphaDeg/TableCl/TableCd (sorted by
// alpha) to use real section data instead.
/**
 * Section lift/drag polar: analytic default (smooth-saturating lift curve +
 * parabolic drag polar), or a table lookup when UseTable is set (fill
 * TableAlphaDeg/TableCl/TableCd, sorted by alpha, with real section data).
 */
struct Polar {
    double cl_alpha = ThinAirfoilClAlpha; ///< Per rad, thin-airfoil default.
    double alpha0Deg = 0.0;     ///< Zero-lift angle.
    double clMax = DefaultClMax;
    double cd0 = DefaultCd0;
    double kCd = DefaultKCd;

    bool UseTable = false;
    std::vector<double> TableAlphaDeg, TableCl, TableCd;
};

/** Evaluate a Polar at alphaDeg, writing the resulting cl/cd. */
void EvalPolar(const Polar& p, double alphaDeg, double& cl, double& cd);

// Blade geometry lives in BEMT/PropGeometry.h so that describing a
// propeller does not require including this solver.

// -------------------------------------------------------------- Solving
/** BEMT solution at one radial blade station. */
struct StationResult {
    double r = 0, vi = 0, wi = 0, phiDeg = 0, alphaDeg = 0, cl = 0, cd = 0;
    double dT_dr = 0, dQ_dr = 0;

    /**
     * Per-station convergence. Result::Converged is the AND of these, which
     * on its own cannot say whether one awkward station failed or the whole
     * blade did -- and those call for very different responses.
     */
    bool Converged = true;
    int Iterations = 0;
    double Residual = 0.0; ///< Last |change| in induced velocity [m/s].
};

/** Full BEMT solve result: per-station breakdown plus integrated thrust/torque/power. */
struct Result {
    std::vector<StationResult> Stations;
    double Thrust = 0, Torque = 0, Power = 0;
    double rpm = 0, omega = 0;
    PropGeometry Geom;
    bool Converged = true;
};

/** Prandtl tip/hub loss factor F at radius rIn and inflow angle phi. */
[[nodiscard]] double TipHubLoss(const PropGeometry& g, double rIn, double phi);

/**
 * Solve the hover-safe BEMT fixed point for every radial station of geom,
 * at rotor speed rpm and freestream Vinf, then integrate thrust/torque/power.
 */
[[nodiscard]] Result Solve(const PropGeometry& geom, const Polar& polar, double rpm, double Vinf,
                           double rho = SeaLevelDensity, int maxIter = DefaultMaxIter,
                           double tol = DefaultTolerance, double relax = DefaultRelaxation);

// ------------------------------------------------------ Performance metrics
// The two hard bounds a propeller result must respect. Both are ratios of
// an ideal power to the power actually absorbed, so both are <= 1 for any
// physically sensible answer -- they are thermodynamic constraints, not
// tuning targets. TestBEMT checks exactly this, and it is how two real
// sign/scale bugs were caught during development.

/**
 * Figure of merit: hover efficiency, ideal momentum-theory induced power
 * over shaft power. Meaningful only at (or near) zero forward speed.
 */
[[nodiscard]] double FigureOfMerit(const Result& result, double rho = SeaLevelDensity);

/**
 * Propulsive efficiency in forward flight: useful power out (thrust times
 * airspeed) over shaft power in. Zero in hover by definition, since a
 * stationary propeller does no useful work however much thrust it makes.
 */
[[nodiscard]] double PropulsiveEfficiency(const Result& result, double Vinf);

/**
 * Disk loading [N/m^2] -- thrust per unit disk area, the quantity that sets
 * induced velocity and therefore hover efficiency.
 */
[[nodiscard]] double DiskLoading(const Result& result);

// --------------------------------------------------------- Slipstream field
// Meant to be passed as Solver::Solve()'s externalField callback: returns the
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
/**
 * Callable slipstream velocity field, meant to be passed as
 * Solver::Solve()'s externalField callback: returns the induced velocity
 * perturbation (axial + swirl) at a global point P, to be added on top of
 * the vane surface's own freestream.
 */
struct SlipstreamField {
    Result BEMTResult;
    Vec3 HubCenter{0, 0, 0};
    Vec3 AxisDir{1, 0, 0}; ///< Unit vector, thrust direction.
    double DevelopmentLength = -1.0; ///< <0 => auto (1 prop diameter).

    /** Induced velocity perturbation at global point P. */
    [[nodiscard]] Vec3 operator()(const Vec3& P) const;
};

} // namespace Aeolion::BEMT
