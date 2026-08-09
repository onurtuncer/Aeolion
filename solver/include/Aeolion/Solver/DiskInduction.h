// Solver/DiskInduction.h
//
// The velocity field a loaded rotor disk induces AROUND itself -- including
// UPSTREAM of it, which is the half that matters for an aft-mounted fan and
// the half the existing slipstream model does not carry.
//
// --- why this exists ---------------------------------------------------------
// Solver::SlipstreamField (ViscousCoupling.h) reconstructs a rotor's wake
// from annular momentum theory and returns zero for any point ahead of the
// disk plane. That is the right economy when the thing being influenced sits
// in the wake -- duct-jet vanes, a downstream tail -- and it is exactly wrong
// when the thing sits in front. On a tail-sitter with a pusher fan the wing's
// trailing edge is a fifth of a chord AHEAD of the duct, so a wake-only model
// reports precisely zero interaction, which is a property of the model rather
// than of the aircraft.
//
// A loaded disk does induce flow upstream: the axial induction rises from
// zero far ahead to v_i at the disk plane and reaches 2 v_i far downstream.
// The wing therefore sits inside an ACCELERATING stream, which is a
// favourable pressure gradient, which delays separation.
//
// --- the model ----------------------------------------------------------------
// A uniformly loaded actuator disk is exactly equivalent to a semi-infinite
// CYLINDRICAL VORTEX SHEET: azimuthal (ring) vorticity of constant strength
// gamma_t per unit length, shed from the disk edge and trailing downstream
// forever. That equivalence is the reason this file is worth having -- it
// turns "what does a disk do upstream" into a Biot-Savart integral over a
// known geometry, with a closed-form answer on the axis to check against.
//
// The sheet is discretized into vortex RINGS, each polygonized into straight
// filaments and evaluated with the ordinary Solver::SegmentVelocity kernel.
// That is the same choice the propeller wake makes (discretized helical legs
// rather than an analytic helix): it reuses a tested kernel, converges under
// refinement, and keeps the one special case -- the axis -- available as an
// exact check rather than as another approximation.
//
// Ring spacing GROWS geometrically downstream. The near field is what sets
// the induction at the disk and just upstream of it; the far tail contributes
// a slowly-varying remainder that does not deserve uniform resolution. This
// is what makes a fifty-radius cylinder affordable.
//
// --- annular disks --------------------------------------------------------------
// A ducted fan around a tail boom is an ANNULUS, not a disk. A uniformly
// loaded annulus carries no trailing vorticity between its radii (dGamma/dr
// is zero there); it sheds only at the two edges, the outer at +gamma_t and
// the inner at -gamma_t. So an annulus is two cylinders superposed, and
// nothing else in this file changes.
//
// --- what is deliberately NOT modelled -------------------------------------------
// SWIRL. The rotor's root vortex and the wake's axial vorticity live
// downstream of the disk, so they contribute nothing ahead of it, which is
// the region this file exists to serve. A consumer wanting the swirl in the
// jet should keep using SlipstreamField, which has it.
//
// NON-UNIFORM LOADING. The cylinder model assumes uniform disk loading. Real
// loading tapers to zero at both edges, which softens the sheet into a
// distributed band rather than a step. Superposing several cylinders at
// intermediate radii would represent that; a single pair is the honest first
// level and is what the closed-form check covers.

#pragma once

#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Solver/Solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numbers>
#include <vector>

namespace Aeolion::Solver {

// --- discretization defaults ---------------------------------------------------
/**
 * How far downstream the cylinder is carried, in disk radii. The sheet is
 * semi-infinite and this truncates it, so the number is a accuracy knob and
 * not a modelling choice.
 *
 * The truncation error is computable on the axis, which is how this default
 * was picked rather than guessed. A cylinder carried to L contributes
 *
 *     u(s) = (gamma_t/2) [ s/sqrt(s^2+R^2) - (s-L)/sqrt((s-L)^2+R^2) ]
 *
 * against the exact semi-infinite value, so one radius UPSTREAM (s = -R,
 * which is roughly where a wing sits) the error is 0.4% at L = 20R and
 * 0.06% at L = 50R. Fifty is cheap here because the spacing grows.
 */
inline constexpr double DefaultCylinderLengthRadii = 50.0;

/** Axial rings per radius AT THE DISK, before the geometric stretch. */
inline constexpr int DefaultRingsPerRadius = 8;

/**
 * Geometric growth of the axial ring spacing going downstream. Uniform
 * spacing over fifty radii would be thousands of rings for no gain: the far
 * tail's contribution varies slowly with position, while the first radius or
 * two sets essentially the whole upstream answer.
 */
inline constexpr double DefaultRingGrowth = 1.08;

/** Straight filaments per ring. A polygon's field converges quickly in this. */
inline constexpr int DefaultSegmentsPerRing = 24;

/**
 * Vortex core, as a fraction of the disk radius, applied to every filament.
 * The sheet is a mathematical singularity that a real slipstream boundary is
 * not, and a control point landing on it would otherwise return a very large
 * number. Only matters to a consumer evaluating ON the cylinder -- for a
 * wing upstream of the disk it never binds.
 */
inline constexpr double DefaultCylinderCoreFraction = 0.02;

/** Minimum radius below which a disk carries no sheet worth building. */
inline constexpr double MinDiskRadius = 1e-9;

// ------------------------------------------------------------------ the disk ---
/**
 * A uniformly loaded actuator disk, stated by its geometry and its induced
 * velocity at the disk plane.
 *
 * `Axis` points DOWNSTREAM -- the direction the disk drives the flow, which
 * for a propeller is the direction its jet leaves. Getting this backwards
 * turns the upstream acceleration into a deceleration, so it is stated
 * rather than inferred from a sign somewhere.
 */
struct ActuatorDisk {
    Vec3 Center{0, 0, 0};
    Vec3 Axis{1, 0, 0};        ///< Unit, pointing downstream.
    double Radius = 0.0;       ///< Outer radius [m].
    double HubRadius = 0.0;    ///< Inner radius of an annular disk [m]; zero for a full disk.
    /**
     * Axial induced velocity AT THE DISK PLANE [m/s], positive along Axis.
     * Momentum theory puts the fully developed wake at twice this, and the
     * sheet strength that produces both is gamma_t = 2 * this.
     */
    double InducedVelocity = 0.0;

    [[nodiscard]] bool Valid() const {
        return Radius > MinDiskRadius && HubRadius >= 0.0 && HubRadius < Radius;
    }
    /** Swept area of the (possibly annular) disk [m^2]. */
    [[nodiscard]] double Area() const {
        return std::numbers::pi * (Radius * Radius - HubRadius * HubRadius);
    }
    /** The cylindrical sheet strength, gamma_t = 2 v_i [m/s]. */
    [[nodiscard]] double SheetStrength() const { return Math::Two * InducedVelocity; }
};

/**
 * Momentum theory's induced velocity for a disk producing `thrust` while
 * flying at `axialSpeed` through it:
 *
 *     T = 2 rho A v_i (V + v_i)   =>   v_i = -V/2 + sqrt(V^2/4 + T/(2 rho A)).
 *
 * The positive root is taken, which is the propeller (accelerating) branch;
 * at V = 0 it reduces to the hover value sqrt(T / (2 rho A)). Returns zero
 * for a non-positive thrust rather than a complex or negative answer --
 * windmilling is a different momentum balance and is not this function's.
 */
[[nodiscard]] inline double InducedVelocityFromThrust(double thrust, double rho, double area,
                                                      double axialSpeed) {
    if (!(thrust > 0.0) || !(rho > 0.0) || !(area > 0.0)) return 0.0;
    const double half = Math::Half * std::max(axialSpeed, 0.0);
    return -half + std::sqrt(half * half + thrust / (Math::Two * rho * area));
}

// ------------------------------------------------- the closed-form axis solution
/**
 * Axial induced velocity ON the disk axis, exactly, at signed distance `s`
 * downstream of the disk plane (negative is upstream).
 *
 * For a single semi-infinite cylinder of radius R,
 *
 *     u(s) = (gamma_t/2) [ 1 + s / sqrt(s^2 + R^2) ],
 *
 * which is v_i at the disk, 2 v_i far downstream, and zero far upstream --
 * the three values momentum theory states, recovered from the vortex system
 * rather than assumed. An annulus superposes the outer cylinder and a
 * negative inner one, and on the axis the two cancel at every s: inside the
 * bore of an annular fan there is no jet, which is correct and is worth
 * having the model say on its own.
 *
 * This is the verification anchor for the discretized field below.
 */
[[nodiscard]] inline double DiskAxisInducedVelocity(const ActuatorDisk& disk, double s) {
    if (!disk.Valid()) return 0.0;
    const double half = Math::Half * disk.SheetStrength();
    const auto cylinder = [&](double radius) {
        return half * (1.0 + s / std::hypot(s, radius));
    };
    double u = cylinder(disk.Radius);
    if (disk.HubRadius > 0.0) u -= cylinder(disk.HubRadius);
    return u;
}

// ------------------------------------------------------- the discretized sheet ---
/** How the cylindrical sheet is discretized into ring filaments. */
struct VortexCylinderOptions {
    double LengthRadii = DefaultCylinderLengthRadii;
    int RingsPerRadius = DefaultRingsPerRadius;
    double Growth = DefaultRingGrowth;
    int SegmentsPerRing = DefaultSegmentsPerRing;
    double CoreFraction = DefaultCylinderCoreFraction;
};

/**
 * The sheet, built once as a list of straight filaments.
 *
 * Held as flat endpoint arrays with one circulation per filament so that
 * evaluating the field is a single loop over segments with no geometry work
 * per point -- an external-field callable is hit once per control point per
 * solve, and per-point trigonometry would dominate.
 */
struct VortexCylinderMesh {
    std::vector<Vec3> Start, End;
    std::vector<double> Gamma; ///< Circulation of each filament [m^2/s].
    double CoreRadius = 0.0;

    [[nodiscard]] std::size_t Count() const { return Gamma.size(); }
    [[nodiscard]] bool Empty() const { return Gamma.empty(); }
};

namespace Detail {

/** An orthonormal pair spanning the plane normal to `axis`. */
inline void DiskFrame(const Vec3& axis, Vec3& e1, Vec3& e2) {
    // Seed away from the axis so the cross product cannot collapse.
    const Vec3 seed = (std::fabs(axis.x) < 0.9) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
    e1 = Cross(axis, seed).Normalized();
    e2 = Cross(axis, e1).Normalized();
}

/**
 * Append one cylinder's rings. `sign` is +1 for the outer sheet and -1 for
 * an annulus's inner one.
 *
 * Ring circulation is gamma_t * dx: the sheet carries gamma_t per unit
 * length, so a ring standing for an axial slab of thickness dx carries that
 * slab's worth. The polygon is traversed in the +theta sense (e1 -> e2, the
 * right-hand rule about the downstream axis), which is the sense that drives
 * flow DOWNSTREAM inside the cylinder -- the solenoid convention. The axis
 * test pins it.
 */
inline void AppendCylinder(VortexCylinderMesh& mesh, const ActuatorDisk& disk, double radius,
                           double sign, const VortexCylinderOptions& options) {
    if (!(radius > MinDiskRadius)) return;

    Vec3 e1, e2;
    DiskFrame(disk.Axis, e1, e2);

    const int sides = std::max(options.SegmentsPerRing, 3);
    const double length = options.LengthRadii * disk.Radius;
    const double growth = std::max(options.Growth, 1.0);
    // Axial spacing keys to THIS cylinder's radius, not the disk's outer
    // one: a cylinder's field varies on the scale of its own radius, so an
    // annulus's inner sheet -- which can be a third the outer radius -- needs
    // proportionally finer rings. Keying both to the outer radius leaves the
    // inner sheet with about three rings across its own scale and costs an
    // order of magnitude in the annular axis field. The LENGTH stays keyed to
    // the outer radius, which only ever over-resolves the inner sheet's tail.
    double step = radius / std::max(options.RingsPerRadius, 1);

    // Precompute the polygon's corner directions once; every ring is the
    // same polygon translated along the axis and scaled by its radius.
    std::vector<Vec3> corner(static_cast<std::size_t>(sides) + 1);
    for (int k = 0; k <= sides; ++k) {
        const double theta = Math::Two * std::numbers::pi * k / sides;
        corner[static_cast<std::size_t>(k)] = (e1 * std::cos(theta) + e2 * std::sin(theta)) * radius;
    }

    for (double s = 0.0; s < length;) {
        const double thickness = std::min(step, length - s);
        const Vec3 center = disk.Center + disk.Axis * (s + Math::Half * thickness);
        const double gamma = sign * disk.SheetStrength() * thickness;

        for (int k = 0; k < sides; ++k) {
            mesh.Start.push_back(center + corner[static_cast<std::size_t>(k)]);
            mesh.End.push_back(center + corner[static_cast<std::size_t>(k) + 1]);
            mesh.Gamma.push_back(gamma);
        }
        s += thickness;
        step *= growth;
    }
}

} // namespace Detail

/** Discretize a disk's cylindrical vortex sheet (both edges, when annular). */
[[nodiscard]] inline VortexCylinderMesh BuildVortexCylinder(const ActuatorDisk& disk,
                                                            const VortexCylinderOptions& options = {}) {
    VortexCylinderMesh mesh;
    if (!disk.Valid() || !(std::fabs(disk.SheetStrength()) > 0.0)) return mesh;

    mesh.CoreRadius = std::max(options.CoreFraction, 0.0) * disk.Radius;
    Detail::AppendCylinder(mesh, disk, disk.Radius, +1.0, options);
    if (disk.HubRadius > 0.0) Detail::AppendCylinder(mesh, disk, disk.HubRadius, -1.0, options);
    return mesh;
}

/** Induced velocity of a built sheet at an arbitrary point. */
[[nodiscard]] inline Vec3 VortexCylinderVelocity(const VortexCylinderMesh& mesh, const Vec3& point) {
    Vec3 v(0, 0, 0);
    for (std::size_t i = 0; i < mesh.Count(); ++i)
        v = v + SegmentVelocity(point, mesh.Start[i], mesh.End[i], mesh.Gamma[i], mesh.CoreRadius);
    return v;
}

/**
 * The disk's induced field as an `externalField` callable for
 * Solver::Solve / SolveWithSystem.
 *
 * The sheet is built ONCE and captured, so the returned callable does no
 * geometry work per evaluation -- which matters, because the solver hits an
 * external field at every control point of every solve in an attitude sweep.
 *
 * This is a ONE-WAY coupling: the airframe sees the disk, the disk does not
 * see the airframe. Closing it needs the partitioned outer fixed point that
 * SolveRotorVaneCoupled already establishes for the rotor-vane problem,
 * because a static airframe and a rotating rotor cannot share one
 * FreestreamConditions.
 */
[[nodiscard]] inline std::function<Vec3(const Vec3&)> DiskInductionField(
    const ActuatorDisk& disk, const VortexCylinderOptions& options = {}) {
    return [mesh = BuildVortexCylinder(disk, options)](const Vec3& point) {
        return VortexCylinderVelocity(mesh, point);
    };
}

} // namespace Aeolion::Solver
