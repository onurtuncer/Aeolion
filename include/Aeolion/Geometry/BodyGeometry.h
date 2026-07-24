// Geometry/BodyGeometry.h
//
// The fuselage/duct body, as an axisymmetric body of revolution: a profile
// of (x, radius) stations that the surface is swept from about the body
// axis. The contract carries no circumferential variation, so the body is
// axisymmetric by construction.
//
// --- reference frame ------------------------------------------------------
// x is in the contract's own frame (aetherion_body_frd: x FORWARD), so the
// stations run from the nose at x = 0 to the tail at x = -Length. That is
// the opposite sense to the solver's x-aft convention; the consumer
// converts explicitly at ingest, exactly as ControlSurface.h describes for
// hinge axes. Nothing here is auto-rotated.
//
// --- what the profile is allowed to be ------------------------------------
// The nose may be sharp (radius 0) or blunt, and the TAIL MAY BE OPEN --
// a nonzero final radius is not a malformed contract, it is a truncated
// base. On this airframe the open tail is the duct exit that the duct-jet
// vanes sit in (see ControlSurface.h), so refusing it would refuse the
// actual vehicle. Any consumer that panels this body has to decide what to
// do about that base; the contract's job is only to state it faithfully.
#pragma once

#include <cstddef>
#include <vector>

namespace Aeolion::Geometry {

/** One axial station of the body-of-revolution profile. */
struct BodyStation {
    double x = 0.0;      ///< [m], contract frame (x forward, so aft is negative).
    double Radius = 0.0; ///< [m], >= 0; zero means the profile closes to a point.
};

/** Fuselage/duct body, as an axisymmetric body of revolution. */
struct BodyGeometry {
    double Length = 0.0;              ///< [m], nose-to-tail.
    std::vector<BodyStation> Stations; ///< Ordered nose -> tail (x strictly decreasing).

    [[nodiscard]] bool IsPresent() const { return !Stations.empty(); }
};

/**
 * The body's radius law: its radius at an arbitrary axial position, in the
 * CONTRACT's frame (x forward, so the body occupies [-Length, 0]). Linear
 * between stations, which is the interpolation the contract itself
 * specifies -- the same rule PlanformStation.h states for the planform.
 *
 * Returns 0 outside the body's extent, so "is this point inside the body?"
 * reduces to comparing a radius against this without a separate bounds
 * test. That question is what trimming a wing lattice at the fuselage
 * surface is made of, which is why this lives here rather than being
 * re-derived by each consumer.
 */
[[nodiscard]] inline double RadiusAt(const BodyGeometry& body, double contractX) {
    const auto& stations = body.Stations;
    if (stations.empty()) return 0.0;
    // Stations run nose to tail, so x DECREASES along the list.
    if (contractX > stations.front().x || contractX < stations.back().x) return 0.0;

    for (std::size_t i = 0; i + 1 < stations.size(); ++i) {
        const double near = stations[i].x;
        const double far = stations[i + 1].x;
        if (contractX <= near && contractX >= far) {
            const double span = near - far;
            const double t = (span > 0.0) ? (near - contractX) / span : 0.0;
            return stations[i].Radius + (stations[i + 1].Radius - stations[i].Radius) * t;
        }
    }
    return 0.0;
}

} // namespace Aeolion::Geometry
