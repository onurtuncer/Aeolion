// Geometry/Propeller.h
//
// The metric propeller a consumer works with: blade count, hub and tip
// radii, and the radial distribution of chord and twist. Pure data
// vocabulary -- it states what the propeller IS, and deliberately says
// nothing about how its aerodynamics get computed (the BEMT solver that
// used to consume this shape is now a separate project; whatever
// calculation method comes next poses itself on this same struct).
//
// Radii here are METRIC, unlike Geometry::PropulsionSpec, which states
// them as fractions of the disk radius because that is the form the blade
// is designed in. ToPropeller() in HandoffContract.h is that conversion,
// and it is the only place it should happen.
#pragma once

#include "Aeolion/Geometry/AirfoilSection.h"

#include <vector>

namespace Aeolion::Geometry {

/** One radial blade station in metric units. */
struct BladeStation {
    double r = 0.0;        ///< Radial position [m].
    double Chord = 0.0;    ///< Local chord [m].
    double TwistDeg = 0.0; ///< Local geometric pitch: blade chordline angle to the rotor plane.
};

/** The metric propeller: blade count, radii, and chord/twist stations. */
struct Propeller {
    int BladeCount = 2;
    double Radius = 0.15;    ///< Tip radius [m].
    double HubRadius = 0.02; ///< Hub radius [m].
    std::vector<BladeStation> Stations; ///< Sorted by increasing r, spanning [HubRadius, Radius].

    /**
     * Blade CST sections (schema >= 1.8.0), reusing the wing's
     * AirfoilSection vocabulary with Eta meaning r/R instead of semi-span
     * fraction -- CstSurface's CamberAt/CamberSlopeAt consume them
     * unchanged. Ordered by increasing Eta; empty means no section data
     * was handed off and the blades stay flat plates at the local twist.
     */
    std::vector<AirfoilSection> Sections;
};

} // namespace Aeolion::Geometry
