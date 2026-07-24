// Geometry/DuctGeometry.h
//
// The duct, as its own body -- disjoint from the fuselage (see
// BodyGeometry.h). Earlier schemas folded the duct into the main body by
// giving it an open tail (see BodyGeometry.h's header comment); schema 1.8.0
// splits it out because the duct is a physically separate structure, not a
// continuation of the fuselage skin.
//
// Unlike BodyGeometry, the duct is not a profile of stations: it is a single
// annular ring, stated as an axial chord and an inner/outer diameter, placed
// at one point in the contract's frame. There is no provision here for a
// duct whose cross-section varies circumferentially or axially.
#pragma once

#include "Aeolion/Math/Vec3.h"

namespace Aeolion::Geometry {

/** The duct, modeled as a single annular ring disjoint from the main body. */
struct DuctGeometry {
    bool IsStated = false;
    double Chord = 0.0;         ///< [m], axial extent of the duct ring.
    double InnerDiameter = 0.0; ///< [m]
    double OuterDiameter = 0.0; ///< [m], must exceed InnerDiameter.
    Math::Vec3 Center{0.0, 0.0, 0.0}; ///< [m], contract frame.
};

} // namespace Aeolion::Geometry
