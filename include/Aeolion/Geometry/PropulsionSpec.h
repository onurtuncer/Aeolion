// Geometry/PropulsionSpec.h
//
// Propeller blade geometry as handed off for a BEMT run. Radial positions are
// normalized (r/R); multiply by DiskRadius to get the metric radii that
// BEMT::PropGeometry wants. The innermost station's RadiusFraction is the hub
// cutout -- the blade is not defined inboard of it.
//
// Blade count arrived in schema 1.4.0 (n_blades); older documents omit it,
// leaving BladeCount zero for the consumer to supply as it did before.
//
// Blade airfoil shape (BladeAirfoilSections) and RotationAxis arrived in
// schema 1.8.0; older documents leave both empty/zero. RotationAxis states
// the physical spin axis the blade stations are defined about -- it is not
// the same thing as BEMT::PropGeometry::RotationSign, which is a CW/CCW
// installation choice the contract does not state (see HandoffContract.h's
// ToPropGeometry bridge).
#pragma once

#include "Aeolion/Math/Vec3.h"

#include <vector>

namespace Aeolion::Geometry {

/** One radial blade station as handed off (fraction-of-disk-radius form). */
struct BladeStationSpec {
    double RadiusFraction = 0.0; ///< r/R, 0..1.
    double Chord = 0.0;          ///< [m], may be exactly 0 at the tip.
    double TwistDeg = 0.0;       ///< Blade chordline angle relative to the rotor plane.
};

/** CST section shape at a blade radial station (schema >= 1.8.0). */
struct BladeAirfoilSection {
    double RadiusFraction = 0.0; ///< r/R, 0..1.
    std::vector<double> CoefficientsUpper;
    std::vector<double> CoefficientsLower;
};

/** Propeller blade geometry as handed off for a BEMT run. */
struct PropulsionSpec {
    double DiskRadius = 0.0;   ///< [m]
    double ReferenceRpm = 0.0;
    int BladeCount = 0;        ///< 0 = not stated by this contract (schema < 1.4.0).
    std::vector<BladeStationSpec> BladeStations; ///< Ordered by increasing RadiusFraction.
    std::vector<BladeAirfoilSection> BladeAirfoilSections; ///< Ordered by increasing RadiusFraction; empty before schema 1.8.0.
    Math::Vec3 RotationAxis{0.0, 0.0, 0.0}; ///< Unit vector, contract frame; zero = not stated (schema < 1.8.0).
};

} // namespace Aeolion::Geometry
