// Geometry/PropulsionSpec.h
//
// Propeller blade geometry as handed off for a BEMT run. Radial positions are
// normalized (r/R); multiply by DiskRadius to get the metric radii that
// BEMT::PropGeometry wants. The innermost station's RadiusFraction is the hub
// cutout -- the blade is not defined inboard of it.
//
// Blade count arrived in schema 1.4.0 (n_blades); older documents omit it,
// leaving BladeCount zero for the consumer to supply as it did before.
#pragma once

#include <vector>

namespace Aeolion::Geometry {

/** One radial blade station as handed off (fraction-of-disk-radius form). */
struct BladeStationSpec {
    double RadiusFraction = 0.0; ///< r/R, 0..1.
    double Chord = 0.0;          ///< [m], may be exactly 0 at the tip.
    double TwistDeg = 0.0;       ///< Blade chordline angle relative to the rotor plane.
};

/** Propeller blade geometry as handed off for a BEMT run. */
struct PropulsionSpec {
    double DiskRadius = 0.0;   ///< [m]
    double ReferenceRpm = 0.0;
    int BladeCount = 0;        ///< 0 = not stated by this contract (schema < 1.4.0).
    std::vector<BladeStationSpec> BladeStations; ///< Ordered by increasing RadiusFraction.
};

} // namespace Aeolion::Geometry
