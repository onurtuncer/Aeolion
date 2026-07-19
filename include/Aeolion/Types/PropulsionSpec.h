// Types/PropulsionSpec.h
//
// Propeller blade geometry as handed off for a BEMT run. Radial positions are
// normalized (r/R); multiply by DiskRadius to get the metric radii that
// BEMT::PropGeometry wants. The innermost station's RadiusFraction is the hub
// cutout -- the blade is not defined inboard of it.
//
// TODO Blade count is not part of the schema; the consumer supplies it.
#pragma once

#include <vector>

namespace Aeolion::Geometry {

struct BladeStationSpec {
    double RadiusFraction = 0.0; // r/R, 0..1
    double Chord = 0.0;          // [m], may be exactly 0 at the tip
    double TwistDeg = 0.0;       // blade chordline angle relative to the rotor plane
};

struct PropulsionSpec {
    double DiskRadius = 0.0;   // [m]
    double ReferenceRpm = 0.0;
    std::vector<BladeStationSpec> BladeStations; // ordered by increasing RadiusFraction
};

} // namespace Aeolion::Geometry
