// BEMT/PropGeometry.h
//
// The propeller a BEMT run is posed on: blade count, hub and tip radii,
// and the radial distribution of chord and twist.
//
// Split out of BEMT.h so that code which only needs to DESCRIBE a
// propeller -- the geometry contract bridge, a CLI, a trim routine -- does
// not have to pull in the whole blade-element solver to do it.
//
// Radii here are METRIC, unlike Geometry::PropulsionSpec, which states them
// as fractions of the disk radius. That conversion is the bridge's job (see
// Geometry::ToPropGeometry).
#pragma once

#include <vector>

namespace Aeolion::BEMT {

struct BladeStation {
    double r = 0.0;        // radial position [m]
    double Chord = 0.0;    // local chord [m]
    double TwistDeg = 0.0; // local geometric pitch: blade chordline angle to the rotor plane
};

struct PropGeometry {
    int NBlades = 2;
    double Radius = 0.15;    // tip radius [m]
    double HubRadius = 0.02; // hub radius [m]
    std::vector<BladeStation> Stations; // sorted by increasing r, spanning [HubRadius, Radius]

    // Sense of rotation about the prop axis, +1 or -1 by the right-hand
    // rule about the thrust direction. Not part of the geometry contract --
    // it is an installation choice, so the consumer states it.
    int RotationSign = 1;
};

} // namespace Aeolion::BEMT
