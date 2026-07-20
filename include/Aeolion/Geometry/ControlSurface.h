// Geometry/ControlSurface.h
//
// A hinged control surface. What it is hinged TO is given by Binding, and
// that choice reinterprets the eta band -- see below.
//
// NOTE: names are NOT unique. A four-vane duct emits four entries all named
// "vane", distinguished only by their hinge axis, so index these by position
// and never key them by name.
//
// --- reference frame ------------------------------------------------------
// HingeAxis is expressed in the contract's own frame. For
// "aetherion_body_frd" that is x-forward / y-right / z-down, which is NOT the
// solver's convention (x aft, y right, z up). Per ADR-0016 the contract states its
// frame once and the consumer converts explicitly at ingest, so nothing here
// is auto-rotated.
//
// The conversion is a 180-degree rotation about y:
//     x_vlm = -x_frd,  y_vlm = y_frd,  z_vlm = -z_frd
// A proper rotation, so handedness is preserved and the axes map cleanly --
// [0,1,0] is invariant, [0,0,+-1] flips sign.
//
// TRAP: deflection sign conventions survive that rotation only if you rotate
// the AXIS and keep the right-hand rule about it. Re-deriving a convention in
// the new frame from an English description ("positive = trailing edge down")
// will silently invert half the controls. Rotate the vector; keep the rule.
#pragma once

#include "Aeolion/Math/Vec3.h"
#include <string>

namespace Aeolion::Geometry {

// Which body a control surface is hinged to. This determines whether the
// surface belongs in the solver's lattice at all, and what its eta band measures.
enum class ControlSurfaceBinding {
    // Hinged to the wing planform. Eta is the semi-span fraction, directly
    // comparable to PlanformStation::Eta.
    Wing,
    // An all-moving plate inside the duct jet at the duct exit -- NOT a
    // lifting surface, and it must never be bound to the wing lattice. Eta is
    // the duct-exit RADIUS fraction, a different quantity that happens to
    // share the name. These are characterized separately (vane CFD ->
    // DAVE-ML); the geometry rides along in the contract so the jet can be
    // modelled directly if the solver chooses to.
    DuctJet,
};

// Travel limits, in degrees about the hinge axis. The hard limits are the
// mechanical stops. The soft limit, where present, is the smaller angle a
// controller should respect in normal operation -- it exists because a vane
// can be driven to its stop but should not be, so a trim/allocation routine
// wanting a usable range should read SoftLimitDeg and not MaxDeg.
struct DeflectionLimits {
    double MinDeg = 0.0;
    double MaxDeg = 0.0;
    bool HasSoftLimit = false;
    double SoftLimitDeg = 0.0; // symmetric: usable range is +/- this
};

struct ControlSurface {
    std::string Name;
    double ChordFraction = 0.0; // fraction of local chord aft of the hinge, 0..1
    double EtaStart = 0.0;      // interpretation depends on Binding
    double EtaEnd = 0.0;
    Math::Vec3 HingeAxis{0.0, 1.0, 0.0}; // unit vector, contract reference frame
    ControlSurfaceBinding Binding = ControlSurfaceBinding::Wing;
    DeflectionLimits Limits;
};

} // namespace Aeolion::Geometry
