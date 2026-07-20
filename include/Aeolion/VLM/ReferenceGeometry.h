// VLM/ReferenceGeometry.h
//
// Normalization constants for force/moment coefficients. A generic panel
// list carries no planform metadata of its own, so these must be supplied
// by the caller (the WingParams convenience overload fills them in
// automatically for a single parametric wing).
#pragma once

namespace Aeolion::VLM {

struct ReferenceGeometry {
    double Area = 0.0;   // S, wing reference area [m^2]
    double Chord = 0.0;  // reference chord (mean aerodynamic chord) -- normalizes Cm [m]
    double Span = 0.0;   // reference span -- normalizes Croll, Cn [m]
};

} // namespace Aeolion::VLM
