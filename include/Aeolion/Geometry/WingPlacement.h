// Geometry/WingPlacement.h
//
// Where the wing sits relative to the body. Without this the planform and
// the body each describe themselves in their own terms and nothing ties
// them together -- both would start at the origin, putting the wing at the
// fuselage nose.
//
// The anchor is the ROOT LEADING EDGE, stated in the contract's own frame
// (aetherion_body_frd: x forward, y right, z down), so a consumer converts
// it exactly as it converts hinge axes and body stations. The leading edge
// is the anchor rather than the quarter chord because it is a point on the
// physical surface, independent of any chordwise convention the solver
// happens to use internally.
//
// Absent before schema 1.5.0. IsStated distinguishes "the wing is at the
// origin" from "nobody said", which matters: the first is a placement and
// the second is a gap that silently buries the wing in the nose.
#pragma once

#include "Aeolion/Math/Vec3.h"

namespace Aeolion::Geometry {

/** Where the wing sits relative to the body, anchored at the root leading edge. */
struct WingPlacement {
    bool IsStated = false;
    Math::Vec3 RootLeadingEdge{0.0, 0.0, 0.0}; ///< [m], contract frame.
};

} // namespace Aeolion::Geometry
