// VLM/Panel.h
//
// One horseshoe-vortex panel: the atomic geometric unit of the lattice.
#pragma once
#include <string>
#include "Aeolion/Math/Vec3.h"

namespace Aeolion::VLM {

struct Panel {
    Vec3 A, B;              // bound vortex endpoints (quarter-chord), A.y < B.y
    Vec3 ControlPoint;      // three-quarter-chord, mid-span
    Vec3 Normal;            // unit outward normal at control point (includes twist+dihedral)
    Vec3 TrailDirA, TrailDirB; // unit direction of trailing legs from A and from B (downstream)
    double Area = 0.0;      // panel planform area (for Cl-per-station bookkeeping)
    double SpanwiseWidth = 0.0;
    std::string Surface;    // which lifting surface this panel belongs to (e.g. "wing", "htail")
};

} // namespace Aeolion::VLM
