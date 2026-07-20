// Lattice/SourcePanel.h
//
// One constant-strength quadrilateral SOURCE panel: the atomic unit of a
// non-lifting body's surface, the way Panel is the atomic unit of a lifting
// surface's lattice.
//
// --- why sources and not horseshoes ---------------------------------------
// A fuselage is a closed non-lifting volume. Its boundary condition is
// "no flow through the surface", which a source distribution enforces
// without introducing bound circulation or shedding a wake. Panelling a
// body with horseshoe vortices instead would have it generate circulation
// and shed lift it physically does not.
//
// Sources produce no net force on a closed body (d'Alembert), which is the
// correct potential-flow answer. What they DO produce is the surface
// pressure distribution, and hence the body's pitching/yawing moment at
// incidence -- the Munk moment, which is the dominant destabilizing
// contribution a fuselage makes to an airframe's stability. They also
// produce the upwash the wing sees, which is the other effect that matters.
//
// --- the base ---------------------------------------------------------------
// PrescribedNormalVelocity exists because this airframe's body does not
// close: its base is a duct exit (see Geometry/BodyGeometry.h). A solid
// wall panel wants zero normal velocity; a base panel passing efflux wants
// the velocity the jet actually delivers there. Rather than special-case
// the base in the solver, every panel carries the normal velocity its
// boundary condition should enforce, and a solid wall simply leaves it at
// zero.

#pragma once

#include "Aeolion/Math/Vec3.h"

#include <array>
#include <string>

namespace Aeolion::Lattice {

using Math::Vec3;

struct SourcePanel {
    // Corners in order around the panel; the winding must be consistent with
    // Normal (right-hand rule), since the influence kernel integrates over
    // the edges in this order.
    std::array<Vec3, 4> Corners;

    Vec3 ControlPoint;    // centroid, where tangency is enforced
    Vec3 Normal;          // unit outward normal (out of the body)
    double Area = 0.0;    // true panel area

    // Normal velocity the boundary condition should enforce at ControlPoint,
    // positive along Normal. Zero is a solid wall; nonzero is transpiration,
    // which is how the open base carries the propeller efflux.
    double PrescribedNormalVelocity = 0.0;

    std::string Surface;  // which body this belongs to (e.g. "fuselage")

    // Axial station index along the body, so a chordwise/axial pressure
    // distribution can be reported without re-deriving the topology.
    int StationIndex = -1;
};

} // namespace Aeolion::Lattice
