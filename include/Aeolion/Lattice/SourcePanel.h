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

/**
 * One constant-strength quadrilateral SOURCE panel: the atomic unit of a
 * non-lifting body's surface, the way Panel is the atomic unit of a lifting
 * surface's lattice.
 */
struct SourcePanel {
    /**
     * Corners in order around the panel; the winding must be consistent with
     * Normal (right-hand rule), since the influence kernel integrates over
     * the edges in this order.
     */
    std::array<Vec3, 4> Corners;

    Vec3 ControlPoint;    ///< Centroid, where tangency is enforced.
    Vec3 Normal;          ///< Unit outward normal (out of the body).
    double Area = 0.0;    ///< True panel area.

    /**
     * Normal velocity the boundary condition should enforce at ControlPoint,
     * positive along Normal. Zero is a solid wall; nonzero is transpiration,
     * which is how the open base carries the propeller efflux.
     */
    double PrescribedNormalVelocity = 0.0;

    /**
     * Whether this panel is a real wall or just a boundary the flow passes
     * through. A permeable panel still imposes its condition -- that is how
     * the efflux shapes the field around the rest of the body -- but it
     * contributes NO pressure force, because there is no surface there for
     * pressure to act on.
     *
     * Getting this wrong is subtle rather than loud: integrating Cp over a
     * permeable face silently converts a fast jet into base suction, which
     * reads as drag. A jet's thrust is momentum flux through the exit plane,
     * which is not a pressure force on the body at all.
     */
    bool Permeable = false;

    std::string Surface;  ///< Which body this belongs to (e.g. "fuselage").

    /**
     * Axial station index along the body, so a chordwise/axial pressure
     * distribution can be reported without re-deriving the topology.
     */
    int StationIndex = -1;

    /**
     * Azimuthal index around the body at this station -- the SECOND surface
     * coordinate, of which StationIndex is the first. Together they make a
     * swept body of revolution's panelling a structured (station x sector)
     * grid.
     *
     * Reporting a pressure distribution needs only StationIndex, which is
     * why it came first. Anything that has to move ALONG the surface needs
     * both: a surface streamline crosses panels, so its integrator has to
     * interpolate the tangential velocity between control points, and that
     * requires knowing which panels neighbour which. Recovering the
     * azimuth from atan2(z, y) of a control point works only for a body
     * whose axis happens to be +x and would silently produce garbage for
     * one that is not, so the builder states the index instead.
     *
     * Negative means the panelling is unstructured (or its topology was
     * never recorded), and a surface-flow analysis must decline rather
     * than guess -- see Solver/SurfaceFlow.h.
     */
    int SectorIndex = -1;
};

} // namespace Aeolion::Lattice
