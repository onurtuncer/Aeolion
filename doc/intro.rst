Introduction and Overview
==========================

Abstract
--------

Aeolion is a small aerodynamics toolkit built around a 3D vortex lattice
method (VLM). It solves for the lift distribution and induced drag of
lifting surfaces using horseshoe vortices, adds a component-buildup
viscous drag estimate on top, and includes a hover-safe blade-element
momentum theory (BEMT) solver for propellers. Aircraft and propeller
geometry is loaded from a versioned, SI-unit JSON contract rather than
from STL/mesh files, so the solver stays coupled to solver-native design
variables instead of a tessellated approximation.

Every module is checked against a closed-form or independently-known
result before being trusted: the VLM core tracks thin-wing lifting-line
theory across a range of aspect ratios, and the BEMT solver is checked
against hard physical bounds (Figure of Merit and propulsive efficiency
must both be :math:`\le 1`).

Modules
-------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Namespace
     - Responsibility
   * - ``Aeolion::Solver``
     - 3D VLM core: horseshoe vortices, LAPACK dense solve, sideslip,
       body rates, moments, and central-difference stability derivatives.
   * - ``Aeolion::PanelBuilder``
     - Builds the solver's lattice from a parsed handoff: CST camber
       surface, spanwise/chordwise discretization, deflected control
       surfaces.
   * - ``Aeolion::BEMT``
     - Propeller blade-element momentum theory. Hover-safe: solves for
       induced velocities directly rather than induction factors.
   * - ``Aeolion::Geometry``
     - Strict parser for the ``aeolion_geometry.json`` handoff contract,
       plus the CST section/planform/control-surface/mesh-topology data
       it describes.
   * - ``Aeolion::DragEstimate``
     - Viscous :math:`C_{D0}` component buildup: flat-plate skin
       friction, form factor, interference factor.
   * - ``Aeolion::Lattice``
     - Shared panel vocabulary (``Panel``) used by the builder, the
       solver, and the viewer alike.
   * - ``Aeolion::Math``
     - ``Vec3`` and shared numeric constants.

Known limitations (read before trusting results on a new geometry)
--------------------------------------------------------------------

- **Thin, flat lattice**: no thickness effects beyond camber-line
  averaging.
- ``CDi`` **is induced drag only.** VLM is inviscid; total drag needs the
  ``Aeolion::DragEstimate`` buildup added on top.
- **Linear, attached-flow method.** No stall, no separation physics.
- **BEMT is mid-fidelity** blade element momentum theory, not a
  substitute for measured prop data or a full rotating-lattice rotor VLM.

Nomenclature
------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Symbol
     - Description
   * - ``Vinf``
     - Freestream airspeed [m/s]
   * - ``alphaDeg``
     - Angle of attack [deg]
   * - ``rho``
     - Air density [kg/m^3]
   * - ``gamma``
     - Vortex ring/segment circulation strength
   * - ``CL``
     - Lift coefficient
   * - ``CDi``
     - Induced drag coefficient (VLM is inviscid -- this is not total drag)
   * - ``Cm``
     - Pitching moment coefficient about the reference point
       (positive: nose up)
   * - ``Vec3.x / .y / .z``
     - Cartesian vector components, body/world axes as documented per call site
   * - ``StabilityDerivatives``
     - Central-difference derivative table (``CL_q``, ``Cm_q``,
       ``Croll_p``, ``Cn_p``, ... and their non-dimensional ``_nd``
       counterparts)
