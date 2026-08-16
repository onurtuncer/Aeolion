.. ------------------------------------------------------------------------------
.. Project: Aeolion
.. Copyright (c) 2025-2026, Onur Tuncer, PhD
..
.. SPDX-License-Identifier: MIT
.. License-Filename: LICENSE
.. ------------------------------------------------------------------------------

.. _api:

Aeolion C++ API Reference
==========================

This page documents the public C++ API extracted from source headers via
Doxygen and rendered by Breathe.

.. contents:: Modules
   :depth: 1
   :local:

Solver
------

The 3D vortex lattice method core: horseshoe vortices, LAPACK dense
solve, sideslip, body rates, moments, and central-difference stability
derivatives. Every other module builds on top of this one, which is why
it is kept as its own top-level component rather than folded into
``include/``.

Key types:

* :cpp:class:`Aeolion::Solver::WingParams` — parametric single-wing planform.
* :cpp:class:`Aeolion::Solver::FreestreamConditions` — flight condition
  (``Vinf``, ``alphaDeg``, sideslip, body rates, ``rho``, reference point).
* :cpp:class:`Aeolion::Solver::ReferenceGeometry` — coefficient
  normalization constants.
* :cpp:class:`Aeolion::Solver::SolveResult` — full solve result (``CL``,
  ``CDi``, ``Cm``, per-surface forces).
* :cpp:class:`Aeolion::Solver::StationResult` — per-spanwise-station output.
* :cpp:class:`Aeolion::Solver::StabilityDerivatives` — central-difference
  derivative table.
* :cpp:class:`Aeolion::Solver::FlowField` — a solved system's velocity
  field, evaluable at any point after the solve.

Stagnation and attachment lines (see :doc:`theory`):

* :cpp:class:`Aeolion::Solver::SurfaceGrid` — a source-panelled surface as
  a structured (station, sector) patch with its solved skin flow.
* :cpp:class:`Aeolion::Solver::CriticalPoint` — a located and classified
  zero of the surface flow: attachment/separation node, saddle, or focus.
* :cpp:class:`Aeolion::Solver::SurfaceStreamline` — a traced surface
  streamline carrying :math:`U_e(s)` and the spreading metric :math:`h(s)`.
* :cpp:class:`Aeolion::Solver::SectionSolution` — Hess-Smith solve of one
  thick section: surface velocity, the stagnation point, and its strain rate.
* :cpp:class:`Aeolion::Solver::SurfaceRun` — one surface's boundary-layer
  run, measured from the stagnation point.
* :cpp:class:`Aeolion::Solver::AttachmentStation` — the wing attachment
  line at one strip, in leading-edge-normal coordinates, with Poll's
  attachment-line Reynolds number.
* :cpp:class:`Aeolion::Solver::SurfaceMarch` — one surface run marched
  from its attachment point to separation: bubble location, transition,
  and the turbulent separation point.
* :cpp:class:`Aeolion::Solver::SeparationSurvey` — the separation picture
  across a lifting surface at one flight condition.
* :cpp:class:`Aeolion::Solver::TrefftzResult` — far-field induced drag and
  span efficiency, from the wake trace rather than the near-field forces;
  the trustworthy ``CDi`` on a coupled configuration.
* :cpp:class:`Aeolion::Solver::BodyAxisCoefficients` --- body-axis (FRD)
  force and moment coefficients, and
  :cpp:class:`Aeolion::Solver::BodyAxisRateDerivatives` --- the
  reduced-rate stability derivatives in that frame. The conversion from
  the solver's working frame is a 180-degree rotation about :math:`y`;
  note that a *rate derivative* does not follow the wrench rule, since
  the flip applies to both the response and the rate and cancels for the
  roll/yaw pairs (see the header, and ``TestBodyAxes``).
* :cpp:class:`Aeolion::Solver::StripSeparationTable` --- the per-strip
  separation point against local incidence that the anchored post-stall
  model is driven by, built once from an inviscid alpha sweep.
* :cpp:class:`Aeolion::Solver::PostStallSectionModel` — the anchored
  post-stall section model: Kirchhoff attenuation from the computed
  separation point, Viterna-Corrigan deep stall with AR-aware
  :math:`C_{d,\max}`, Rayleigh centre of pressure for the section
  :math:`c_m` (see :doc:`theory`).
* :cpp:class:`Aeolion::Solver::ViternaConstants` — the Viterna
  extension's matching constants, fixed by continuity at the emergent
  stall junction.
* :cpp:class:`Aeolion::Solver::DiscreteVortexResult` — the unsteady
  discrete-vortex cross-check's converged mean/RMS loads and shedding
  counts (see :doc:`theory`).
* :cpp:class:`Aeolion::Solver::ParticleWakeResult` — the
  three-dimensional particle-wake cross-check's mean/RMS loads, with
  the circulation-lift diagnostic that separates bound-solve defects
  from impulse-accounting ones (see :doc:`theory`).
* :cpp:class:`Aeolion::Solver::ActuatorDisk` — a uniformly loaded rotor
  disk, stated by geometry and disk-plane induced velocity.
* :cpp:class:`Aeolion::Solver::VortexCylinderMesh` — its semi-infinite
  cylindrical vortex sheet, discretized into ring filaments; the field it
  induces is valid upstream of the disk as well as in the wake.

.. doxygennamespace:: Aeolion::Solver
   :content-only:

PanelBuilder
------------

Builds the solver's lattice from a parsed geometry handoff: CST camber
surface, spanwise and chordwise discretization, deflected control
surfaces, with breakpoints on every surface edge. A compiled STATIC
library (``aeolion_panelbuilder``) rather than headers alone; consumers
link it explicitly in addition to ``aeolion``.

.. doxygennamespace:: Aeolion::PanelBuilder
   :content-only:

BEMT (external project)
-----------------------

Propeller blade-element momentum theory now lives in its own repository,
`onurtuncer/BEMT <https://github.com/onurtuncer/BEMT>`_ -- it is a
momentum method, not a panel method, so it was excised from this toolkit,
which carries no dependency on it. The handoff contract's
``propulsion_bemt`` block (``Geometry::PropulsionSpec``) remains part of
the schema as propulsion vocabulary for whatever calculation method
consumes it next. See that repository for the solver's API and theory.

Geometry
--------

Strict parser for the ``aeolion_geometry.json`` handoff contract.

* :cpp:class:`Aeolion::Geometry::HandoffContract` — top-level parser and
  contract invariants.
* :cpp:class:`Aeolion::Geometry::CstSurface` — CST evaluation: camber
  mean line, its slope, and camber-line arc length.
* :cpp:class:`Aeolion::Geometry::PlanformStation` — one spanwise planform
  station.
* :cpp:class:`Aeolion::Geometry::AirfoilSection` — CST section shape at a
  station.
* :cpp:class:`Aeolion::Geometry::ControlSurface` — hinged surface and
  which body it binds to.
* :cpp:class:`Aeolion::Geometry::MeshTopology` — requested lattice
  discretization.
* :cpp:class:`Aeolion::Geometry::PropulsionSpec` — propeller blade
  geometry for a BEMT run.
* :cpp:class:`Aeolion::Geometry::SectionContour` — the THICK closed
  contour of a CST section, with its leading-edge radius; what a
  stagnation-point calculation needs and a camber line cannot supply.

.. doxygennamespace:: Aeolion::Geometry
   :content-only:

DragEstimate
------------

Viscous :math:`C_{D0}` component buildup: flat-plate skin friction, form
factor, and interference factor.

.. doxygennamespace:: Aeolion::DragEstimate
   :content-only:

Lattice
-------

Panel vocabulary shared between the builder, the solver, and the viewer
(``Panel`` — one horseshoe-vortex panel). Re-exported into
``Aeolion::Solver`` for source compatibility.

.. doxygennamespace:: Aeolion::Lattice
   :content-only:

Math
----

``Vec3`` (3D double vector with dot/cross/axis-rotation helpers) and
shared numeric constants used throughout the toolkit.

.. doxygennamespace:: Aeolion::Math
   :content-only:
