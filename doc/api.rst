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

.. doxygennamespace:: Aeolion::Solver
   :content-only:

PanelBuilder
------------

Builds the solver's lattice from a parsed geometry handoff: CST camber
surface, spanwise and chordwise discretization, deflected control
surfaces, with breakpoints on every surface edge. The one component with
compiled sources rather than headers alone.

.. doxygennamespace:: Aeolion::PanelBuilder
   :content-only:

BEMT
----

Propeller blade-element momentum theory. Hover-safe: solves for induced
velocities directly rather than induction factors, and exposes the
resulting slipstream field for downstream control vanes.

.. doxygennamespace:: Aeolion::BEMT
   :content-only:

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
