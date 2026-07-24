---
title: 'Aeolion: A 3D vortex-lattice, source-panel, and BEMT toolkit for coupled wing-body-propeller aerodynamic analysis'
tags:
  - C++
  - aerodynamics
  - vortex lattice method
  - blade element momentum theory
  - aircraft design
authors:
  - name: Onur Tuncer
    orcid: ''
    affiliation: 1
affiliations:
  - name: '' # institution
    index: 1
date: '' # %e %B %Y
bibliography: paper.bib
---

## Summary

Aeolion is a C++23 aerodynamics toolkit combining a 3D vortex lattice
method (VLM) for lifting surfaces, source panels for bodies of
revolution, and a blade-element-momentum-theory (BEMT) solver for
propellers, coupled together into a single blocked linear system so
that wing, fuselage, and propeller effects can be analyzed together
rather than in isolation. Geometry is loaded from a versioned JSON
handoff contract rather than an STL/mesh file, which keeps the
description parametric (planform stations, CST-defined camber, control
surface deflections, blade geometry) and traceable back to the CAD
source that generated it.

## Statement of need

<!-- Who needs this, what problem does it solve, how does it relate to
other tools? Fill in with reference to conceptual-design-stage
aero tools that either (a) only handle an isolated lifting surface,
or (b) require full CFD meshing for wing-body-propeller interaction. -->

## State of the field

<!-- Compare against other open VLM/BEMT tools (e.g. AVL, OpenVSP's
VSPAERO, XROTOR/XFOIL-adjacent tools, MachUp/MachUpX) and state the
build-vs-contribute justification: JSON-contract geometry input,
permeable-panel propeller-to-fuselage coupling, wing circulation
carried through the fuselage. -->

## Software design

Aeolion is organized as one namespace per header folder
(`Aeolion::Math`, `Aeolion::Lattice`, `Aeolion::BEMT`,
`Aeolion::Geometry`, `Aeolion::DragEstimate`), with larger components
(`Solver`, `PanelBuilder`, `viewer`) promoted to top-level directories
with their own build targets. The wing-body-propeller coupling is
implemented as a single blocked linear system solved once per
geometry via LAPACK's `dgetrf`/`dgetrs`, with permeable source panels
carrying the propeller's slipstream efflux onto the fuselage base
without imposing spurious pressure there.

<!-- Expand: trade-offs behind the JSON-contract geometry boundary
(vs. STL/mesh), why discretization is fixed per solve (CppAD future
direction), why LAPACK's Fortran ABI is called directly rather than
through LAPACKE. -->

## Research impact statement

<!-- Needs concrete evidence: has this been used for the VBAT UAV
program, any resulting internal reports/design decisions, citations,
or other adopters? Aspirational statements ("could be used for...")
are explicitly insufficient for JOSS review. -->

## AI usage disclosure

<!-- State plainly which parts (if any) of the software, docs, or this
paper were produced with AI assistance, and how correctness was
verified (e.g. the closed-form/physical-bound validation approach
documented in the repository README). -->

## Acknowledgements

<!-- Funding, collaborators, VBAT program acknowledgements, etc. -->

## References
