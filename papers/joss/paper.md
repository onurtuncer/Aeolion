---
title: 'Aeolion: coupled vortex-lattice and source-panel analysis of wing-body-propeller aircraft, including ducted fans with jet vanes'
tags:
  - C++
  - aerodynamics
  - vortex lattice method
  - panel method
  - propeller
  - ducted fan
  - aircraft design
authors:
  - name: Onur Tuncer
    orcid: '' # TODO
    affiliation: 1
affiliations:
  - name: '' # TODO: institution
    index: 1
date: '' # %e %B %Y
bibliography: paper.bib
---

## Summary

Aeolion is a C++23 aerodynamics toolkit for the conceptual and
preliminary design of propeller-driven aircraft. Lifting surfaces are
modeled with a 3D vortex lattice method (VLM) and bodies of revolution
with constant-strength source panels [@katzPlotkin2001], assembled
into a single blocked linear system solved with LAPACK
[@lapack1999], so that wing, fuselage, and control-surface effects
are analyzed together rather than in isolation: the wing's bound
circulation is carried through the fuselage, and permeable base
panels let a propeller's efflux cross the body without imposing
spurious pressure.

Propellers are solved with the same vortex-lattice machinery posed in
the blade-fixed rotating frame: blades are meshed on their CST camber
surfaces [@kulfan2008] wrapped on radius cylinders, trailing legs
follow the local relative wind (which keeps the hover limit
well-posed), and thrust and shaft torque are read off as the induced
drag and rolling moment of the rotating solve. On top of the inviscid
lattice, a leveled viscous-inviscid coupling adds sectional lift and
profile-drag feedback per radial strip, up to a
transpiration-coupled integral boundary-layer section solver
(Thwaites [@thwaites1949], Head [@head1958], Squire-Young
[@squireYoung1938]) whose displacement effect enters the lattice as
transpiration boundary conditions; the outer fixed point is
stabilized with Anderson acceleration [@anderson1965].

The same solve extends to ducted propulsors: a shroud is paneled into
the rotor solution as a two-way interaction (the duct's boundary
condition sees the propwash, the blades see the bore's constriction),
and all-moving vanes at the duct exit are meshed as lifting surfaces
in the inertial frame, reading the rotor through a slipstream
reconstructed from the converged radial loading by annular momentum
theory. A block Gauss-Seidel alternation couples rotor and vanes
two-way, with the recovered swirl constrained by the jet's
angular-momentum budget — the vanes cannot extract more torque than
the shaft delivers. The result is the complete propulsive wrench of a
ducted fan with jet vanes, including control moments at zero airspeed,
where the propwash is the vanes' only dynamic pressure.

Geometry is loaded from a versioned JSON handoff contract rather than
an STL/mesh file, keeping the description parametric (planform
stations, CST camber, control-surface hinges, blade sections) and
traceable to the CAD source that generated it. An OpenGL viewer
renders the solved lattices with interactive control-surface and vane
deflection, and scripted multi-angle screenshots for documentation.

## Statement of need

Conceptual-design aerodynamics tools mostly sit at two extremes:
classical VLM/panel codes that treat an isolated lifting surface (or
uncoupled components), and full CFD, which requires volume meshing and
hours per configuration. The gap is widest for the vehicle class
Aeolion targets — VTOL aircraft whose propulsor and airframe cannot be
analyzed separately. A ducted fan with jet vanes at hover is the
sharpest case: at zero airspeed a conventional VLM has no freestream
dynamic pressure to put on a control surface, so control authority,
swirl recovery, and the thrust cost of vectoring are simply outside
such tools' model, while CFD is too slow for the design loop. Aeolion
resolves the whole interaction — rotor, duct, vanes, and airframe — in
seconds per operating point, at the fidelity level appropriate to
sizing and control-layout decisions.

## State of the field

AVL [@avl] and OpenVSP's VSPAERO [@openvsp2022] are the standard
open VLM tools for airframes; propeller effects enter, if at all, as
prescribed actuator-disk models rather than solved rotor aerodynamics.
XROTOR [@xrotor] and blade-element-momentum codes solve the isolated
rotor but carry no airframe. Aeolion's contribution is not a new
singularity element — its horseshoe vortices and source panels are
textbook [@katzPlotkin2001] — but the coupling: one blocked system
for wing and body, a rotating-frame lattice for the rotor sharing the
same kernels, a two-way ducted-fan solve, and a
momentum-budget-constrained rotor-vane alternation, none of which the
established tools provide. The JSON-contract geometry boundary
(instead of mesh files) and the header-only C++23 core are secondary
design points aimed at embedding the solver in design pipelines.

## Software design

The library is organized as one namespace per header folder
(`Aeolion::Math`, `Aeolion::Lattice`, `Aeolion::Geometry`,
`Aeolion::Solver`, `Aeolion::DragEstimate`), with compiled components
(`PanelBuilder`, the viewer) promoted to top-level directories with
their own build targets and test suites. The solver core is
header-only; the dense blocked system is factorized once per geometry
via LAPACK's `dgetrf`/`dgetrs` and reused across right-hand sides.

Two design decisions carry most of the physics architecture. First,
the frame split: blades (and the duct, a body of revolution about the
rotation axis) solve in the rotating frame, exit vanes solve static,
and the two communicate through a reconstructed time-mean slipstream
rather than instantaneous induction — blade-passing effects are
deliberately outside the model. Second, the leveled viscous coupling:
the 2-D section model is an explicit interface seam, so the analytic
polar and the integral boundary-layer solver are interchangeable
behind the same lattice coupling loop.

Every module was validated against closed-form or independently known
results before being trusted: thin-wing lifting-line CL across aspect
ratios, Oswald efficiency, the analytic source-sheet jump, a
source-panelled sphere's exact surface speed, d'Alembert zero force
and the Munk moment on the fuselage, wing-body lift round-trips, and
hard physical bounds (figure of merit and efficiency below one,
vane torque recovery bounded by the jet's angular-momentum flux).
Known limitations are documented in the README and theory
documentation.

## Research impact statement

<!-- TODO: concrete evidence required by JOSS — VBAT UAV program use,
internal reports/design decisions, citations, or other adopters.
Aspirational statements are explicitly insufficient for JOSS review. -->

## AI usage disclosure

Aeolion was developed with substantial AI assistance (Anthropic's
Claude, via Claude Code), with the AI contribution recorded per commit
in the repository history. Correctness was established by the
human-directed validation methodology described above: every module is
checked against closed-form results, external reference data, or hard
physical bounds in the automated test suite, and several defects were
caught precisely by those bounds. <!-- TODO: review wording; state
any AI role in this manuscript itself. -->

## Acknowledgements

<!-- TODO: funding, collaborators, program acknowledgements. -->

## References
