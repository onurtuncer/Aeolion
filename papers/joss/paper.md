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
  - name: Caglar Ucler # TODO: confirm diacritics (Çağlar Üçler?) and ORCID
    orcid: ''
    affiliation: 2
  - name: Ahmet Gunes # TODO: confirm diacritics (Güneş?) and ORCID
    orcid: ''
    affiliation: 3
affiliations:
  - name: 'Department of Aeronautical Engineering, Istanbul Technical University, Türkiye'
    index: 1
  - name: 'Özyeğin University, Istanbul, Türkiye' # TODO: department
    index: 2
  - name: 'Faculty of TODO, Istanbul Technical University, Türkiye' # TODO: which faculty
    index: 3
date: '' # %e %B %Y
bibliography: paper.bib
---

## Summary

Aeolion is a C++23 aerodynamics toolkit for the conceptual and
preliminary design of propeller-driven aircraft. Lifting surfaces are
modeled with a 3D vortex lattice method (VLM) and bodies of revolution
with constant-strength source panels [@katzPlotkin2001], assembled
into a single blocked linear system solved with LAPACK [@lapack1999],
so that wing, fuselage, and control-surface effects are analyzed
together rather than in isolation: the wing's bound circulation is
carried through the fuselage, and permeable base panels let a
propeller's efflux cross the body without imposing spurious pressure.

Propellers are solved with the same vortex-lattice machinery posed in
the blade-fixed rotating frame: blades are meshed on their CST camber
surfaces [@kulfan2008] wrapped on radius cylinders, and thrust and
shaft torque are read off as the induced drag and rolling moment of
the rotating solve. The wake model is leveled. The trailing legs
depart along the local relative wind with a radially-varying pitch
iterated to self-consistency with the solved loading, and each leg's
first revolutions follow the contracted, accelerating helix a shed
vortex particle traces — prescribed from annular momentum kinematics
in the lineage of the classical prescribed-wake hover analyses
[@landgrebe1972; @kocurekTangler1977] and discretized as polylines of
finite-core (Vatistas
[@vatistas1991]) Biot-Savart segments. This curved near wake carries
the induced power a straight-leg wake cannot: it brings the hover
figure of merit from above the momentum-theory ideal [@leishman2006]
down into the
physically admissible range. On top of the inviscid lattice, a
leveled viscous-inviscid coupling adds sectional lift and profile-drag
feedback per radial strip, up to a transpiration-coupled integral
boundary-layer section solver (Thwaites [@thwaites1949], Head
[@head1958], Squire-Young [@squireYoung1938]) in the
viscous-inviscid interaction tradition of XFOIL [@drela1989], with
the outer fixed
point stabilized by Anderson acceleration [@anderson1965;
@walkerNi2011].

The same solve extends to ducted propulsors: a shroud is paneled into
the rotor solution as a two-way interaction (the duct's boundary
condition sees the propwash, the blades see the bore's constriction),
and all-moving vanes at the duct exit are meshed as lifting surfaces
in the inertial frame, reading the rotor through a slipstream
reconstructed from the converged radial loading by annular momentum
theory [@glauert1935]. Vane loads are closed, by default, with the
cascade momentum
closure a stator row obeys — annulus-by-annulus angular-momentum
bookkeeping, the Euler turbomachine relation posed sector by sector
[@dixon2014], each strip's tangential force exactly the angular
momentum it removes from its sector's mass flow, so every force is
bounded by construction — with a lifting-surface vane closure and its
swirl-budget root solve retained for attached-incidence jets. A block
Gauss-Seidel alternation (\autoref{fig:flowchart}) couples rotor and
vanes two-way — remeshing blades and vanes each pass as the wake
pitch and local flow evolve — under the jet's angular-momentum
budget: the vanes cannot extract more torque than the shaft delivers.
The result is the complete propulsive wrench of a ducted fan with jet
vanes, including control moments at zero airspeed, where the propwash
is the vanes' only dynamic pressure.

![The coupled solution procedure: an outer block Gauss-Seidel
alternation between the rotating and inertial frames. Under the
default cascade closure the swirl budget is satisfied by construction
and the root-solve loop is skipped; under the lattice vane closure the
budget is solved each pass as a bracketed root problem on the swirl
factor $s$ with the rotor state
frozen.\label{fig:flowchart}](figures/solution-flowchart.png)

Geometry is loaded from a versioned JSON handoff contract rather than
an STL/mesh file, keeping the description parametric (planform
stations, CST camber, control-surface hinges, blade sections) and
traceable to the CAD source that generated it. An OpenGL viewer
renders the solved lattices with interactive control-surface and vane
deflection, live propulsive-wrench readouts, and scripted multi-angle
screenshots for documentation.

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
rotor but carry no airframe. The configuration ingredients themselves
are classical — shrouded propulsors and their static thrust
augmentation [@kuchemannWeber1953; @pereira2008], jet-vane thrust
vectoring [@sutton2016], swirl-recovery vanes [@gazzaniga1992] — but
their coupled analysis at hover is not something the established
tools provide. Aeolion's contribution is not a new
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

Three design decisions carry most of the physics architecture. First,
the frame split: blades (and the duct, a body of revolution about the
rotation axis) solve in the rotating frame, exit vanes solve static,
and the two communicate through a reconstructed time-mean slipstream
rather than instantaneous induction — blade-passing effects are
deliberately outside the model. Second, leveled fidelity on explicit
seams: the wake ladder (straight legs, self-consistent pitch,
prescribed curved near wake — with force-free relaxation and roll-up
as documented next rungs) and the 2-D section model (analytic polar
or integral boundary-layer solver behind one interface, each assigned
where its physics is honest) can be climbed independently without
touching the coupling loops. Third, degradation is explicit rather
than silent: deep stall hands off to a saturated polar, reversed-flow
strips are excluded from convergence criteria, and every such
boundary is documented and test-pinned.

Every module was validated against closed-form or independently known
results before being trusted: thin-wing lifting-line CL across aspect
ratios, Oswald efficiency, the analytic source-sheet jump, a
source-panelled sphere's exact surface speed, d'Alembert zero force
and the Munk moment on the fuselage, wing-body lift round-trips, and
hard physical bounds — figure of merit and efficiency below one, vane
torque recovery bounded by the jet's angular-momentum flux. Several
defects were caught precisely by such bounds. Known limitations are
documented in the README and theory documentation.

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

This work was supported by the Istanbul Technical University
Scientific Research Projects Coordination Unit (İTÜ BAP) under
project MGA-2026-48282, "JETTAIL: Design, Control, and Autonomous
Mission Capability Development of a Jet-Flow-Vectored EDF-Based
Tail-Sitter VTOL UAV."

<!-- TODO: any further funding, collaborators, program
acknowledgements. Verify the English rendering of the project title
against the official translation, if one exists. -->

## References
