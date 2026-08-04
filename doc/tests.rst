Test Cases and Validation
===========================

Every module is checked against a closed-form or independently-known
result before being trusted, not just eyeballed for plausibility. Each
top-level component owns its own test folder rather than sharing one
central directory -- ``solver/tests/``, ``panelbuilder/tests/``, and
``tests/`` for the core header-only library
(Math/Lattice/Geometry/DragEstimate) -- registered through the shared
``aeolion_add_test()`` in ``cmake/AeolionTest.cmake``. The whole suite is
still one ``ctest`` invocation regardless of which folder a test lives in;
run it with ``ctest --output-on-failure`` from the build directory. Bugs
found and fixed along the way are documented in the relevant code comments
rather than scrubbed from history -- several were caught only by checking
a *hard* physical bound rather than a "looks about right" check. (The
propeller BEMT solver and its physical-bound tests moved to the external
`onurtuncer/BEMT <https://github.com/onurtuncer/BEMT>`_ repository, which
this project does not depend on.)

TestSolverCore
--------------

Validates ``Aeolion::Solver``'s core solver against closed-form thin-wing
theory for a rectangular, unswept, untwisted wing:

- ``CL`` should track the lifting-line-corrected formula
  :math:`C_L = 2\pi\alpha / (1 + 2/AR)` within a few percent across aspect
  ratios 6, 8, 12, 20 (VLM is more accurate than the 2-term formula, so
  convergence in the right direction as AR increases is the real check,
  not an exact match).
- Oswald efficiency for that same wing should sit close to 1.0.
- ``CL``/``CDi`` should converge monotonically as panel count increases.

TestDenseSolve
---------------

Sanity check for the LAPACK/OpenBLAS dense linear-algebra backend: run a
representative solve plus a full stability-derivative sweep and confirm
every result is finite and physically sane.

TestSourcePanel
-----------------

Validates the constant-strength source panel kernel
(``Solver::SourceInfluence``) against results known in closed form, not
against itself:

- **Far field**: at range, a panel of area :math:`A` and unit strength
  must look like a point source of strength :math:`A`, i.e. radial flow of
  magnitude :math:`A/(4\pi r^2)` -- pins both the overall scale factor and
  the sign.
- **On the sheet**: the normal velocity just off a panel's own face must
  be :math:`\pm 1/2`, the analytic source-sheet jump, approached from
  both sides.
- **Sphere**: a source-panelled sphere in uniform flow must reproduce
  potential flow's exact answer, surface speed
  :math:`= \tfrac{3}{2} U \sin\theta`, equivalently
  :math:`C_p = 1 - \tfrac{9}{4}\sin^2\theta`. This exercises the kernel,
  the influence matrix, and the solve together, and is sensitive to any
  error in the solid-angle sign convention.

TestBodyPanels
----------------

Validates the fuselage panelling against physics that holds regardless of
the body's shape, so the checks do not encode the builder's own
arithmetic:

- **Closure**: a capped body of revolution must enclose its stated volume
  and every panel normal must point outward -- a single inverted winding
  would flip that panel's pressure contribution and is otherwise very
  hard to see.
- **d'Alembert**: a closed body in steady potential flow feels no net
  force, at any incidence -- exact, independent of the body's shape, and
  sensitive to the frame conversion, the winding, the kernel sign, and the
  solve all at once.
- **Munk moment**: the same body does feel a pitching moment at
  incidence, and it is destabilizing (nose-up for nose-up incidence) --
  the dominant thing a fuselage does to airframe stability, and the
  reason for modelling it at all.
- **Frame**: the contract states x forward; the solver is x aft. The nose
  must end up ahead of the tail in solver axes.

TestDuctPanels
----------------

The same physics-based checks as TestBodyPanels, applied to the duct
(``Geometry::DuctGeometry``, disjoint from the fuselage -- see
``BuildDuct()``): closure and outward-normal winding across all four
faces (outer wall, inner bore wall, leading/trailing caps), the enclosed
volume matching the annular ring's exact :math:`\pi (r_{outer}^2 -
r_{inner}^2) \times \text{chord}`, d'Alembert (zero net force on the
closed shell), and the frame conversion.

TestPanelBuilder
------------------

Validates that the lattice built from a geometry handoff is a genuine
cambered, chordwise-discretized surface:

- **Reduction**: with one chordwise row and a symmetric (zero-camber)
  section, the builder must reproduce the trusted
  ``ToWingParams()``/``BuildWing()`` lattice essentially to floating
  point, in variants that isolate twist and dihedral separately (the two
  builders' section frames agree exactly when only one of the pair is
  nonzero).
- **Camber physics**: a positively-cambered wing must lift at zero
  geometric incidence, and a symmetric one must not -- the whole point of
  consuming CST sections, and the check a flat lattice structurally cannot
  pass.
- **Chordwise discretization**: requesting M rows must produce M panels
  per spanwise strip, still report one station per strip, conserve total
  planform area, and converge in CL as M increases.
- The on-disk JSON fixture (``tests/Data/AeolionGeometryHandoff-*.json``)
  must build and solve.

TestHandoffContract
----------------------

Parsing, contract invariants, surface binding, and trapezoid reduction for
``Geometry::HandoffContract``: schema version gating, unit checks, the
eta-spanning-sequence rule, CST coefficient-order bounds, control-surface
binding resolution (explicit ``surface`` field vs. the legacy implicit
name convention, which the parser still recovers from on a document that
omits ``surface`` even though no fixture exercises that shape anymore),
deflection-limit consistency, and the body length/station consistency
check -- each one a rejection path a malformed handoff must hit with a
located error message rather than silently produce a wrong lattice.

TestAirframe
-------------

The whole chain on the real airframe: parse the 1.5.0 handoff, build the
wing lattice and the fuselage panels, and solve them coupled. This is the
first place every piece meets, so it is also where integration problems
that no single-component test can see are expected to show up.

Also covers the pattern the viewer's ``Application`` runs every time a
control-surface slider moves: one ``PanelBuilder::LatticeBuilder`` kept
alive across repeated ``ClearDeflections()``/``Deflect()``/``Build()``
calls, re-solved coupled with the SAME cached fuselage + duct panels each
time (the 1.8.0 handoff, which carries a duct). Asserts an antisymmetric
aileron command rolls the coupled airframe without changing its lift, and
that clearing the deflection reproduces the original solve exactly.

TestSectionBoundaryLayer
------------------------

The 2-D transpiration-coupled boundary-layer section solver
(``Solver::BoundaryLayerSectionModel``), checked against external
references rather than itself: a flat plate at zero incidence must not
lift and must land its drag at the Blasius/flat-plate friction level; a
lifting flat plate must sit below but near the :math:`2\pi` thin-airfoil
slope (a boundary layer can only DEcamber); positive parabolic camber
must lift at zero incidence below its inviscid value; drag must fall as
Reynolds number rises; and the converged coupled lift must sit below the
same solver's zero-transpiration (inviscid) pass. See theory.rst,
"Level 3", for the method these checks guard.

TestPropellerLattice
--------------------

The rotating-frame propeller VLM
(``PanelBuilder::BuildPropellerLattice`` + ``Solver::Solve`` with
``FreestreamConditions::p = Omega``), checked against physics that must
hold regardless of the blade: a positively-twisted schedule spun at
:math:`+\Omega` must thrust upstream (:math:`-x`), absorb torque against
its own rotation (:math:`M_x < 0`), cancel its in-plane forces exactly by
blade symmetry, lose thrust as axial inflow rises at fixed rpm, and --
the blade analogue of the wing's camber-lifts-at-zero-incidence check --
thrust at zero twist when its CST sections carry positive camber, which
pins the camber sign convention. Also pins the lattice structure and the
wake direction: every trailing leg leaves along the local relative wind
plus the prescribed momentum-theory hover inflow (a purely axial wake
destroys the lattice geometry at hover; a purely tangential one never
leaves the disk plane).

TestPropellerDuct
-----------------

The ducted-fan interaction (``PanelBuilder::BuildPropellerDuct`` sharing
the viscous-coupled solve as its source body): the shroud mesh must
enclose exactly its annular volume with outward normals
(divergence-theorem closure, the TestDuctPanels style); the coupled solve
must converge with the shroud present; axisymmetry must cancel the
ducted system's in-plane forces; the duct must carry a nonzero share of
the axial force through its pressure field; and the blade circulation
must shift measurably when the shroud is present -- the two-way
interaction being the point of panelling the duct at all. See
theory.rst, "The ducted propeller".

TestPropellerVanes
------------------

The downstream duct-jet vanes (``PanelBuilder::BuildDuctVanes`` +
``Solver::SlipstreamField`` + the static-frame viscous-coupled solve):
mesh/strip alignment; exact cruciform symmetry undeflected (side forces
and control moments at machine zero); undeflected counter-torque
opposing the rotor's :math:`M_x` (the vanes recover swirl); opposite
commands pulling the side force and control moment opposite ways about
the swirl-biased neutral, with a differential that is a real fraction
of the thrust -- the hover control authority thrust-vectoring vanes
exist for; and drag cost for commanding against the swirl. See
theory.rst, "Downstream control vanes and the propulsive wrench".

The mechanics diagnostics run the analytic polar (under the Level-B
wake the hover jet parks the vane tip in reversed jet-edge flow and a
deflected vane's span deep post-stall -- outside any section model's
matching contract; the coupling's reversed-flow residual exclusion and
the deflected cases' coarse plateau tolerance are both documented in
the test). The PROPER boundary-layer vane treatment is pinned by the
mild-jet case: the transpiration-coupled section model on the full
cruciform in a brisk, mildly swirling jet (chord Reynolds number near
3e5, every strip attached), converging tight with exact four-fold
symmetry and counter-torque of the correct sign. Below roughly
Re :math:`2 \times 10^5` the flat plate's laminar-separation polar is
non-monotonic and the coupled cruciform genuinely bifurcates into
asymmetric equilibria -- real physics, noted in the test, but not a
unique fixed point a symmetry assertion can pin.

TestRotorVaneCoupling
---------------------

The two-way rotor-vane coupling (``Solver::SolveRotorVaneCoupled``):
outer convergence of the block Gauss-Seidel alternation; the momentum
budget binding on the high-solidity cruciform (recovered vane torque at
but not beyond the jet's angular-momentum flux -- the one-way scheme
was measured at three times it); the rotor's loading shifting when the
vanes are present (the two-way proof) while remaining a perturbation;
cruciform symmetry surviving; and the control response surviving the
budget, now nearly antisymmetric about neutral since the budgeted swirl
bias is small. See theory.rst, "Two-way rotor-vane coupling".

TestSelfConsistentWake
----------------------

The Level-A wake (theory.rst, the self-consistent trailing-leg pitch):
the outer iteration converges; the converged pitch is a FIXED POINT of
the loading (rebuilding the lattice from the final solved
:math:`v_i(r)` reproduces the final mesh's trailing legs); the pitch
varies radially with the loading rather than following one global
constant; and a heavier disk loading steepens the wake helix -- the load
dependence a hardcoded :math:`\lambda` could not express.

The Level-B wake (theory.rst, "The curved near wake") is pinned here
too: a polyline leg whose vertices lie exactly on the straight
direction reproduces the straight-leg kernel's induction to machine
precision (the kernel-consistency check), and the helical contracted
cored wake must pull the hover figure of merit DOWN versus straight
legs while the rotor still thrusts -- the induced power a straight-leg
wake cannot carry, and the reason the straight-leg hover solve read
optimistic.

TestVaneCascade
---------------

The cascade momentum closure for vane loads (theory.rst, "The cascade
momentum closure"), the coupled driver's default: recovered torque
positive but never beyond the jet's angular-momentum flux -- the bound
that exists BY CONSTRUCTION, where the lattice's capped strips could
still mutually inflate a runaway force; a contracted jet's dead-air
strips carrying exactly nothing (force follows mass flow, closing the
jet-edge pathology); bounded, sign-correct control response in the
swirl-rich jet that bifurcated the lattice closure; and the two-way
coupled solve converging under the closure with no swirl-budget
iteration at all, a single 8-degree vane never out-forcing the
propulsor.

Fixture data
------------

``tests/Data/AeolionGeometryHandoff-{1.5.0,1.8.0}.json`` carries one real
geometry handoff per supported schema revision, so contract-parsing tests
exercise the actual wire format rather than a hand-rolled stand-in. 1.5.0
is the baseline every general invariant test parses; 1.8.0 additionally
covers the duct, moment reference point, blade airfoil sections, and
propeller rotation axis introduced since (see TestSchema180Blocks).
Earlier schema revisions (1.0.0, 1.1.0, 1.4.0) are no longer shipped or
tested against -- this project does not carry a back-compat burden for
retired schema versions, only for the currently-supported ones.
