Test Cases and Validation
===========================

Every module is checked against a closed-form or independently-known
result before being trusted, not just eyeballed for plausibility. The
suite is wired into ``ctest`` (``tests/CMakeLists.txt``); run it with
``ctest --output-on-failure`` from the build directory. Bugs found and
fixed along the way are documented in the relevant code comments rather
than scrubbed from history -- several were caught only by checking a
*hard* physical bound rather than a "looks about right" check, which the
BEMT tests below rely on deliberately.

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

TestBEMT
--------

Validates ``Aeolion::BEMT`` against **hard physical bounds**, not
plausibility checks:

- Figure of Merit (ideal hover power / actual power) must be :math:`\le
  1.0` -- a real thermodynamic constraint, and a strong bug detector: this
  caught two real sign/exponent bugs during development.
- Propulsive efficiency (:math:`T V_\infty / P`) must be :math:`\le 1.0`
  for every forward-flight point where the prop is actually absorbing
  power.
- Every blade station must converge.

TestPropVane
------------

Integration test tying ``BEMT::SlipstreamField`` into
``Solver::Solve``'s ``externalField`` hook: a deflected vane sitting
behind a running propeller, at essentially zero airspeed (hover), must
produce real side force / yaw moment purely from the propwash, while the
same vane with no propwash present produces ~nothing. This is the literal
physical mechanism that gives a thrust-vectoring-vane aircraft hover
control authority, so it is kept as a standing regression check.

TestPropContract
-----------------

The handoff-to-propeller bridge (``Geometry::ToPropGeometry``): the
contract states blade radii as fractions of the disk radius and BEMT
wants metres, and unit-conversion code is exactly the kind that looks
obviously right and silently isn't -- a factor-of-R error would shift
every radius and still produce a plausible-looking thrust. Also covers
the hub station, which broke the solve the first time this bridge was
wired up.

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
binding resolution (explicit ``surface`` field vs. the 1.0.0 name
convention), deflection-limit consistency, and the body length/station
consistency check -- each one a rejection path a malformed handoff must
hit with a located error message rather than silently produce a wrong
lattice.

TestAirframe
-------------

The whole chain on the real airframe: parse the 1.4.0 handoff, build the
wing lattice and the fuselage panels, and solve them coupled. This is the
first place every piece meets, so it is also where integration problems
that no single-component test can see are expected to show up.

Fixture data
------------

``tests/Data/AeolionGeometryHandoff-{1.0.0,1.1.0,1.4.0,1.5.0}.json``
carries one real geometry handoff per schema revision, so contract-parsing
tests exercise the actual wire format rather than a hand-rolled stand-in,
and the schema's evolution (body, control-surface binding, wing
placement) is tested against the version it was introduced in.
