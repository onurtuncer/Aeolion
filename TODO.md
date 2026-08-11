# TODO — stagnation points and attachment lines

Branch: `feature/stagnation-attachment-lines`.

Goal of this work: identify stagnation streamlines on the wing and on the
body so a boundary-layer method can be coupled to them later, across a
sweep of angle of attack **and sideslip**.

What remains is section 4 (the boundary-layer coupling itself, which
deliberately wants its own branch) and section 5's smaller items; the
build, the docs and the paper are all current — see "How this was
verified" at the bottom.

---

## What is done and verified

All 21 test suites pass under the real toolchain.

### Plumbing

- **`Lattice::SourcePanel::SectorIndex`** — new field. Together with the
  existing `StationIndex` it makes a body-of-revolution's panelling a
  structured *(station, sector)* grid. Set by `PanelBuilder::BuildBody`
  and by the duct's `AnnularRingPanels`; deliberately left unset on the
  fuselage base cap, which is a disc stacked in rings at one station and
  so has no unique key (documented in place).
- **`SolveResult::sigma`** — converged source strengths are now kept, not
  discarded. `gamma` + `sigma` is the full state of the solved field.
- **`Solver::FlowField`** (in `Solver.h`) — the velocity field as a
  first-class object, evaluable anywhere after the solve. Three
  evaluations: `Velocity(P)`, `BoundMidpointVelocity(i)` (own bound
  segment excluded), `SourceSurfaceVelocity(k)` (analytic ½ sheet jump
  substituted). **`SolveWithSystem` was refactored to go through it**, so
  there is no second definition left to drift — this is the change most
  worth re-reviewing, and `TestSolverCore` / `TestSourcePanel` /
  `TestDenseSolve` all still pass against it.
- **`FreestreamVelocity(fc)`** — moved the duplicated freestream
  expression into `FreestreamConditions.h`.

### Body — `solver/include/Aeolion/Solver/SurfaceFlow.h`

Skin-flow topology of a source-panelled surface, in the panelling's own
index space (never `atan2(z,y)`, which is silently wrong for an offset
duct).

- `BuildSurfaceGrid` — structured patch + metric + contravariant velocity.
  Declines an incompletely indexed surface rather than guessing.
- `FindCriticalPoints` / `AnalyzeSurfaceFlow` — locates and classifies
  attachment/separation nodes, saddles, foci.
- `TraceSurfaceStreamline` / `TraceFromCriticalPoint` — RK4 in index
  space, carrying `U_e(s)` and the spreading metric `h(s)`.
- `PrimaryAttachment` / `PrimarySeparation`, `AttachesUpstream` for the
  zero-incidence apex case.

### Wing

- **`include/Aeolion/Geometry/SectionContour.h`** — the thick closed
  contour from the same CST coefficients the camber line uses, plus
  `r_LE/c = A0²/2` and `ToLeadingEdgeNormal`.
- **`solver/include/Aeolion/Solver/SectionPanelMethod.h`** — Hess–Smith
  source+vortex solve on that contour. Returns `U_e(s)`, the stagnation
  point, its strain rate, and the two `SurfaceRun`s a BL march consumes.
- **`solver/include/Aeolion/Solver/AttachmentLine.h`** — per-strip
  attachment line in **leading-edge-normal** coordinates, with effective
  sweep, Poll's `Rbar`, and root-kink detection.

### Docs

- `doc/theory.rst` — new section "Stagnation points and attachment lines"
  (inserted before "Viscous drag buildup", ~line 1184).
- `doc/tests.rst` — entries for the three new tests.
- `doc/api.rst` — new types listed under Solver and Geometry.
- `doc/references.bib` — Poll, Hess & Smith, Moran, Lighthill, Tobak &
  Peake, Cebeci & Cousteix.

### Paper

- `papers/journal-of-aircraft/paper.tex` — full draft, AIAA JoA format
  (`\documentclass[submit]{aiaa-tc}`), plus its own `references.bib`.
  Sections I–V are written against implemented, tested code.

---

## Two physics results worth not re-deriving

1. **Stagnation-point offset scales with `√r_LE`, not `r_LE`.**
   `s_stag/c ~ √(2 r_LE/c) · α_e`. The intuitive "the nose is a cylinder
   in a stream at angle α" reading gives `r_LE·α` and is wrong by an
   order of magnitude (0.001c vs a true 0.013c for a NACA 0012 at 4°).
   The nose sits inside the *outer* thin-airfoil leading-edge
   singularity, `u ~ α√(c/s)`; matching against the parabolic nose's own
   `√(s/r_LE)` gives the square root. Measured collapse constant ≈ √2,
   verified across a 12× range of nose radius.

2. **Surface strain rates must be differenced in a local orthonormal
   frame**, not from the contravariant components `(u,v)`. Those carry a
   `1/r` chart factor near a nose that varies tens of percent per cell;
   differencing them gives a sphere's *isotropic* stagnation node a 2:1
   eigenvalue split that does **not** converge. With the orthonormal
   frame the ratio converges 0.28 → 0.67 → 0.85 → 0.94 for N = 24, 32,
   48, 64.

---

## Remaining work

### 1. ~~Verify the build under the real toolchain~~ — DONE (2026-08-08)

Full `cmake --preset windows` build and all **21** ctest suites pass on
the vcpkg machine against real OpenBLAS/LAPACK; `PanelBuilder.cpp`
compiles with real nlohmann/json and the whole panelbuilder suite is
green. One warning fixed along the way: dead `cosAlpha`/`sinAlpha`
locals in `SectionPanelMethod.h`.

### 2. ~~Sphinx build~~ — DONE (2026-08-08)

Doxygen XML + `sphinx-build -b html` succeed; all 12 `references.bib`
entries parse, every new `:cite:` key resolves (Poll renders as a linked
`[Pol79]`), and the new theory section is present in the output. The
only warnings are graphviz-missing (local machine; CI installs it).

### 3. ~~Section VI "Application"~~ — DONE (2026-08-08)

Written against a real α/β sweep: `aeolion_attachment_sweep`
(app/AttachmentSweepExport.cpp) runs the coupled solve on the 1.8.0
handoff across α ∈ {0,4,8}° × β ∈ {−10,0,+10}° and exports
`papers/journal-of-aircraft/figures/attachment-sweep.json`;
`render-attachment-figures.py` draws the four figures; the paper
compiles to 18 pages with everything resolved. Findings worth knowing
when rereading:

- The fixture wing is exactly rectangular/unswept, so all effective
  sweep is flow-induced: ~30° at the root junction at β = 0 (the body's
  crossflow deflected around the wing root), decaying to <1° outboard;
  `Rbar` peaks at the root (86–104) and stays a factor 2.4 below the
  contamination threshold everywhere.
- Sideslip splits the roots (windward keeps the spike, leeward passes
  through a spanwise-flow null one station out) and leaves a ~2°
  outboard asymmetry that grows with CL — induced, not geometric.
- The body's station list is nose-refined in the driver (consumer-side
  resampling of the contract's own radius law); with it the attachment
  node resolves at 8 of 9 conditions (α = β = 0 honestly
  AttachesUpstream at the apex).
- MEASURED: raising BodyCircumferentialPanels to 24 walks flank control
  points onto the no-carry-through wing's root trailing-leg line —
  pivot ratio ×5 worse, CL inflated 80%. Keep 16, or keep control
  points off that line (comment in AttachmentSweepExport.cpp).

Still outstanding in the paper: author block (co-authors? the SciTech
draft has three), AIAA member grades, acknowledgments/funding,
external validation anchors, and regenerating `style/aiaa-tc.cls` +
`aiaa.bst` via `latex aiaa.ins` before any real submission (see
`papers/journal-of-aircraft/README.md`, whose status/open items are now
current).

### 3a. Separation and the coefficient/derivative matrix — DONE (2026-08-09)

Section VI now reports what the whole method was for: a matrix of
coefficients and stability derivatives over α ∈ [−4, 16]° × β ∈ [−10,
10]°, restricted to attitudes where the flow is still attached, plus the
separation boundary itself.

- **New module** `solver/include/Aeolion/Solver/AttachmentBoundaryLayer.h`
  marches each strip from its *real* attachment point (Thwaites → Michel
  → Head/Ludwieg-Tillmann) and reports separation. Additive: it does not
  touch `SectionBoundaryLayer.h`, so the tested fixed point is unchanged.
  It may latch the first crossing precisely because nothing iterates on
  it. `TestAttachmentBoundaryLayer` (new, 22nd suite) pins it on Blasius,
  Howarth (s/L = 0.123, mesh-converged), and an unseeded Hiemenz θ₀.
- **Separation criterion.** At Re_n ≈ 3e5 the laminar layer reaches
  λ = −0.09 at *every* attitude including negative α, so "laminar
  separation" draws no boundary. The bubble is treated as the transition
  trigger and the verdict is turbulent separation (H ≥ 2.4). Onset is
  α = 6° at every β tested (34% of span, x/c ≈ 0.88), reaching the whole
  span by α = 8°. Bubble *bursting* is not modelled — stated in the paper
  as making this an upper bound on usable incidence, not a stall
  prediction.

**Two bugs found and fixed on the way** (both pre-existing, neither
pinned by any test):

1. `StabilityDerivatives`' reduced-rate fields multiplied by
   `length/(2V)` where the convention Cl_p = ∂Cl/∂(pb/2V) requires
   `2V/length` — the reciprocal, wrong by (2V/b)² ≈ 2000 at these
   numbers. It reported Cl_p = −2e−4 for an AR=6 wing whose textbook
   value is −0.45; with the factor corrected it reads −0.452.
   `TestSolverCore` now pins both the value and the identity, and also
   pins that Cm_q is *exactly* zero for a single-row wing about its own
   quarter chord (correct, not a missing term) while being negative about
   a point two chords aft.
2. `AttachmentSweepExport.cpp` passed the contract's moment reference
   point into `FreestreamConditions::RefPoint` without the contract →
   solver frame flip (x_solver = −x_frd). That put the reference point an
   equal distance the wrong side of the origin — a spurious moment arm of
   ~0.48 m, about one body length, inflating Cm_α by an order of
   magnitude. Mine, introduced in this branch.

**One real defect found — NOW FIXED on `fix/trefftz-induced-drag`:**

`SolveResult::CDi` is not trustworthy on a coupled (wing + closed body)
configuration. Forces are integrated in the near field (Kutta-Joukowski
at each bound-vortex midpoint, see the header comment at the top of
`Solver.h`), and near-field induced drag is a small difference of much
larger quantities — the streamwise component of forces dominated by lift,
so its relative error scales with L/D. Adding a body changes the induced
velocity at the wing's bound vortices and adds a pressure integration over
the body in a non-uniform field:

- CDi is **negative** at zero lift (−0.006 at α = −4° where CL ≈ −0.013);
- fitting CDi = CDi0 + k·CL² over the α sweep gives CDi0 = −0.0037 and
  k = 0.0347, i.e. an apparent Oswald e of **1.53**, which is impossible
  (e ≤ 1 for any planar wing).

CL and the moments are unaffected — the lift slope (4.86) and roll damping
(−0.458) both check out against theory.

NOT a d'Alembert violation — that is the tempting explanation and it is
measurably false. Closed bodies here carry zero net force in uniform flow
to **machine precision** (|F|/qA ~ 1e-16, TestBodyPanels/TestDuctPanels),
and the body's force in a coupled solve is physical: it sits in the wing's
upwash and carries ~8% of the lift.

**Fixed** by `solver/include/Aeolion/Solver/TrefftzPlane.h`: a far-field
integration over the wake trace, which never evaluates the body at all
because a closed body sheds no wake. On the coupled airframe the
zero-lift intercept falls −0.0037 → −0.0003 and CDi is positive at every
attitude. `TestTrefftzPlane` pins it on elliptic loading (e = 1 and
CDi = CL²/(πAR)), the e ≤ 1 bound across five distributions, exact
cancellation of chordwise stacks, and agreement with the near-field
method to 0.4% on a wing alone.

**The gotcha that cost an hour, so it is written down:** span efficiency
must be formed with the LIFTING SYSTEM's lift, not the configuration's.
Only the lifting system sheds a wake, so the far-field drag is its alone.
The fuselage here carries ~9% of the lift, and comparing the whole CL
against a wing-only CDi gives e = 1.17 — still above the bound, and
looking exactly like a defect in a perfectly sound integral. Against the
wing's own lift (`SolveResult::LiftBySurface["wing"]`) it is 0.965, with
per-condition values 0.96–0.99.

Still open: `SolveResult::CDi` remains the near-field number, with the
Trefftz result computed alongside rather than replacing it. Making the
far-field value the default is a behaviour change to a tested field and
wants its own decision.

### 3b. Third paper — aft-fan inflow induction (scoped 2026-08-09)

Scoped in `papers/journal-of-aircraft-fan-induction/README.md`; no draft
and no code yet. Target: a second Journal of Aircraft article asking **how
far the aft ducted fan delays wing separation during tail-sitter
transition**.

Two findings from scoping that are worth not rediscovering:

- **The fan is AFT of the wing** on the 1.8.0 handoff — duct LE at solver
  x = 0.422, wing TE at 0.384, a gap of 0.21 chord. So the mechanism is
  the fan's *upstream induction* (a favourable gradient over the wing),
  **not** propwash blowing. Anyone scoping this as a blown-wing study is
  describing a different aircraft.
- **`Solver::SlipstreamField` cannot be used for it**: it returns zero for
  `point.x < 0` by construction (momentum-theory wake only), so it would
  report exactly zero effect on a wing that sits upstream of the disk —
  an artifact, not a result. An upstream induction model has to be built;
  the semi-infinite vortex cylinder has a closed-form upstream field and
  is what would make the paper verifiable the way the first two are.

The alignment that makes it a paper: duct outer radius / semi-span = 0.21,
and the second paper's worst separation station is |2y/b| ≈ 0.15–0.23. The
fan sits over exactly the span the wing sheds first.

Deliberately NOT written: post-stall coefficients to 90 deg. No
closed-form verification exists past separation, and the answer would be
set by the four hand-tuned constants of `AnalyticSectionModel`'s
deep-stall blend rather than by the method. Recorded here so the decision
is not relitigated.

### 3c. Post-separation study — Phase 0 (2026-08-10)

A separate study from the papers: where does the configuration's behaviour
converge to flat-plate scaling in alpha AND sideslip? Phase 0 built the
scaffolding and the convergence metrics, ran the Level-2 coupled solve far
past anything it had seen (alpha in [-4, 90] x beta in [0, 30]), and fixed
what broke. The deep-stall NUMBERS are still set by AnalyticSectionModel's
hand-tuned constants — the 3b caveat stands; Phase 1 replaces them.

**New driver** `aeolion_poststall_sweep` (app/PostStallSweepExport.cpp):
Level-2 coupled solve on the CLEAN lattice (a carry strip's bound midpoint
is inside the fuselage where the source field is the interior continuation),
warm-start continuation up each alpha column, body+duct sources coupled.
Exports per condition: CN/CC/CL/CD/Cm, total and wing-only force/moment
vectors, force angle off the chord-plane normal, sigma (total inclination),
xcp, per-strip [eta, alpha_eff, cl, cd, residual], coupling diagnostics.
CLI: handoff, out, Vinf, relaxation, andersonDepth, maxIterations, [beta].
Reproduce the map:
`aeolion_poststall_sweep tests/Data/AeolionGeometryHandoff-1.8.0.json out.json 25 0.05 0 1000`

**Three failure modes found, all understood, two fixed:**

1. **Spanwise checkerboard multistability.** The Anderson-accelerated fixed
   point on 44 tightly packed wing strips lands on sawtooth equilibria
   (alpha_eff alternating +-5 deg strip to strip; alpha=4 gave CL 0.35 or
   0.93 depending on start). These are the classic spurious equilibria of
   collocation nonlinear lifting-line. Plain damped iteration (omega=0.05,
   AndersonDepth=0) converges the whole attached range to residual < 1e-4
   in ~110 iterations; the propeller consumers (12 strips, kinematic-
   dominated) never see this. NOT fixed in the solver — driver passes the
   options. A spanwise-smoothed or Newton update is the real cure if the
   wing becomes a first-class Level-2 consumer.
2. **Camber double-count in strip frames** (driver-side, fixed). Building
   ChordDir/LiftDir from the cambered panel geometry absorbs the lattice's
   zero-lift shift, and Alpha0Deg then subtracts camber again: measured
   coupled zero-lift at -7.8 deg = lattice -3.9 + thin-airfoil -4.2. The
   driver now uses the true chord frame; coupled CL matches the inviscid
   lattice to 4 digits at zero lift. The exact trap the ViscousCoupling.h
   header warns about; a future PanelBuilder wing-strip builder must carry
   the section plane for swept/twisted wings.
3. **cl->Gamma inversion degeneracy at |alpha_eff| ~ 90 deg** (solver,
   fixed in ViscousCoupling.h). The circulatory force is perpendicular to
   the lift direction there, k passes through zero, and cl/k slammed the
   target between +-gammaCap with the sign of k's noise (CN ~ 5 at
   alpha=90, artifact), while the |k|~0 branch froze stale continuation
   circulation. Fix: targets ramp to zero beyond the residual's own
   contract edge (ResidualIncidenceLimitDeg, TargetDecayRampDeg=10). A
   Tikhonov-damped inversion was tried first and REJECTED: its 0.25% bias
   floors the residual above tolerance — TestViscousCoupling and
   TestPropellerDuct caught it. With the decay, all suites pass and the
   rotor-vane suites run ~20x faster (reversed vane tips stop chasing the
   degenerate inversion): TestRotorVaneCoupling 137 s -> 7 s,
   TestVaneCascade 191 s -> 17 s.

**Deep-stall limit cycles are inherent, and the cycle means are
reproducible**: independent iteration paths (Anderson vs plain damped)
agree to 4+ digits on the cycle-mean loads at alpha 25..80. The map's
deep-stall values are cycle means, exported with Converged=false and the
residual — by design, not laundering.

**The Phase-0 map** (donated constants and all): cross-beta collapse in
total inclination sigma (sin sigma = sin alpha cos beta) holds to <5%
spread from alpha ~ 6 deg through 75 deg — sideslip up to 30 deg only
rescales the loads through cos beta. CN/sin sigma decays from ~7.5
(attached) to the plate plateau ~2.1-2.3 by alpha ~ 55-65. CLmax = 1.39 at
alpha = 20 (the analytic blend stalls 14 deg later than the computed
separation onset at 6 — the gap Phase 1's Kirchhoff bridge closes).
CN(90) = 2.11 vs Viterna CDmax(AR=6) = 1.22: +73%, the quantified cost of
the AR-blind PlateNormal=1.8. xcp is structurally pinned at ~0.25c:
SectionCoefficients has no cm, so the strip force acts at the quarter
chord and the plate's walk to mid-chord CANNOT be represented — Phase 1
must add cm to the section interface. Above alpha ~ 80 (and beta >= 15)
the strip contract itself dies (most strips beyond the incidence limit);
that corner of the map is scaffolding, not physics.

**Phases agreed** (chat, 2026-08-10): 1 — anchored post-stall section
models (Viterna AR-aware CDmax + Hoerner CN as the deep anchor, Kirchhoff
attenuation driven by AttachmentBoundaryLayer's computed separation point,
section cm, validation against Sheldahl & Klimas Re=3.6e5 / Ostowari-Naik;
plus hysteresis map via up/down continuation). 2 — Maskew-Dvorak double
wake on SectionPanelMethod's Hess-Smith solve. 3 — vortex particles (2D
LESP discrete-vortex sections first, 3D particle wake from the computed
separation line as an unsteady spot-check).

### 3d. Phase 1 + the strip-frame paper correction (2026-08-11)

**Phase 1 landed.** `Solver/PostStallSection.h`: Kirchhoff attenuation
K(f) on the exact thin-airfoil cn/cc pair (both classical limits exact:
f=1 is d'Alembert-clean, f=0 the plate quarter-slope), f interpolated
per strip from AttachmentBoundaryLayer tables built over an inviscid
alpha sweep; Viterna-Corrigan beyond the EMERGENT stall with AR-aware
CdMax (2.01 at the AR=50 edge vs Hoerner's 1.98); Rayleigh's
free-streamline xcp for the section cm. `SectionCoefficients` gained cm
(additive; every existing consumer bit-identical) and the coupling
applies it as a pure quarter-chord couple
(`ViscousCoupledResult::SectionMoment`). `TestPostStallSection` (25th
suite) pins the anchors, the exactness of the attached limit, emergent
stall, junction continuity, and the couple identity. Docs: theory.rst
section, api.rst, tests.rst, viternaCorrigan1982 in references.bib.

Anchored map vs Phase-0 baseline (beta=0): CN(90) 2.11 -> 1.52 against
the Viterna anchor 1.218 (+25% residual: local dynamic pressure at the
inboard strips + cycle-mean circulation); xcp walks 0.25c -> 0.52c and
LANDS ON Rayleigh's mid-chord — a metric the model was never fitted to;
CLmax 1.76 at alpha=18, emergent from the separation tables and an
upper bound (bubble bursting, see 3a); sigma-collapse across beta
within 13% (baseline 6% — the sharper stall cycles harder).

**Strip-frame bug, and the papers.** The Phase-0 camber double-count
(3c, finding 2) also lived in AttachmentSweepExport: ComputeAttachmentLine
measures alpha_n in the strip frame, the driver supplied camber-tilted
panel axes, and the section contour carries the camber again — alpha_n
biased HIGH by ~4.3 deg. Fixed (true chord frame). Consequences,
verified by full-matrix diff: CL/Cm/derivatives BIT-IDENTICAL; the
separation onset moves 6 -> 12 deg (55% span at 12, 82% at 14, ~full at
16 with forwardmost x/c = 0.80); worst Rbar 104 -> 66 within the
attached envelope (contamination margin 2.4x -> 3.7x), 135 over the
whole matrix; stagnation offsets roughly halve at low alpha. SciTech
lattice-solution.json: bit-identical after all solver changes (verified
numerically) — no paper impact. Fan-induction fine sweeps regenerated
with the corrected frames (delay table to be re-rendered; the Delta
alpha result is a same-frame difference and needs re-measuring, not
assuming).

**Papers restructured as Part I / Part II** (user direction,
2026-08-11): Part I = papers/journal-of-aircraft, "...Across the
Separation Boundary — Part I: Attachment Lines, the Separation March,
and the Attached-Flow Envelope" — everything up to the separation
boundary, coefficient tables truncated there, all prose numbers
corrected. Part II = papers/journal-of-aircraft-poststall (new) — the
post-separation study as a paper: anchored section model, the
alpha x beta map to 90/30, sigma-collapse, plate convergence, the two
anchors. Both compile clean (26 pp / draft). Author blocks must stay
synchronized.

### 4. The actual boundary-layer coupling

This work deliberately stopped at the *prerequisite*. Everything a march
needs is now produced but nothing consumes it yet:

- `SectionSolution::UpperRun()` / `LowerRun()` give `U_e(s)` from the
  stagnation point; `StagnationMomentumThickness` gives `θ₀`.
- `SurfaceStreamline` gives `U_e(s)` and `h(s)` on the body.
- `AttachmentStation` gives `Rbar`, so the march knows whether it may
  start laminar at all.

`SectionBoundaryLayer.h` currently starts its march at the camber-line
leading edge with `θ = 0` on both surfaces — i.e. it assumes the
stagnation point is at `x/c = 0`, which is exactly the approximation this
work removes. Wiring it to start at the real attachment point with the
real `θ₀` is the natural next step, and is a behaviour change to an
existing tested module, so it wants its own branch.

### 5. Smaller items

- Body streamline tracing assumes the *lateral* fuselage surface; the
  base cap is (correctly) declined. If the attachment analysis should
  ever cross onto the base, that needs a second patch and a join.
- `MinSurfaceStations`/`MinSurfaceSectors` and the nose-resolution issue:
  on a slender body at low incidence the stagnation point can fall inside
  the first panel ring, where the honest answer is `AttachesUpstream`.
  **Worked around consumer-side** in `AttachmentSweepExport.cpp`
  (`RefineNoseStations` cosine-clusters the contract's station list toward
  the apex, keeping every original breakpoint so the shape is unchanged),
  which resolves the node at 8 of the paper's 9 conditions. Whether
  `BuildBody` should offer this as a `LatticeOptions` knob rather than
  leaving every caller to reimplement it is still open — the argument for
  the knob is that axial spacing is exactly as much a consumer choice as
  `BodyCircumferentialPanels` already is.
- `AttachmentLine` takes one Weissinger row per strip (same contract as
  `SolveViscousCoupled`). Multi-row chordwise lattices need the caller to
  pass the leading-edge row.

---

## How this was verified

Current, on the vcpkg machine (real OpenBLAS/LAPACK, nlohmann/json):

```
cmake --preset windows && cmake --build --preset windows && ctest --test-dir build/windows
```

**21/21 tests pass** — the seven solver suites plus the whole panelbuilder
suite (`TestBodyPanels`, `TestAirframe`, `TestVaneCascade`,
`TestRotorVaneCoupling`, …), the logger, and the contract tests. Note the
build must run inside the MSVC dev environment with `VCPKG_ROOT` re-set
after `vcvars64.bat`, which overrides it to VS's bundled vcpkg.

Docs: Doxygen XML + `sphinx-build -b html` succeed against
`doc/requirements.txt` (graphviz is the only local gap; CI installs it).

Paper: `pdflatex && bibtex && pdflatex && pdflatex` in
`papers/journal-of-aircraft/` with `TEXINPUTS=./style;` — 18 pages, no
undefined references. (`latexmk` does not work on this machine; no perl.)

The three new solver tests check against closed-form answers, not against
themselves — sphere stagnation points exact at any α/β, cylinder surface
velocity exact, the `√r_LE` collapse, and the swept-vs-unswept sideslip
asymmetry with its control case.
