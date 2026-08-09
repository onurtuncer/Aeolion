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

**One real defect found and NOT fixed — worth its own branch:**

`SolveResult::CDi` is not trustworthy on a coupled (wing + closed body)
configuration. Forces are integrated in the near field (Kutta-Joukowski
at each bound-vortex midpoint, see the header comment at the top of
`Solver.h`), and on a coupled solve that integration also collects a
contribution from the discretized closed source bodies. A closed body
carries no net force exactly, but a *panelled* one does, and here the
residual acts as a thrust:

- CDi is **negative** at zero lift (−0.006 at α = −4° where CL ≈ −0.013);
- fitting CDi = CDi0 + k·CL² over the α sweep gives CDi0 = −0.0037 and
  k = 0.0347, i.e. an apparent Oswald e of **1.53**, which is impossible
  (e ≤ 1 for any planar wing).

CL and the moments are unaffected — the lift slope (4.86) and roll damping
(−0.458) both check out against theory. The fix is a **Trefftz-plane
integration** over the wake, which measures the wake's kinetic energy and
is blind to the body's near-field residual. `Solver.h` explicitly notes it
does induced drag "without a separate Trefftz-plane integration"; that is
fine for a wing alone and not fine once a closed body is in the system.
Documented as a stated limitation in `papers/journal-of-aircraft/paper.tex`
Sec. VI.E rather than papered over.

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
