# TODO — stagnation points and attachment lines

Handoff notes for continuing on another machine. Branch:
`feature/stagnation-attachment-lines`.

Goal of this work: identify stagnation streamlines on the wing and on the
body so a boundary-layer method can be coupled to them later, across a
sweep of angle of attack **and sideslip**.

---

## What is done and verified

All seven solver tests pass locally (see "How this was verified" below for
the caveat about the local toolchain).

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

### 1. Verify the build under the real toolchain — do this first

The local machine had **MSVC but no vcpkg, LAPACK, or nlohmann/json**, so:

- Solver tests were built and run against a **hand-written reference
  `dgetrf_`/`dgetrs_` shim**, not real LAPACK. The shim is a
  straightforward column-major partial-pivoting LU and all seven tests
  pass against it, but **re-run under OpenBLAS** to confirm.
- **`panelbuilder/src/PanelBuilder.cpp` was never compiled** — it needs
  `nlohmann/json.hpp`. The change there is three added lines, all of the
  form `panel.SectorIndex = j;`. Low risk, but unverified.
- The panelbuilder test suite (`TestBodyPanels`, `TestAirframe`, etc.)
  was likewise not run. `SourcePanel` gained a defaulted field, which
  should not disturb them.

```
cmake --preset windows && cmake --build build/windows && ctest --test-dir build/windows
```

### 2. Sphinx build

`doc/theory.rst`, `tests.rst`, `api.rst` and `references.bib` were edited
but **Sphinx was never run** (no Python doc toolchain locally). Check the
new math renders and the new `:cite:` keys resolve.

### 3. Finish the paper — Section VI "Application"

The only substantive gap. It needs an α/β sweep on the real geometry
handoff (`tests/Data/AeolionGeometryHandoff-1.8.0.json`), which needs
`panelbuilder` and therefore vcpkg. Planned content is listed as a
comment block in `paper.tex`:

- attachment line vs. span at α = 0, 4, 8°, β = 0
- the same at β = ±10°, showing the left/right asymmetry
- `Rbar` vs. span against the 245 / 583 thresholds
- fuselage skin-flow topology: nose attachment node migrating around the
  windward meridian with combined α and β, plus the leeward separation
  line
- figure: traced surface streamlines on the fuselage, coloured by Cp

Also outstanding in the paper: author block (co-authors? the SciTech
draft has three), AIAA member grades, acknowledgments/funding, and
regenerating `style/aiaa-tc.cls` + `aiaa.bst` via `latex aiaa.ins` before
any real submission (see `papers/journal-of-aircraft/README.md`).

`papers/journal-of-aircraft/README.md` still says **Status: planned** and
describes a different scope (wing–body–propeller application study). It
should be updated to match what `paper.tex` now is.

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
  Nose-clustered station spacing in `BuildBody` would resolve it; noted
  in the test as the reason for its body/incidence choices.
- `AttachmentLine` takes one Weissinger row per strip (same contract as
  `SolveViscousCoupled`). Multi-row chordwise lattices need the caller to
  pass the leading-edge row.

---

## How this was verified

```
# from the repo root, with MSVC on PATH
cl /nologo /EHsc /O2 /std:c++latest /I include /I solver/include \
   solver/tests/TestSurfaceFlow.cpp <lapack-shim>.cpp /Fe:TestSurfaceFlow.exe
```

Tests passing locally: `TestSolverCore`, `TestDenseSolve`,
`TestSourcePanel`, `TestSectionBoundaryLayer`, `TestSurfaceFlow`,
`TestSectionPanelMethod`, `TestAttachmentLine`.

The three new tests check against closed-form answers, not against
themselves — sphere stagnation points exact at any α/β, cylinder surface
velocity exact, the `√r_LE` collapse, and the swept-vs-unswept sideslip
asymmetry with its control case.
