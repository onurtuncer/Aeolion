# TODO

Branch: `study/post-separation`.

## How to read this

This file is **forward-looking**: what is still open, roughly in the
order it is worth doing. It replaces the running study log that occupied
it through 2026-08-19. That log is preserved in git —
`git show cfbfc75:TODO.md` — and its durable content has been moved to
where it is actually read:

- **measured results and traps** now live in the module headers that own
  them (`Solver/BodyAxes.h`, `Solver/ViscousCoupling.h`,
  `Solver/ParticleWake.h`, `app/AileronEffectivenessExport.cpp`, …) and
  in the suites that pin them;
- **the flight model's decisions of record** are in `models/README.md`,
  which is normative;
- **the study's narrative** is in the Journal of Aircraft drafts.

A short index of settled findings is at the bottom, so nothing is
re-derived by accident.

---

## A. Promises not kept

Places where a shipped artifact claims more than it delivers. First
because they are cheap, and because the cost of leaving them is a
consumer trusting something that is not there.

- [x] **A1. Uncertainty bounds carry real data** — DONE 2026-08-20
  (`c5f680b`). The solver measures the cycle WIDTH
  (`CycleFluctuation()`), the drivers export it, and the assembler emits
  12 DAVE-ML `uncertainty`/`normalPDF` elements with per-cell additive
  bounds. 65 converged rows correctly get no bounds; 110 cycle-mean rows
  get them, median 1.14e-7 but reaching **0.591 at α = 90, β = 30 against
  CZ = −1.485** — a 40% bound, which is the honest number there.

  The bound is deliberately the cycle width and NOT the residual: the
  residual is a max over strips of a section mismatch and sits at 0.2–0.6
  post-stall, while the load-level fluctuation is ~1e-6, so quoting it
  would overstate uncertainty by five orders of magnitude.

- [x] **A2. DTD check could silently skip in CI** — DONE:
  `libxml2-utils` added to `sanitizers.yml`. *Original:* `verify-daveml.py`
  uses `xmllint` when present and says so when it is not, but
  `sanitizers.yml` — the only workflow that runs `ctest` — does not
  install `libxml2-utils`. So the grammar validation that caught six
  classes of real deviation may not be running on any push. Add the
  package.
- [x] **A3. `coupling*` is β = 0 only** — DONE: declared in the file
  header alongside the other limits. *Original:* The interaction tables are
  indexed by α and Tc and were swept at zero sideslip, but the buildup
  applies them at every β. Either sweep β or declare the restriction in
  the file header, as every other limit is declared.

## B. Measurements that would settle an open question

- [x] **B1. Powered vs power-off separation points, directly** — DONE
  2026-08-20, and it changed a conclusion rather than confirming one.
  125 conditions. The fan holds the separation point aft at every
  attitude from −4° to 70°, monotone in thrust at each, peaking at
  α = 30 with +0.065 chord at Tc = 8 — and it reaches the section
  tables only through local incidence, since those tables contain no
  representation of the fan. **The headline is negative**: Part II's
  threefold lift step at 16° is *not* the separation-delay signature.
  Across it the load grows ×3.38 while the delay moves ×1.09. The step
  marks the onset of limit-cycle behaviour, exactly as Part II
  suspected but could not show. Two traps recorded in
  `InductionMapExport.cpp` and in the paper: `fMin` saturates at 0 from
  α = 20 (differencing saturated values is not a measurement), and the
  shift's *ratios* are not sensitivities (f(α) is steep at the knee;
  take sensitivity from `dCZ`, which is flat at 1.50–1.65 per thrust
  doubling). At 80–90° the sign reverses, corroborated independently by
  `dCZ` reversing across the same two attitudes. Pinned by
  `TestSeparationDelay`. Data: `separation-map.json` (force-added).
  **Carried forward — now closed:** f(γ̄) ≠ mean f(γ) is bounded and
  negligible, 1e-15 to 1e-10 against f of 0.7–0.86. The cycle is wide in
  circulation and nearly stationary in the local incidence f is posed in.
  `ViscousCoupledResult` now carries per-strip incidence mean+variance so
  any consumer evaluating something nonlinear at the mean can estimate
  its own error. **But see B5**, which bounding it uncovered.
- [x] **B2. Duct separated drag at incidence** — DONE: the ring carries
  a bluff-body crossflow term on its side-projected area, CD0(90°)
  0.195 → 0.208. Only body-wake/wing interference now keeps `aeroCD0` a
  lower bound. *Original:* `aeroCD0` covers body
  friction, duct friction, and the body's slender-body crossflow branch.
  An annular ring at 90° is a bluff body the slender-body form does not
  describe, so the table is a declared lower bound at high α.
- [x] **B3. The aileron's section pitching-moment increment** — DONE:
  Δcm = −(δ/2)sinθ_h(1−cosθ_h), putting the flap load 0.21 chords aft of
  the quarter chord. Gap leakage and viscous decay remain. *Original:* The flap
  is carried as a zero-lift-angle shift, which is thin-airfoil theory's
  lift result only. The same theory supplies the moment. Gap leakage and
  viscous decay at large deflection are separate omissions; all three
  push the tabulated authority the same way, upward.
- [x] **B5. The alpha grid is a continuation path** — RESOLVED
  2026-08-21 by repeating the whole sweep at half the incidence step.
  **113 of 125 shared conditions agree to better than 1e-6**: the map is
  grid-converged over ninety percent of its extent. The twelve that are
  not are isolated bistable conditions, all post-stall, each a single
  solve settling into a different limit cycle according to the attitude it
  was reached from — α 20 (all Tc, power-off), α 24 (all Tc, power-off;
  powered too at Tc ≤ 1), α 26 Tc 2, α 30 Tc 2. Named, so they can be
  carried as an uncertainty rather than as a caveat over the whole map.

  **The question B1 left open is answered, and the answer is "sometimes".**
  Where only the power-off solve moves (α 20) its 0.011 shift passes
  through every thrust column undiminished; where powered and power-off
  move together (α 24, Tc ≤ 1) the increment cancels to 0.0005. Both
  happen in one map, so cancellation is real but must not be assumed.
  Worst case 0.027 against a 0.044 signal.

  **B1's headline is not an artifact of the step.** The 14° and 16° rows
  agree between grids to 1e-9…1e-7, and the step across them is ×3.39
  against ×1.09 on the fine grid where it was ×3.38 against ×1.09 on the
  coarse. Data: `separation-map-fine.json`; table generated into
  `figures/tables/gridconv.tex`.

- [x] **B4. Hysteresis — RESOLVED 2026-08-21, and a static table is
  adequate.** Swept 90 → −4 against the shipped −4 → 90 at β = 0, same
  speed and relaxation. Attached flow is path-independent (worst 6.6e-4 in
  CZ below α = 14). **Only 7 of 25 attitudes differ by more than 1e-3**;
  largest 1.65% of CZ at α = 24, 2.6% of Cm at α = 26.

  **Not a classical hysteresis loop**, which is the part that decides the
  model's shape. Two branches would part across the whole post-stall range
  and rejoin at its ends. Instead six of seven cluster in 18–35, the
  branches *agree* at 16 and 20 inside that band, and agree exactly from
  40 through 70 — isolated bistable conditions, the same signature B5
  produced by halving the incidence step. **The DAVE-ML tables need no
  branch axis**; the named attitudes need wider uncertainty, and 1.65% is
  well inside the ±10–20% already carried on parasite drag.

  Unexplained: α = 80 is the only departure outside the stall band and the
  largest in CZ (1.76%). It sits beside α = 90, where the fan's mechanism
  inverts (B1), but no mechanism is offered. Data:
  `models/data/aero-map-descending.json`.

  Enabled by giving `rates` its own block selector — the rate block never
  warm-starts, so it carries no branch and a reversed α list changes
  nothing in it; paying for it in a hysteresis study buys data that cannot
  answer the question.

## C. Method extensions — each wants its own branch

- [x] **C1. The section march starts from the stagnation state** — DONE
  2026-08-21/22, both halves.

  The item's premise needed correcting first. `theta = 0` is *not* the
  defect: Thwaites started at a real stagnation point, where Ue ~ a·s,
  produces θ₀ = √(0.075ν/a) unaided — which is why `MarchSurfaceRun` is
  already right without a seed. The defect is that
  `BoundaryLayerSectionModel` is a **camber line**: no thickness, so no
  stagnation region exists, `ue` is finite at station 0, and starting
  there discards the upstream run.

  **Mechanism.** `StagnationStrain` supplies the nose strain; empty
  reproduces the old march bit for bit. The seed had to enter the
  **Thwaites integral**, not the variable — Thwaites is an integral
  formula, so an assigned `theta` is overwritten at the first station. The
  equivalent start is `I₀ = θ₀²Ue₀⁶Re/0.45`. The first version assigned it
  and the test caught that it changed nothing.

  **Supplier.** `Solver/StagnationStrainTables.h`, mirroring
  `SeparationTables.h`: one Hess–Smith solve per (section, α) over ±20° at
  1°, built once. It must not run per call — this model iterates inside
  the coupling's own iteration, so a panel solve there would be paid tens
  of thousands of times per condition for a number depending only on
  (η, α). Two rules differ from the separation tables deliberately: the
  table is **signed** in α (camber makes a nose asymmetric — measured
  42.6 at +6° against 22.4 at −6°), and the edge rule **clamps** rather
  than extrapolating, because extrapolated strain goes negative and
  √(0.075ν/a) then has no real value.

  **Nothing shipped moved, and that is not luck.** `BoundaryLayerSectionModel`
  is Level-3 and appears only in its own test; every driver uses the
  anchored `PostStallSectionModel`. So C1 had no destination to be wired
  to — the tier is not in the production pipeline. Proven end to end
  instead by `TestSectionBoundaryLayer`: real CST geometry → panel solve →
  table → Thwaites → drag, cd 0.00840 → 0.01160.

  Measured en route: the seed acts mostly **through transition**. At
  Re = 5e5, a = 10 it is 1.2e-4 chords and moves cd 83%, because a thicker
  leading-edge layer raises Re_θ, trips Michel earlier, and turns a longer
  run turbulent. The magnitude check is a **vanishing-seed limit** rather
  than a bound — a bound would only be a tripwire on where transition sits.

- [~] **C2. Non-axial propulsor inflow — BOUNDED, not modelled.**
  2026-08-22. The model itself still needs once-per-revolution loading and
  that is genuinely a different solver: `BuildPropellerLattice` takes a
  scalar `axialSpeed` and bakes an axial helical wake, so a disk at
  incidence would need skewed wake legs *and* azimuthal averaging, and the
  rotor–vane path is axisymmetric end to end besides. Fabricating an
  alphaDisk sweep from either would be an artifact.

  What *was* wrong is that `alphaDiskDeg` shipped as a "validity monitor"
  with **no scale**: a consumer saw 25° and had no way to judge it. It now
  states both errors, which are not the same size. **First order, and
  correctable:** the propulsor tables are indexed by `advanceRatio` on the
  full free stream while a propeller advances on the axial component only,
  so they are read at a J high by 1/cos(alphaDisk) — 1.5% at 10°, 6.4% at
  20°, 15.5% at 30°. `advanceRatioAxial` now carries the corrected value.
  **Second, not correctable here:** the in-plane force and hub moment are
  absent entirely.

  **Decided 2026-08-22 by the user: take the better approximation.** All 26
  propulsor function references now index on `advanceRatioAxial`. No
  tabulated value changes — alphaDisk is zero at every solved condition, so
  the two agree exactly there; it changes only where an off-axis consumer
  lands.

  That switch also exposed a coverage hole worth more than the switch: **no
  staticShot exercised a propulsor table at all.** All six pinned `aero*`
  quantities, so the entire prop path — breakpoints, ordering,
  interpolation, and the variable the tables are indexed *by* — was
  unpinned, and the rebinding verified green without any check having
  looked at it. `propEncodingAtBreakpoint` closes that, at zero disk
  incidence where the map was solved. An off-axis shot cannot be generated:
  there is no solve at disk incidence to generate one from, which is C2
  itself.

  Corrected en route: I first reported the monitor did not exist. It does
  — I had checked `PropulsionMapExport.cpp`, where alphaDisk is only a
  comment, and not the assembler, which emits it.

- [x] **C3. `SolveResult::CDi` does NOT become the Trefftz value** —
  DECIDED 2026-08-21, against the change, for three independent reasons.
  Propeller thrust is `-Di` (`PanelBuilder.h`) and `CDi = Di/(qS)`, so
  changing one breaks the identity and changing both breaks the rotor. A
  rotating-frame rotor sheds a *helical* wake, which is not what a plane
  at downstream infinity models. And the coupled solver calls `Solve`
  ~1000× per condition without ever reading `CDi`, so an O(strips²) wake
  integral on every call would be paid entirely in sweeps that discard it.
  The real defect was never the default — it was that `SolveResult::CDi`
  documented itself as "induced drag" with no hint that it goes negative
  at zero lift and fits e = 1.53 on a coupled configuration. That trap is
  now documented at the point of use, pointing to `TrefftzPlane.h`, which
  stays opt-in. A consumer who reads the field now learns the trap; one
  who wants the far-field number asks for it by name.

## D. Assumptions to revisit

- [x] **D1. The rate-derivative taper** — RESOLVED 2026-08-19, and the
  answer was that the taper was never implemented AND was wrong in
  principle. Measured at two perturbation amplitudes: attached flow and
  deep stall agree to a tenth of a percent, while α ≈ 20–35 disagrees by
  34–252% with the sign not surviving. Across that band **no linear
  coefficient exists**, so taper, clamp and measured value are equally
  fabrications. `alphaRateBp` now omits those six breakpoints and a
  lookup interpolates across an acknowledged gap.
- [x] **D2. The sign of ΔC_l after roll authority collapses** —
  RESOLVED 2026-08-20. The wobble (positive 40–45°, negative 50–80°,
  positive at 90°) is **not** cycle-mean scatter: an independent
  iteration path — damped relaxation 0.02 against 0.05, both running to
  the ceiling as cycle means — reproduces every value to 0.1% and every
  sign. The cycle means are iteration-path independent for the control
  increments exactly as they are for the baseline map. Determinate
  within the model; past stall the model still rests on the anchored
  section model, so this is not yet a statement about the vehicle.

- [x] **D3. Vane names settled** — 2026-08-22, by the user: the existing
  `vaneDeflection_{Pitch,Yaw,Roll}` / `_{Bottom,Left,Top,Right}` coinage
  stands. It follows Annex A's compound pattern (the same
  `Name_Qualifier` shape as the genuine `bodyAngularRate_Roll`), and Annex
  A has no propulsive-vane concept to align to, so any alternative would
  be a different invention with no better claim.

  **Not done, and needs the purchased text:** confirming the other Annex A
  spellings against the published ANSI/AIAA S-119 standard rather than the
  reference documentation's examples. Without the standard I can only
  re-read the examples, which is the thing this item existed to stop.

## E. Solver housekeeping

- [x] **E1. A declined surface grid now says why** — DONE 2026-08-21
  (`275eb83`). The original note was that streamline tracing assumes the
  lateral fuselage surface and correctly declines the base cap. True, and
  it hid a usability defect: `SurfaceGrid::Valid()` is derived from sample
  counts, so **six distinct causes collapsed into one `false`** — and the
  base cap being declined *by design* was indistinguishable from a
  misspelled surface name. `SurfaceGridStatus` separates `NoSystem`,
  `NoSuchSurface`, `Unindexed` (the cap), `TooSmall`, `IncompleteGrid`
  and `DuplicateKey` across all six return paths, pinned by
  `TestSurfaceFlow`.

  **Still true and still conditional:** if the attachment analysis should
  ever need to cross onto the base, that wants a second patch and a join.
  Nothing needs it today.

- [x] **E2. Axial nose resolution is a `LatticeOptions` knob** — DONE
  2026-08-21. `BodyNoseRefineStations` / `BodyNoseRefineFraction`, default
  zero so no existing mesh moves. The argument was the stated one: the
  contract states shape, not mesh, and azimuthal resolution was already a
  consumer choice while axial resolution was not. Refinement only — every
  original station survives and added ones sit on the contract's own
  radius law, so a coarsening knob is deliberately absent. Pinned by
  `TestBodyNoseRefine`, which also recorded two real mesh facts: caps put
  corners on the axis, and the base cap is panelled as **concentric
  annuli**, so a single-valued-in-x radius law cannot describe its corners.

- [x] **E4. The 380x condition was the MACHINE, not the solver** —
  RESOLVED 2026-08-21, and the note it asked for is this one.

  Re-running `alpha = 21, Tc = 0.5` on an idle machine took **113 s**
  against the original 40,741 s. That alone was weak evidence, since a
  cold start takes a different continuation path (B5) and `fMean` confirmed
  the state differed (0.6737 vs 0.6787). The decisive evidence arrived
  from the environment instead: the machine runs at **254 MB free of
  8 GB** (3%), `Get-Process` itself threw `OutOfMemoryException`, and
  **two long sweeps were killed mid-run** with truncated JSON and empty
  stderr — the same silent kill that ended the first B4 attempt.

  A machine thrashing at 3% free memory produces exactly the observed
  signature: identical arithmetic, identical iterations, identical
  results, orders of magnitude of wall time, because `seconds` is
  `steady_clock` and measures wall rather than CPU. Denormal arithmetic —
  the other candidate — would have slowed all five Tc rows at that
  incidence, since they share the flow state. It slowed exactly one.

  **Operational consequence, which is the part worth keeping:** long
  sweeps on this machine are not reliable. Two of the last four were
  killed. A multi-hour sweep should be chunked, or run when memory is
  actually free, and a truncated JSON with empty stderr should be read as
  a kill rather than as a solver fault.

- [x] **E3. A malformed attachment-line call is now reported** — DONE
  2026-08-21. `ComputeAttachmentLine` takes one Weissinger row per strip;
  a multi-row chordwise lattice trips that and used to return an empty
  station list, which is byte-for-byte what a wing with no resolvable
  attachment line returns. `AttachmentLine::Status`
  (`Ok`/`TooFewStations`/`SizeMismatch`) and `Valid()` separate them,
  pinned by `TestAttachmentLine`.

  **A worse trap surfaced while documenting it**, and no status can catch
  it: passing the LE row is not just a slice. `StripLeadingEdge` steps a
  quarter of `strip.Chord` ahead of the bound segment, which is the
  leading edge only if that row spans the whole chord. On a multi-row
  stack the LE row's bound vortex sits at a quarter of *its own* panel's
  chord, so a caller who slices the LE row but keeps section-chord strips
  gets a plausible wrong answer. Stated at the contract.

## F. Papers

Two items block **every** submission and should be settled once:

- [ ] **F1. Author blocks.** Departments and AIAA member grades are
  `TODO` in all four drafts. The blocks must stay synchronised across
  Part I, Part II, SciTech and the fan-induction paper.
- [ ] **F2. Regenerate `style/aiaa-tc.cls` and `style/aiaa.bst`** with
  `latex aiaa.ins` before any real submission. The current files were
  extracted by a hand-written docstrip equivalent.

Then per paper:

- [ ] **F3. Part I** (`journal-of-aircraft`): validation beyond
  closed-form verification; the `CDi` note in §VI.E.
- [ ] **F4. Journal of Propulsion and Power**: validation anchors are the
  schedule gate for the whole paper — isolated-propeller thrust/torque
  data, and a ducted-fan-with-vanes experiment or RANS comparison. Also
  how the swirl momentum budget is presented, and the figure set.
- [ ] **F5. SciTech**: abstract deadline for the target year, the
  demonstration figure set, the validation anchor for the conference
  version, and the scope split against the JPP article.
- [~] **F6. Fan-induction paper** — operating points CONFIRMED and the
  paper corrected, 2026-08-22.

  **Internally consistent:** μ recomputed from the contract's own disk
  geometry reproduces the published 0.26 / 0.62 / 1.04 / 1.86 to within
  1%. The non-obvious part, now stated in the caption: `v_h` is referred
  to the **annular** disk area, excluding the blade root at r/R = 0.42
  (the motor hub). Recomputing with the full disk gives μ 8% low and looks
  like an error.

  **Not a transition, and the paper said it was.** Confirming these
  bracket a *real* transition needs the thrust at each speed checked
  against a trimmed condition, and **no mass is stated anywhere in the
  repo** — the contract carries shape alone. The paper's own header
  already said "the transition operating line… does not exist yet", while
  its caption called these "four points along a transition" and its body
  "the transition point". Corrected: they are a parametric sweep spanning
  the transition's μ range. Thrust does fall monotonically with airspeed,
  which is the right qualitative shape and is not a trim schedule.

  **Lip suction — SETTLED 2026-08-22 by the user: not reported
  separately, one paper carries it.** That is already the state: the
  SciTech paper reports it quantitatively as one line of its thrust split
  (+0.26 N of bore-lip suction), and nothing else reports it. The
  fan-induction paper contains no lip force at all and cannot — it models
  the fan as a bare actuator disk and panels no duct — so a scope note now
  says so and points to the single account, since a reader who sees
  "ducted fan" will otherwise wonder where the duct's own force went.

  **Still open:** the AIAA duplicate-submission position against Part II.
  The user's call.

- [ ] **F7.** Cite the JOSS paper's DOI for the software once minted.

---

## Settled — do not re-derive

Short index. Full detail lives where the work does.

| Finding | Written down in |
|---|---|
| Stagnation offset scales as √r_LE, not r_LE | Part I; `SectionPanelMethod.h` |
| Surface strain rates need a local orthonormal frame | `SurfaceFlow.h` |
| Near-field `CDi` is untrustworthy on a coupled configuration; span efficiency must use the *lifting system's* lift | `TrefftzPlane.h`, `TestTrefftzPlane` |
| Reduced-rate derivatives: the factor is 2V/ℓ, and inverting it is wrong by ~2000 | `StabilityDerivatives.h`, `TestSolverCore` |
| Camber double-count: strip frames must be the **true chord frame** | `ViscousCoupling.h` — bit this repo twice |
| Spanwise-checkerboard spurious equilibria; plain damped iteration suppresses them | `PostStallSweepExport.cpp` |
| Deep-stall limit-cycle means are iteration-path independent **at fixed continuation** — warm-start history matters more than relaxation (B5) | Part II |
| A rate derivative does **not** follow the wrench frame rule | `Solver/BodyAxes.h`, `TestBodyAxes` |
| `SolveViscousCoupled` never populates coefficient members — read the dimensional fields | `Solver/BodyAxes.h` — bit this repo twice |
| Vane mode-sum buildup is wrong (33%); per-vane summation verified (1.4%) | `models/README.md`, `PropulsionMapExport.cpp` |
| Parasite drag cannot be a constant: CD0(90°) is 17.5× the friction term | `ParasiteDragExport.cpp`; Part I, Part II |
| A hinge cannot live on one chordwise row; the flap belongs in the **section** | `ViscousCoupling.h`, `TestFlapSection` |
| Past ~30°, the aileron is predominantly a **yaw** effector | `models/README.md`; technical report |
| The fan's upstream induction needed no new physics — `DiskInduction.h` already had it | `InductionMapExport.cpp` |
| The fan's separation delay is real and thrust-ordered, but the 16° lift step is **not** its signature (load ×3.38, delay ×1.09) | Part II; `InductionMapExport.cpp` |
| A separation-point *location* is measurable where a lift *increment* across two cycle means is not | Part II; `TestSeparationDelay` |
| `fMin` saturates at 0 past α = 20 — a difference of saturated values is not a measurement | `InductionMapExport.cpp` |
| Overlap-resolved VPM is a separate project ([onurtuncer/VPM](https://github.com/onurtuncer/VPM)); five recorded rules plus the RK2 midpoint-source rule | `ParticleWake.h` |
| BEMT is a separate project; no dependency either way | `CLAUDE.md` |
| This machine runs at ~3% free memory; long sweeps get silently killed (truncated JSON, empty stderr) and wall-clock timings can be absurd | TODO E4 |
| `assert()` is compiled out — the `windows` preset is Release, so NDEBUG. Tests must use the suites' own `CHECK` macro | `TestSeparationDelay.cpp` header note |
| The base cap is panelled as concentric annuli; `RadiusAt` is single-valued in x and cannot describe its corners | `TestBodyNoseRefine` |

---

## How this was verified

Current, on the vcpkg machine (real OpenBLAS/LAPACK, nlohmann/json):

```
cmake --preset windows && cmake --build --preset windows && ctest --test-dir build/windows
```

**34/34 suites pass.** The build must run inside the MSVC dev environment
with `VCPKG_ROOT` re-set after `vcvars64.bat`, which overrides it to VS's
bundled vcpkg.

The flight model additionally validates against `models/DAVEfunc.dtd`
(DAVE-ML 2.0.1) and passes 407 checks through `models/verify-daveml.py`,
which runs as `TestDaveMLModel` — see A2 for the CI caveat.

Papers: `pdflatex && bibtex && pdflatex && pdflatex` in each paper's
folder with `TEXINPUTS=./style;`. `latexmk` does not work on this machine
(no perl). Part I 29 pp, Part II 13 pp, technical report 29 pp.
