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

- [~] **A1. Uncertainty bounds — plumbing DONE, data pending.**
  The solver now measures the cycle WIDTH (`CycleFluctuation()`), the
  drivers export it, and the assembler emits it as a DAVE-ML
  `uncertainty`/`normalPDF` with per-cell additive bounds — verified end
  to end against the DTD on synthesised data before committing to the
  sweeps. The bound is deliberately the cycle width and NOT the residual:
  the residual is a max over strips of a section mismatch and sits at
  0.2–0.6 post-stall, while the load-level fluctuation is ~1e-6, so
  quoting the residual would overstate uncertainty by five orders of
  magnitude. Remaining: the sweeps that populate it (running).
  *Original statement:*
  `models/README.md` decisions 9 and 10 say post-stall tables carry
  DAVE-ML uncertainty bounds derived from the cycle RMS, and the
  technical report repeats it. `build-daveml.py` emits none — the
  standard's `uncertainty` element appears nowhere in the file. Either
  emit them or strike the claim from both documents. Emitting is the
  better answer: 110 of 175 airframe conditions and 72 of 125
  interaction conditions are cycle means, and a consumer currently
  cannot tell which. The sweeps already export per-condition residuals;
  the cycle RMS itself would have to be collected.
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
- [ ] **B5. The alpha grid is a continuation path, and the path
  sensitivity is the size of the signal.** Found 2026-08-20 while closing
  B1's carried-forward item. Every solve warm-starts from the preceding
  incidence, so *which* attitudes are run is part of the specification of
  a post-stall condition, not just the order of computation. Arriving at
  α = 20, Tc = 0.5 from 16 rather than 18 moves fMean 0.7189 → 0.7061:
  both in 1000-iteration cycles, differing only in history, by **0.013
  against a fan effect of 0.018 at that condition**. Two consequences.
  (a) B1's increments difference powered against power-off at a common
  attitude on a common grid, so they are internally consistent — but
  whether the path dependence *cancels* in the difference is **untested**,
  and settling it needs a second full sweep on a different alpha grid
  (~2 h). (b) It narrows D2: cycle means are iteration-path independent
  *at fixed continuation*; warm-start history matters more than the
  relaxation does. Documented on the alpha selector in
  `InductionMapExport.cpp` so the next subset run does not repeat it.

- [ ] **B4. Hysteresis.** The maps are the ascending-α branch by
  construction (warm-start continuation up each column). Whether the
  descending branch differs materially decides whether a static gridded
  table can represent this vehicle near stall at all.

## C. Method extensions — each wants its own branch

- [ ] **C1. Wire the boundary-layer march to the real attachment point.**
  The item the attachment-line work stopped short of on purpose.
  Everything a march needs is produced and nothing consumes it:
  `SectionSolution::UpperRun()`/`LowerRun()` give U_e(s) from the
  stagnation point, `StagnationMomentumThickness` gives θ₀,
  `SurfaceStreamline` gives U_e(s) and h(s) on the body, and
  `AttachmentStation` gives Rbar so the march knows whether it may start
  laminar. `SectionBoundaryLayer.h` still starts at the camber-line
  leading edge with θ = 0 — exactly the approximation the attachment
  work removes. A behaviour change to a tested module.
- [ ] **C2. Non-axial propulsor inflow.** The rotor–vane machinery is
  axisymmetric end to end (slipstream bands, azimuthal-mean vane
  feedback, `AxialInflowFromBands`), so disk incidence is not
  representable and `alphaDiskDeg` is exported as a validity monitor
  rather than faked as a table axis. Representing it needs
  once-per-revolution loading, which is a different solver.
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

- [ ] **D3. Confirm the S-119 Annex A spellings** against the published
  standard text rather than the reference documentation's examples, and
  settle the vane command names, which have no Annex A counterpart.

## E. Solver housekeeping

- [ ] **E1.** Body streamline tracing assumes the lateral fuselage
  surface; the base cap is correctly declined. If the attachment
  analysis should ever cross onto the base, that needs a second patch
  and a join.
- [ ] **E2.** `RefineNoseStations` is a consumer-side workaround in
  `AttachmentSweepExport.cpp` for the stagnation point falling inside
  the first panel ring. Whether `BuildBody` should offer it as a
  `LatticeOptions` knob is open — the argument for the knob is that
  axial spacing is exactly as much a consumer choice as
  `BodyCircumferentialPanels` already is.
- [ ] **E4. One condition costs 380x normal, with no outer-loop symptom.**
  Found 2026-08-21 during B5's fine-grid sweep. `alpha = 21, Tc = 0.5`
  took **40,741 s (11.3 h)**; every other row at the same incidence took
  107-163 s. The output is *normal in every respect* — same 1000
  iterations as its neighbours, residual 0.6307 against 0.5527 at 20 and
  0.6627 at 22, `fMean` 0.6787 interpolating smoothly between them, cycle
  fluctuation 5.45e-08. Same iteration count and same quality for 380x
  the wall time, so the cost is **inside** the iterations, not in more of
  them: almost certainly the inner section root-find burning its full
  budget on some strip every outer sweep where neighbouring conditions
  resolve it in a few steps.

  Two reasons this matters more than a slow row. It is **invisible to the
  shipped 2-degree map**, which skips alpha = 21 entirely — so it was
  never going to be found without a finer sweep. And because there is no
  outer-loop signal, a consumer who hits it sees an apparent *hang*, not a
  diagnosis; `Converged`, `Iterations` and the residual all look ordinary.

  Wants an inner-iteration budget that is *counted and reported* — a
  per-solve tally of inner steps alongside `Iterations`, so a pathological
  condition announces itself instead of merely taking a long time. Not a
  correctness defect: the row's numbers are sound and B5's comparison is
  unaffected.

- [ ] **E3.** `AttachmentLine` takes one Weissinger row per strip. A
  multi-row chordwise lattice needs the caller to pass the leading-edge
  row.

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
- [ ] **F6. Fan-induction paper**: confirm the fan operating points that
  bracket a real transition, decide whether the duct's own lip suction
  is reported separately, and check the AIAA duplicate-submission
  position against Part II. **B1 is done** — `separation-map.json` is this paper's headline data, and the existing Δα = +0.80° horizontal-shift result corroborates it from an independent driver.
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

---

## How this was verified

Current, on the vcpkg machine (real OpenBLAS/LAPACK, nlohmann/json):

```
cmake --preset windows && cmake --build --preset windows && ctest --test-dir build/windows
```

**33/33 suites pass.** The build must run inside the MSVC dev environment
with `VCPKG_ROOT` re-set after `vcvars64.bat`, which overrides it to VS's
bundled vcpkg.

The flight model additionally validates against `models/DAVEfunc.dtd`
(DAVE-ML 2.0.1) and passes 407 checks through `models/verify-daveml.py`,
which runs as `TestDaveMLModel` — see A2 for the CI caveat.

Papers: `pdflatex && bibtex && pdflatex && pdflatex` in each paper's
folder with `TEXINPUTS=./style;`. `latexmk` does not work on this machine
(no perl). Part I 29 pp, Part II 13 pp, technical report 29 pp.
