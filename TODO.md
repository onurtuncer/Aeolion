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

- [ ] **B1. Powered vs power-off separation points, directly.** The
  fan-induction question — how far the aft fan delays separation — is
  presently answered by a lift increment differenced across two limit
  cycles, which cannot support the claim. The coupled solve already
  computes a separation point on every strip; reporting the powered and
  power-off locations outright settles it. **This is the proper close of
  the fan-induction study**, and its headline measurement.
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
- [ ] **C3. Decide whether `SolveResult::CDi` becomes the Trefftz value.**
  The far-field integral is the trustworthy one on a coupled
  configuration; the near-field number is still what the field returns,
  with the Trefftz result computed alongside. Making the far-field value
  the default is a behaviour change to a tested field and wants its own
  decision, not a drive-by.

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
  position against Part II. **B1 is this paper's headline measurement.**
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
| Deep-stall limit-cycle means are iteration-path independent | Part II |
| A rate derivative does **not** follow the wrench frame rule | `Solver/BodyAxes.h`, `TestBodyAxes` |
| `SolveViscousCoupled` never populates coefficient members — read the dimensional fields | `Solver/BodyAxes.h` — bit this repo twice |
| Vane mode-sum buildup is wrong (33%); per-vane summation verified (1.4%) | `models/README.md`, `PropulsionMapExport.cpp` |
| Parasite drag cannot be a constant: CD0(90°) is 17.5× the friction term | `ParasiteDragExport.cpp`; Part I, Part II |
| A hinge cannot live on one chordwise row; the flap belongs in the **section** | `ViscousCoupling.h`, `TestFlapSection` |
| Past ~30°, the aileron is predominantly a **yaw** effector | `models/README.md`; technical report |
| The fan's upstream induction needed no new physics — `DiskInduction.h` already had it | `InductionMapExport.cpp` |
| Overlap-resolved VPM is a separate project ([onurtuncer/VPM](https://github.com/onurtuncer/VPM)); five recorded rules plus the RK2 midpoint-source rule | `ParticleWake.h` |
| BEMT is a separate project; no dependency either way | `CLAUDE.md` |

---

## How this was verified

Current, on the vcpkg machine (real OpenBLAS/LAPACK, nlohmann/json):

```
cmake --preset windows && cmake --build --preset windows && ctest --test-dir build/windows
```

**31/31 suites pass.** The build must run inside the MSVC dev environment
with `VCPKG_ROOT` re-set after `vcvars64.bat`, which overrides it to VS's
bundled vcpkg.

The flight model additionally validates against `models/DAVEfunc.dtd`
(DAVE-ML 2.0.1) and passes 407 checks through `models/verify-daveml.py`,
which runs as `TestDaveMLModel` — see A2 for the CI caveat.

Papers: `pdflatex && bibtex && pdflatex && pdflatex` in each paper's
folder with `TEXINPUTS=./style;`. `latexmk` does not work on this machine
(no perl). Part I 29 pp, Part II 13 pp, technical report 29 pp.
