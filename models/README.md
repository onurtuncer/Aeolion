# Aetherion flight model — DAVE-ML specification

**Target artifact:** `models/AetherionFlightModel.dml` — a single DAVE-ML
2.0.1 (AIAA S-119) file carrying the airframe aerodynamics, the ducted
propulsor with its vane cruciform, and the fan-on-airframe interaction
increments. Generated, never hand-edited; committed alongside its source
JSONs (the house pattern used for paper figures).

**Status (2026-08-19):** COMPLETE end to end — five generators, assembler, and an in-repo verifier running as `TestDaveMLModel` with DTD validation. 35 tables, 407 checks, validates against DAVEfunc.dtd 2.0.1. Every table family is generated; nothing remains blocked.
Regenerate with `python models/build-daveml.py`; check with
`python models/verify-daveml.py` (or `ctest -R TestDaveMLModel`).

---

## Decisions of record (2026-08-14)

Settled in design review; do not relitigate without new evidence.

1. **Single file.** The power-on increments straddle the aero/propulsion
   boundary and DAVE-ML has no cross-file reference mechanism; one file
   keeps the shared inputs and breakpoints internally consistent and
   lets MathML calculations assemble the total wrench inside the
   exchanged model. varIDs are namespaced `aero*` / `prop*` /
   `coupling*` so a later split into two files is mechanical.
2. **No Mach axis.** The entire toolkit is incompressible — Mach is
   accepted and ignored (`ViscousCoupling.h`, `SectionBoundaryLayer.h`).
   Validity `M < 0.3` is declared in the fileHeader instead.
3. **No Reynolds axis.** All source sweeps ran at one condition
   (V = 25 m/s, Re_n ≈ 3×10⁵). Declared as a fixed reference condition,
   not pretended as an axis.
4. **Speed and RPM collapse to advance ratio `J = V/(nD)`** for the
   propulsor (inviscid similarity). Revisit only if Level-3 blade
   viscous coupling introduces a measured tip-Reynolds dependence.
5. **One antisymmetric aileron axis** (`aileronDeg`). No flaperon mode.
6. **Vane command space, not per-vane angles**: roll = common mode,
   pitch/yaw = pair differentials (mode shapes below). Mixing-matrix
   signs are pinned by checkData, never derived in prose.
7. **Part II (coupled post-stall solve) is the single table truth**,
   including the attached range; Part I's inviscid matrix is demoted to
   checkData cross-checks. One solver path, no seam.
8. **Fan state variable for the interaction tables is
   `Tc = T/(q∞·S)`**, computed inside the file from the propulsor's own
   thrust output (no algebraic loop: the prop wrench does not depend on
   Tc).
9. **Hysteresis is deferred.** Static tables carry up-branch cycle
   means; post-stall tables carry DAVE-ML uncertainty bounds from the
   cycle RMS. Representation to be revisited with the hysteresis study.
10. **Deep-stall values are cycle means** exported with their RMS as
    uncertainty bounds — the file must not launder the `Converged=false`
    honesty of the source sweeps.

---

## Governing standards

- **ANSI/AIAA S-119-2011** (approved 22 Mar 2011, reaffirmed 2016;
  doi:10.2514/4.867965) — the exchange standard. Annex A = standard
  variable names, Annex B = the DAVE-ML reference.
- **DAVE-ML Reference 2.0.1** (31 Mar 2011) / the `DAVEfunc` DTD —
  element syntax and `checkData`/`staticShot` semantics.
- **ANSI/AIAA R-004-1992** — body axes and sign conventions, adapted
  from **ISO 1151-1:1988**. The contract's FRD frame *is* R-004's
  standard aeronautical body axis system, which is exactly why the
  solver-side flip has to be explicit rather than assumed.
- **XML 1.1** and **MathML 2.0** — the two open standards DAVE-ML is
  built on.

**Variable naming — RESOLVED, no conflict.** The standard separates the
two roles: `varID` is "an internal identifier that is unique within the
file" and is *unconstrained in form*; `name` "should correspond to the
standard AIAA parameter name". So the namespaced `aero*`/`prop*`/
`coupling*` varIDs stay, and the Annex A name rides alongside:

| varID | name (Annex A) |
|---|---|
| `alphaDeg` | `angleOfAttack` |
| `betaDeg` | `angleOfSideslip` |
| `trueAirspeedMps` | `trueAirspeed` |
| `airDensityKgpm3` | `airDensity` |
| `rollRateRadps` / `pitchRateRadps` / `yawRateRadps` | `bodyAngularRate_Roll` / `_Pitch` / `_Yaw` |
| `propSpeedRevps` | `propellerSpeed` |
| `aileronDeg` | `aileronDeflection` |
| `vane{Pitch,Yaw,Roll}Deg` | no Annex A counterpart; formed to the same pattern |

**Three attributes to populate** (previously unused): `axisSystem`
("body" on every force/moment/rate), `sign` (positive-direction token —
`+UP`, `+RWD`, `TED`), and `symbol`. The `sign` attribute REFINES the
"no prose glosses" rule above: the objection is to prose standing *in
place of* a verified convention, not to a machine-readable declaration.
Populate `sign` AND pin it with a generated staticShot; a mismatch
between the two is itself a defect worth catching. `alias` is available
and deliberately unused (the reference discourages it for portability).

Validate the assembled file against the `DAVEfunc` DTD as a build step.

## Conventions (normative)

### Frames

All file quantities are **body-frame FRD** (`aetherion_body_frd`:
x forward, y right, z down), the contract's own frame. The solver frame
(x aft, y right, z up) is related by one proper rotation, 180° about y:

```
x_frd = -x_vlm    y_frd = y_vlm    z_frd = -z_vlm
```

Being a proper rotation it applies identically to force and moment
vectors: **CX, CZ, Cl, Cn flip sign between frames; CY and Cm are
invariant.** The assembler applies this map in exactly two places —
reference-point ingest and wrench output — and nowhere else. The ingest
line carries its own unit test independent of any table (this is where
the 2026-08-09 Cm_α bug lived).

### Moment reference

All moments are about the contract's `moment_reference_point`,
**(x, y, z) = (−0.2423, 0, 0) m, FRD**, stated twice in the file: prose
in the fileHeader and as constant variableDefs `XmrpM`, `YmrpM`,
`ZmrpM` so a consumer can transfer to its own CG mechanically.

### Deflection signs

Positive deflection is the **right-hand rule about the hinge-axis
vector as stated in the contract (FRD)**. English glosses (TE-up/down)
are banned from normative text; every sign is defined by citation to a
checkData staticShot. `aileronDeg` is the right-surface angle;
antisymmetry (left = −right) is applied by the generator.

### Sideslip mirroring

Power-off airframe generation runs β ≥ 0 only; the assembler writes the
full signed β range using airframe symmetry: **CY, Cl, Cn odd in β;
CX, CZ, Cm even.** Propulsor tables are parameterized by disk incidence
and crossflow azimuth instead and are never mirrored.

### Validity envelope (fileHeader)

M < 0.3; Re fixed at the reference condition; `coupling*` tables valid
for V ≥ 10 m/s (q∞ normalization degenerates toward hover — defensible
because the fan is aft of the wing and its upstream induction decays
with distance); propulsor tables valid for powered operation
(`propSpeedRevps` well above zero — the ρn²D⁴ normalization excludes
windmilling); α_disk azimuth assumption for vane tables (below).

---

## Inputs

| varID | units | meaning |
|---|---|---|
| `alphaDeg` | deg | angle of attack |
| `betaDeg` | deg | sideslip |
| `trueAirspeedMps` | m/s | V |
| `airDensityKgpm3` | kg/m³ | ρ |
| `rollRateRadps` `pitchRateRadps` `yawRateRadps` | rad/s | body rates p, q, r |
| `propSpeedRevps` | rev/s | n (RPM/60) |
| `aileronDeg` | deg | antisymmetric aileron, right surface |
| `vanePitchDeg` `vaneYawDeg` `vaneRollDeg` | deg | vane mode commands |

## Constants (variableDefs, values filled by assembler from the contract)

`WingAreaM2` (≈0.1883), `WingSpanM` (1.0629), `WingChordM` (0.1771),
`DiskDiameterM` (0.203), `XmrpM` (−0.2423), `YmrpM` (0), `ZmrpM` (0).

**Drag accounting — read before touching CX.** Three components, and
only two are computed anywhere:

| Component | Where it lives | Status |
|---|---|---|
| induced | `aeroCX`/`aeroCZ` (near-field, in `Base.Di`) | computed |
| wing section profile | the same tables — the section model's `cd` is integrated per strip into `Base.Di` (`PostStallSweepExport.cpp:390`, "induced + profile") | computed |
| body + duct parasite | nothing | **MISSING** |

So the body-axis force tables already carry induced *and* wing profile
drag; they are not CDi tables. What they omit is the body/duct parasite
term. Two consequences:

1. **Never add a whole-airframe CD0 buildup on top of these tables** — it
   would double-count the wing's skin friction, which is already inside
   the section `cd`. The parasite term must be **body + duct wetted area
   only**.
2. **A constant CD0 is wrong for this model's α range.** The tables run
   to α = 90°, where a body's crossflow drag dominates and is strongly
   attitude-dependent (a slender body broadside carries a crossflow drag
   of order 1.2 on projected area). A single constant would underpredict
   axial force in exactly the deep-stall regime the model exists to
   cover.

Therefore parasite drag is a **table, `aeroCD0(alphaDeg)`**, not a
constant, generated by a new path: `DragEstimate`'s component buildup
(body + duct wetted areas from the contract) for the attached branch,
blended into a crossflow-drag estimate on the body's projected area at
high incidence. Neither half is wired today — `DragEstimate` has never
been called, and no crossflow model exists. Until that driver is built
the table ships as an all-zero placeholder with the omission declared in
the fileHeader, so a consumer sees a stated gap rather than a silently
optimistic drag polar.

## Derived variables (MathML calculations)

- `qbarPa = ½ ρ V²`
- `J = V / (n · DiskDiameterM)`
- `alphaDiskDeg`: freestream in FRD is
  u = (cosα cosβ, sinβ, sinα cosβ); rotation axis is +x, so
  `cos(alphaDisk) = cosα·cosβ`. **Validity monitor with a stated scale.**
  Two errors grow with it and they are not the same size. *First order,
  correctable:* the propulsor tables are indexed by `advanceRatio` on the
  FULL free stream, while a propeller advances on the axial component
  only, so they are read at a J high by 1/cos(alphaDisk) — 1.5% at 10°,
  6.4% at 20°, 15.5% at 30°. *Second, not correctable here:* the in-plane
  force and hub moment of a disk at incidence are absent entirely, since
  the generating solver is axisymmetric end to end.
- `advanceRatioAxial = advanceRatio · cos(alphaDisk)` — **the variable the
  propulsor tables are indexed by** (decision of record, 2026-08-22). A
  propeller advances on the component of the free stream along its own
  axis, so `V/(nD)` is the wrong argument off-axis, high by 1/cos(alphaDisk).
  At the conditions the tables were generated at alphaDisk is zero and the
  two are identical, so **no tabulated value depends on this** — it changes
  only where an off-axis consumer lands in them. `advanceRatio` remains
  defined and output, as the geometric advance ratio.
- `phiWDeg = atan2(sinβ, sinα·cosβ)` — crossflow azimuth, 0 in the
  pitch plane
- `phat = p·b/(2V)`, `qhat = q·c̄/(2V)`, `rhat = r·b/(2V)` (the
  StabilityDerivatives `_nd` convention)
- `Tc = T / (qbarPa · WingAreaM2)` with T from the propulsor's own
  thrust output — feeds the `coupling*` tables
- Total wrench buildup: aero terms dimensionalized by q∞S(b,c̄), prop
  terms by ρn²D⁴(D⁵), prop disk-frame components rotated into FRD by
  `phiWDeg`, coupling increments added; outputs both the coefficient
  sets and dimensional `FXTotalN … MZTotalNm`.

## Vane mode shapes

Contract vane order and outward radial hinge axes: bottom [0,0,1],
left [0,−1,0], top [0,0,−1], right [0,1,0].

| Mode | (δ_bottom, δ_left, δ_top, δ_right) | Net effect |
|---|---|---|
| `vaneRollDeg` | (+δ, +δ, +δ, +δ) common mode | pure Mx couple |
| `vanePitchDeg` | (0, −δ, 0, +δ) y-pair differential | Fz + My |
| `vaneYawDeg` | (+δ, 0, −δ, 0) z-pair differential | Fy + Mz |

Numeric signs of the ± entries are pinned by the vane-mode checkData
shots. Allocation saturates per-vane sums against the contract's
**soft limit ±15°** (hard stops ±20° are not for normal operation).

---

## Breakpoint sets

| bpID | values |
|---|---|
| `alphaBp` (deg) | −4, −2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 35, 40, 45, 50, 60, 70, 80, 90 |
| `betaBp` (deg) | ±30, ±20, ±15, ±10, ±5, 0 (generated β ≥ 0, mirrored) |
| `aileronBp` (deg) | −20, −10, −5, 0, 5, 10, 20 |
| `jBp` | 0, 0.1, … , 1.0 |
| `alphaDiskBp` (deg) | 0, 15, 30, 45, 60, 75, 90 |
| `vaneBp` (deg) | −15, −10, −5, 0, 5, 10, 15 |
| `tcBp` | 0, 0.5, 1, 2, 4, 8 |

α is dense (2°) through the stall break — CLmax is emergent at 18° —
and coarsens in the plate regime. Refine from generated data if the
verifier's interpolation-error spot-checks demand it.

---

## Table inventory

### `aero*` — airframe, power-off (source: Part II coupled solve)

| Tables | Axes | Notes |
|---|---|---|
| `aeroCX aeroCY aeroCZ aeroCl aeroCm aeroCn` | α × β | baseline, aileron 0; **carry induced + wing profile drag**, not CDi alone; post-stall entries carry uncertainty bounds from cycle RMS |
| `aeroCD0` | α | body + duct parasite only (never the wing — see the drag-accounting note above). BLOCKED: placeholder zeros until the buildup + crossflow driver exists |
| `aeroClp aeroClr aeroCmq aeroCnp aeroCnr` (+ any further `_nd` fields StabilityDerivatives provides) | α | β = 0; **tapered linearly to zero over α ∈ [20°, 40°]**, held zero beyond — a documented assumption, revisit with the hysteresis work |
| `aeroDCl aeroDCn aeroDCm aeroDCX` | α × aileron | increments from baseline, β = 0; the coupled solve's Kirchhoff attenuation provides the post-stall effectiveness decay — verify it decays plausibly before shipping, clamp only if it does not |

### `prop*` — rotor + duct + vanes (source: SciTech configuration, rotor-vane coupled solve)

Body-axis (FRD) components directly — the vanes are body-fixed and the
inflow is axial, so there is no disk-frame rotation to apply:

| Tables | Axes | Notes |
|---|---|---|
| `propCT propCQ` | J | vanes neutral; thrust = −Di, shaft torque = −Mx of the rotating-frame solve, converted per the frame rule |
| `propDCXvane … propDCnvane` + `propDCQvane` | J × δ_i | **PER-VANE, not per-mode** — one vane's response to its own angle, serving all four positions by rotational symmetry. Full 6-component wrench + torque back-pressure; ±δ is not antisymmetric (jet dynamic pressure), so no mirroring in δ |

**MEASURED (2026-08-15), and it changed the structure.** The obvious
buildup — tabulate three command modes, add them — was tested against
directly-solved simultaneous commands and **fails**:

| Pair | worst error |
|---|---|
| pitch+yaw | 0.09 % — superposes |
| pitch+roll, yaw+roll | 6–33 % — does **not** superpose |

Pitch and yaw act on *disjoint* vane pairs, so they add. Roll is the
common mode and re-deflects the *same* vanes, and a vane's load is
nonlinear in its own angle. The pitch+yaw result (all four vanes
deflected, still 0.05 %) proves vane-to-vane interference is
negligible, so the failure is purely per-vane nonlinearity.

Replacement, **verified** against every command in the study:

| Buildup | worst error | tables |
|---|---|---|
| sum of 3 command modes | 32.7 % | 21 |
| sum of 4 per-vane responses | **1.4 %** | **7** |

So the buildup is `ΔC = Σ_i propDC*vane(J, δ_i)` with `δ_i` the per-vane
angles from the mixing matrix — 20× more accurate and 3× cheaper.
Reproduce: `aeolion_propulsion_map <handoff> <out> <rpm> combo 0.4 3`
and `... single 0.4 3`.

**Non-convergence, also measured:** positive roll at J = 0.6 does not
converge, and neither cure works — 24 vs 48 outer passes give identical
residuals (limit cycle, not slow convergence), and damping relaxation
0.35 → 0.15 still fails. Sign-asymmetric (negative roll converges in
2–3 passes). Two iteration paths bound the value: thrust reproducible to
0.004 % / 0.06 % at δ = +5 / +15, but 3.5 % apart at +10. Carry as
DAVE-ML uncertainty bounds.

### `coupling*` — fan-on-airframe increments (BLOCKED, see Open items)

| Tables | Axes | Notes |
|---|---|---|
| `couplingDCX … couplingDCn` | α × Tc | β = 0; valid V ≥ 10 m/s; Δ from the power-off baseline at matched α |

---

## checkData (staticShots — every one generated from a solve, none hand-written)

1. **Moment-transfer pin**: one condition exported about the contract
   reference point and about a point shifted +0.10 m in x_frd; the Cm
   difference must equal the hand-derivable transfer arm. A wrong frame
   flip errors by 2·0.2423 m ≈ two chord lengths.
2. **Aileron sign pin**: `aileronDeg = +5` at α = 4°; the solved Cl is
   the sign definition, cited by the `aileronDeg` variableDef.
3. **Rate-derivative pin**: Cl_p = −0.452 (the `TestSolverCore` value).
4. **Vane-mode pins**: pure pitch, pure yaw, pure roll commands — each
   pinning the dominant components and recording the off-axis residuals
   at their solved values (they are physics, not zeros). These pin the
   mixing signs.
5. **Propulsor pin**: the SciTech Table 1 reference row (the
   `TestVaneCascade` configuration; its thrust/pass count is already
   the paper's own consistency anchor).
6. **Part I cross-checks**: attached-range CL/Cm at a few (α, β) points
   from the Part I matrix, with tolerance loose enough for the
   viscous-vs-inviscid delta — the recorded justification for
   decision 7.
7. **Summed-wrench pin**: one full-up condition (α, β, rates, aileron,
   vane commands, powered) through the MathML buildup to dimensional
   `FXTotalN … MZTotalNm` — pins normalization, rotation, and summation
   end to end.

---

## Generation pipeline

```
app/AeroMapExport.cpp        -> models/data/aero-map.json        (baseline + rates + aileron)
app/PropulsionMapExport.cpp  -> models/data/propulsion-map.json  (prop baseline + vane modes)
app/InductionMapExport.cpp   -> models/data/coupling-map.json    (BLOCKED on upstream-induction model)
models/build-daveml.py       -> models/AetherionFlightModel.dml  (assembler: cached JSONs -> XML + checkData)
models/verify-daveml.py      -> ctest suite: DTD-validate, evaluate every staticShot with an
                                independent gridded-table/MathML interpreter, compare within tolerance
```

Each driver exports JSON independently (the `PostStallSweepExport`
precedent), so regenerating the vane map never reruns the post-stall
map. The `.dml` and the JSONs are committed together; the verifier runs
in CI so the file cannot drift from its generators. No external DAVE-ML
dependency (no Janus) — the verifier's evaluator is deliberately small
and in-repo.

Rough solve budgets: aero baseline 25α × 6β = 150 coupled solves;
aileron 25α × 6δ = 150; rates: central differences at the attached α
points; prop baseline 11J × 7α_disk = 77 rotor-vane coupled solves;
vane modes 11J × 7α_disk × 6δ × 3 ≈ 1400 (coarsen α_disk to
{0, 45, 90} for the first pass ≈ 600 if runtime bites — record the
coarsening here if taken).

---

## Declared limits inherited from how the sweeps are run

Two properties of the generation, not of the aerodynamics. Both apply to
every post-stall row of `aero*` and `coupling*`, because both sweeps are
built the same way.

**The sweeps are continuations, so the path is part of the condition.**
Every solve above the stall is warm-started from the preceding incidence
(`ViscousCouplingOptions::InitialGamma`). Which attitudes were visited
therefore helps determine which limit cycle a condition settles into --- it
is not merely the order in which the tables were filled. Measured
2026-08-20: arriving at alpha = 20, Tc = 0.5 from 16 rather than 18 moves
the span-mean separation point by 0.013, against a fan effect of 0.018 at
the same condition; at alpha = 26, Tc = 1 a coarser path reverses the sign
of the measured delay outright. Most conditions are insensitive --- every
row at alpha = 16 reproduces to four digits --- and the finer continuation
is the better approximation, so the shipped tables are the better-resolved
of the two compared.

That check is now done. Repeating the sweep at half the incidence step,
**113 of 125 shared conditions agree to better than 1e-6** — the tables
are grid-converged over ninety percent of their extent, including every
condition the fan-induction comparison rests on. The twelve that disagree
are isolated bistable conditions, all post-stall, each a single solve
settling into a different limit cycle according to the attitude it was
reached from: alpha 20 (all Tc, power-off moves), alpha 24 (all Tc,
power-off; powered too at Tc <= 1), alpha 26 Tc 2, alpha 30 Tc 2. Where
only the power-off solve moves its shift passes through every thrust
column undiminished; where both move together the increment nearly
cancels. Cancellation is therefore real but not reliable. Those rows carry
an uncertainty of order 0.01 in separation location, 0.027 at worst. This also narrows the earlier finding
that limit-cycle means are iteration-path independent: that holds at
*fixed continuation*, and warm-start history matters considerably more
than the relaxation does.

**The tables are the ascending-alpha branch, and that is adequate.** Every
sweep is a continuation started below the stall, so the tables are one
branch of a system that could in principle have two. Sweeping the same grid
downward from 90 degrees settles the cost. Below alpha 14 the branches
agree to 6.6e-4 in CZ — attached flow is path-independent. Only **7 of 25
attitudes differ by more than 1e-3**, the largest being 1.65% of CZ at
alpha 24 and 2.6% of Cm at alpha 26.

This is **not** a classical hysteresis loop, and the distinction decides
how the model is built. A two-valued system would show the branches parting
across the whole post-stall range and rejoining at its ends. Instead six of
the seven cluster in 18–35, the branches agree at 16 and 20 *inside* that
band, and agree exactly from 40 through 70 — isolated bistable conditions,
the same signature the incidence-step refinement produced. So the tables
need **no branch axis**; the named attitudes need wider uncertainty, and
1.65% sits well inside the ±10–20% already carried on the parasite buildup.
One departure of 1.76% in CZ at alpha 80 falls outside the stall band and is
recorded without an explanation. Data: `aero-map-descending.json`.

**Cycle-mean nonlinearity is present and negligible.** Anything reported
past stall is evaluated on a final sweep at the cycle-mean circulation, so
for nonlinear `g`, `g(gammabar)` is not `mean g(gamma)`. The leading term
is second order in the cycle width, and `ViscousCoupledResult` now carries
the per-strip incidence mean and variance needed to evaluate it. For the
separation point it lands at 1e-15 to 1e-10 against values of order unity:
the cycle is wide in circulation and nearly stationary in local incidence.
Recorded because it is the approximation a reader would reasonably suspect,
and it is not the one that matters --- the continuation is.

## Open items

- [x] ~~Reconcile varIDs against S-119 Annex A~~ — RESOLVED 2026-08-15:
      no conflict, `varID` is unconstrained and `name` carries the
      standard name (mapping table above). Remaining: confirm the
      spellings against the published standard text rather than the
      reference docs' examples, and settle the vane names.
- [x] ~~Validate against the `DAVEfunc` DTD in CI~~ — DONE 2026-08-19; `verify-daveml.py` runs xmllint when present. It found six classes of real deviation on first run (see the technical report).
- [x] ~~`aeroCD0` parasite-drag driver~~ — DONE 2026-08-15,
      `app/ParasiteDragExport.cpp`. Body+duct wetted areas from the
      contract through `DragEstimate` (previously never called by
      anything) plus an Allen-Perkins/Jorgensen crossflow branch.
      MEASURED: friction CD0 = 0.01058 (body 0.00521 ≈ duct 0.00507 —
      the duct's shorter chord raises its Cf), crossflow coefficient
      0.1848, so **CD0(90°) = 0.1954, 17.5× the friction term**. That
      settles the constant-vs-table question: a constant would omit 95%
      of parasite drag at 90°. Excludes the wing by construction.
      Reflected in Part I §Parasite drag and Part II §The parasite
      branch.
- [x] ~~`coupling*` tables blocked on the upstream-induction model~~ —
      DONE 2026-08-19. The blocker was **STALE**: `DiskInduction.h`
      already implemented the semi-infinite vortex cylinder (pinned by
      `TestDiskInduction`, already used by the attachment sweep), so no
      new physics was needed — only `app/InductionMapExport.cpp`.
      MEASURED, 125 conditions (25 α × 5 Tc): below α=14 both solves
      converge and the fan adds **2–4% lift**, growing smoothly with Tc —
      solid. From 16 the increment jumps to **7–9%, peaking at α=20**,
      just past CLmax(18), then decays and **reverses sign by 90°**
      (plate broadside: the axial induction is perpendicular to the
      freestream, so it reduces effective incidence).
      CAVEAT stated in all three documents: the step at 16 coincides
      exactly with limit-cycle onset in BOTH solves (residuals 1e-4 →
      0.2–0.6), so post-stall increments are differences of two cycle
      means and the step size is confounded. What survives is the
      ORDERING — clean monotone Tc scaling (ratios 1.49/1.49/1.49 at
      α=20), which scatter would not produce. Real thrust-ordered effect;
      attribution to separation delay is consistent but NOT established.
      SIMILARITY justifying the single-speed sweep: 4(A/S)(1+vi/V)(vi/V)
      = Tc, so vi/V depends on Tc alone and the coefficient increments
      are speed-independent. Induced velocities satisfy momentum theory
      to every printed digit.
- [ ] **Report powered vs power-off SEPARATION POINTS directly** — the
      coupled solve already computes f per strip, so the fan-induction
      question can be answered outright instead of inferred from a lift
      increment taken across two limit cycles. The proper close of 3b.
- [x] ~~`aeroDC*` (aileron) tables~~ — DONE 2026-08-17. Two fixes, in
      order. FIRST the blocker: a hinge cannot live on one chordwise row
      (`MinRowsToResolveHinge = 2` vs the coupling's one-row-per-strip
      contract), so a deflection was EXACTLY ZERO at every attitude,
      silently. Fixed by carrying the flap in the SECTION —
      `Geometry::FlapEffectiveness` (thin-airfoil tau) +
      `StripSection::{FlapChordFraction, FlapDeflectionDeg,
      EffectiveAlpha0Deg()}`. Additive, so existing consumers are
      bit-identical. THEN the sweep: `aeolion_aero_map ... aileron`
      writes `models/data/aero-aileron.json`, 100 conditions over 25
      alphas.
      VALIDATED at alpha=0 against the RESOLVED-HINGE 8-row lattice to
      **1.4%** — an independent geometric representation of the same
      flap. MIRROR EXACT: the negative-deflection probe reproduces the
      mirrored positive row at all 25 alphas, odd components to 0%, even
      to 2.6e-5%, so the one-sided sweep costs nothing.
      PHYSICS WORTH KNOWING: roll authority falls to 72% by alpha=20 and
      under 10% by 30, while |dCn/dCl| climbs 0.09 → above 1. **Past
      ~30 deg the aileron is predominantly a YAW effector** — the classic
      pre-departure signature. An inviscid-sourced table would have shown
      gentle decay and small adverse yaw: an aircraft that rolls
      obediently out of a stall.
      tau is strongly concave — the 12%-chord aileron is worth **0.432**,
      not 0.12 (linear guess understates 3.5x). `TestFlapSection` (31st
      suite) pins the exact limits, textbook half/quarter-chord values,
      additive inertness, and the stall decay.
      CAVEATS: deep-stall rows are limit-cycle means; once |dCl| collapses
      its sign changes more than once (positive 40-45, negative 50-80,
      positive at 90). RESOLVED 2026-08-20 (D2): an independent
      iteration path (relaxation 0.02 vs 0.05, both cycle means)
      reproduces every value to 0.1% and every sign, so the wobble is
      determinate rather than scatter. It is a property of the MODEL,
      not yet of the vehicle. The flap
      model is LIFT-ONLY: no section cm increment, no gap leakage, no
      viscous decay at large deflection, so tabulated authority is an
      upper bound.
- [x] ~~Simultaneous-command superposition~~ — MEASURED 2026-08-15;
      mode-sum rejected, per-vane summation adopted and verified (above).
- [x] ~~Vane-mode azimuth error at φ_w = 45°~~ — moot: the propulsor
      model is axial-inflow only (decision below), so there is no
      azimuth to sweep. Non-axial modelling is the real open item.
- [ ] Non-axial propulsor inflow: the rotor–vane machinery is
      axisymmetric end to end (slipstream bands, azimuthal-mean vane
      feedback, `AxialInflowFromBands`), so disk incidence is not
      representable. `alphaDiskDeg` is exported as a validity monitor,
      never as a table axis.
- [ ] Duct separated drag at incidence — the one parasite contribution
      the slender-body crossflow form cannot supply; `aeroCD0` is a
      lower bound at high α until it lands.
- [ ] Rate-derivative taper [20°, 40°]: assumption, revisit with the
      hysteresis study.
- [ ] Hysteresis representation (deferred by decision 9).
