# Journal of Aircraft submission — aft-fan inflow induction

**Venue:** AIAA Journal of Aircraft
**Type:** full research article
**Status:** skeleton draft written (`paper.tex`, 6 pages, compiles). The
upstream induction model (`Solver/DiskInduction.h`, verified by
`TestDiskInduction`) exists and is wired into `aeolion_attachment_sweep`,
which now takes optional thrust and airspeed. The headline number is
**Δα ≈ 0.7–0.9° of incidence** bought by the fan at the transition point
V = 12 m/s, T = 25 N — see "The onset that did not move" below, which is
also a cautionary tale about how that number was nearly missed.

Still to build: the two-way coupling (the rotor does not yet see the
airframe) and the Level-B blade-lattice induction that would verify the
Level-A cylinder used here.

**What the draft does and does not contain.** Sections I--IV are written
against implemented, tested code and against the sweeps recorded below;
every quantitative claim traces to one of them. Section II now carries the
full derivation of the vortex-cylinder model --- the disk/cylinder
equivalence, the ring integral and its closed form, the three momentum
limits recovered rather than assumed, the inverse-square upstream decay,
and the annulus by superposition --- with three figures and two generated
tables. Section V (Conclusions) is deliberately a stub, because it needs
the two-way coupling and the transition operating line, and writing it
from the one-way numbers would present a lower bound as a result.

**Regenerating the figures and tables:**

```
aeolion_fan_induction tests/Data/AeolionGeometryHandoff-1.8.0.json     papers/journal-of-aircraft-fan-induction/figures/fan-induction.json
# the two sweeps the separation figure compares (V = 12 m/s):
aeolion_attachment_sweep tests/Data/AeolionGeometryHandoff-1.8.0.json     papers/journal-of-aircraft-fan-induction/figures/fine-off.json 0.001 12.0 3.0 16.0 0.5
aeolion_attachment_sweep tests/Data/AeolionGeometryHandoff-1.8.0.json     papers/journal-of-aircraft-fan-induction/figures/fine-on.json  25.0  12.0 3.0 16.0 0.5
cd papers/journal-of-aircraft-fan-induction/figures && python render-fan-figures.py
```

The `.json` files are untracked (the repo ignores `*.json`); the figures,
tables and renderer are tracked, the same convention the other papers
follow.

The second Journal of Aircraft article
([../journal-of-aircraft/](../journal-of-aircraft/)) establishes where a
coupled wing–body solve attaches and where it separates. This one spends
that machinery on the question the JETTAIL airframe actually poses.

## The question

**How far does an aft-mounted ducted fan delay wing separation during
tail-sitter transition, and where on the span?**

## Why the mechanism is suction, not propwash

This was very nearly scoped as a blown-wing study, which would have been
wrong for this airframe. On the 1.8.0 handoff the fan is **behind** the
wing:

| | solver x [m] |
|---|---|
| wing LE → TE | 0.207 → 0.384 |
| duct LE → TE | 0.422 → 0.514 |

The duct leading edge sits 38 mm — 0.21 chord — aft of the wing trailing
edge. Its slipstream goes downstream, away from the wing. There is no
propwash over the wing at all.

What the fan does do is induce flow *upstream* of its disk. Actuator-disk
theory puts the axial induction at zero far upstream, rising to `vi` at
the disk; the wing sits inside that accelerating field. Acceleration is a
**favourable pressure gradient**, and a favourable gradient delays
separation. Same conclusion a blown wing would give, opposite mechanism,
and a good deal less studied.

**Note for anyone reaching for the existing slipstream code:**
`Solver::SlipstreamField` returns zero induced velocity for `point.x < 0`
by construction — it is a momentum-theory *wake* reconstruction and models
nothing upstream of the disk. Used as-is here it would report exactly zero
effect on the wing, which is a modelling artifact, not a result.

## Why this airframe makes the point sharply

The alignment is the reason this is worth a paper rather than a note:

- duct outer radius / semi-span = 0.1125 / 0.5314 = **0.21**
- the second JoA paper's worst separation station = **|2y/b| ≈ 0.15–0.23**

The fan sits directly over the span station that separates first. The part
of the wing most prone to letting go is the part the fan protects, and
that is a configuration statement, not a coincidence of this particular
mesh.

## How big is the effect likely to be?

Worth establishing before committing, because a second-order effect would
not carry an article. The annular disk (bore radius 0.1045 m less the
0.0406 m tail boom) is `A = 0.0291 m^2`, so the hover induced velocity
`v_h = sqrt(T / (2 rho A))` runs:

| thrust [N] | `v_h` [m/s] |
|---|---|
| 15 | 14.5 |
| 25 | 18.7 |
| 40 | 23.7 |

against the 25 m/s cruise the second paper solves at. So `v_h` is
**comparable to the flight speed across the whole of transition**
(`mu = Vinf/v_h` runs 0 to roughly 1.5), and the induction is a
first-order term rather than a correction. It is also strongly radial: the
disk bore is 0.10 m against a 0.53 m semi-span, so the induction is
concentrated exactly where the separation is — inboard — and is negligible
outboard. That radial falloff is not a nuisance, it is the result.

## First numbers from the induction model (2026-08-09)

`Solver/DiskInduction.h` now exists and is verified (see
`TestDiskInduction`). Probing it on the real geometry — annular disk at
solver x = 0.468, bore 0.1045 m, boom 0.0406 m, wing from x = 0.207 to
0.384 — gives the following **increase in axial velocity across the wing
chord**, as a percentage of freestream:

| 2y/b | hover-ish (V=5) | transition (V=12) | late (V=18) | cruise (V=25) |
|---|---|---|---|---|
| 0.10 | +66% | +23.0% | +11.2% | +4.5% |
| 0.15 | +57% | +19.7% | +9.6% | +3.8% |
| 0.20 | +42% | +14.5% | +7.1% | +2.8% |
| 0.30 | +15% | +5.1% | +2.5% | +1.0% |
| 0.60 | −1.8% | −0.6% | −0.3% | −0.1% |
| 1.00 | −1.1% | −0.4% | −0.2% | −0.1% |

Three things fall out, and together they are the paper:

1. **The gradient, not the magnitude, is the mechanism.** At the quarter
   chord the induction is a modest few percent of freestream. Across the
   *chord* it is 15–23% in transition, because the trailing edge sits
   three times closer to the disk than the leading edge. That is a
   favourable streamwise gradient, and it is concentrated at the trailing
   edge — which is exactly where turbulent separation begins and from
   which it marches forward.

2. **It is concentrated inboard**, falling by a factor of four between
   2y/b = 0.10 and 0.30. The second paper found separation most advanced
   at |2y/b| ≈ 0.15–0.23. The fan protects the span it sheds first.

3. **It REVERSES outboard.** Beyond roughly 2y/b ≈ 0.4 the chordwise
   change goes slightly negative: outside the vortex cylinder the induced
   axial velocity is a return flow, so the fan marginally *penalizes* the
   outer wing. Small (under 2%), but it is a real sign change and not an
   artifact — and it is the sort of result that only appears once the
   upstream field is modelled at all.

The magnitudes are inviscid, one-way-coupled induction. They set the
expectation; the deliverable is still Δα_sep from the separation march.

## First coupled result (2026-08-09) — and a problem

`aeolion_attachment_sweep` now takes optional thrust [N] and airspeed
[m/s] arguments and applies `DiskInductionField` as the solve's
`externalField`. At the transition point V = 12 m/s, T = 25 N
(v_i = 14.3 m/s), against the same sweep unpowered:

| α | ΔC_L | separated span off → on | separation x/c_n off → on |
|---|---|---|---|
| 6 | +0.026 | 64% → 59% | 0.853 → 0.866 |
| 8 | +0.028 | 100% → 100% | 0.810 → 0.828 |
| 12 | +0.032 | 100% → 100% | 0.653 → 0.698 |
| 16 | +0.036 | 100% → 100% | 0.468 → 0.498 |

The sign is right everywhere — lift up, separation aft — and at cruise
speed the same comparison gives roughly a third of this, which is the
expected scaling with v_i/V.

### The onset that "did not move" — resolved

A first pass reported the separation onset unchanged at α = 6° powered and
unpowered, which looked like Δα_sep ≈ 0. It was a measurement artifact,
and the way it failed is worth keeping.

"Onset" was defined as a **threshold crossing** — the first α at which any
station separates forward of x/c = 0.90 — evaluated on a **2° grid**. A
shift smaller than the grid spacing moves the crossing *within* an
interval without ever moving it *across* one, so the reported onset is
quantized to the grid and a sub-degree shift is invisible by construction.

The right measurement is the **horizontal shift of the separation curve**
x_sep(α), not the crossing of one threshold. Measured that way:

Measured on a **0.5° grid** (α = 3° to 16°, 27 stations), V = 12 m/s,
T = 25 N, v_i = 14.3 m/s:

| criterion x/c_n | α (fan off) | α (fan on) | Δα |
|---|---|---|---|
| 0.90 | 4.52° | 5.05° | **+0.52°** |
| 0.85 | 6.18° | 6.91° | **+0.73°** |
| 0.80 | 8.36° | 9.17° | **+0.81°** |
| 0.75 | 9.96° | 10.75° | **+0.79°** |
| 0.70 | 11.09° | 11.95° | **+0.86°** |
| 0.60 | 13.04° | 13.89° | **+0.85°** |
| 0.55 | 14.08° | 14.88° | **+0.80°** |

So the fan buys **Δα = +0.81°** (mean over x/c = 0.85 to 0.55), and the
consistency across criterion levels is what makes it a property of the
flow rather than of the definition. The shift is smaller near onset
(+0.52° at x/c = 0.90), where the separation point sits at the trailing
edge and the layer is least sensitive to the gradient.

The same shift computed by interpolating the original **2°** grid gave
+0.68 to +0.93° — agreeing with the fine grid to within 0.08° everywhere
separation is established, and differing most (0.30° against 0.52°) right
at onset, which is exactly where a coarse grid should be least trusted.
The coarse data always contained the answer; only the threshold *statistic*
destroyed it.

**Lesson for the paper:** report Δα as the shift of the curve at several
criterion levels, and state the criterion. A single threshold on a coarse
grid would have reported "no effect" from data that contains a 0.8° one.

## Scope

- upstream induction of the ducted fan, applied to the coupled
  wing–body–duct solve through the existing `externalField` hook
- separation march (the second paper's `AttachmentBoundaryLayer.h`) run
  with and without the fan field
- the deliverable: **Δα_sep as a function of fan operating point**,
  resolved spanwise — how many degrees of incidence the fan buys, and
  over which stations
- the transition sweep parameterized by velocity ratio
  `mu = Vinf / v_h`, with `v_h = sqrt(T / (2 rho A))` the hover induced
  velocity: `mu -> 0` is hover, large `mu` is cruise, and transition is
  the interesting middle

Deliberately **out of scope**: post-stall coefficients to 90 deg. That is
the paper this repo decided not to write, for the reasons recorded in
`TODO.md` — no closed-form verification exists past separation, and the
answer would be set by the four hand-tuned constants in
`AnalyticSectionModel`'s deep-stall blend rather than by the method.

## What has to be built

1. ~~**An upstream induction model, Level A**~~ — **DONE**.
   `Solver/DiskInduction.h`: the semi-infinite vortex cylinder, with a
   closed-form axis solution as its verification anchor, exposed as an
   `externalField` callable. `TestDiskInduction` pins it against momentum
   theory, refinement convergence, the annulus, and direction.
2. **Level B — the blade lattice's own field.** Azimuthal mean of the
   solved rotor's Biot–Savart field, to verify the uniform-loading
   cylinder against a real spanwise load distribution.
   `Solver::MeanInducedField` already does this construction in the
   vane→rotor direction; the rotor→airframe direction needs the same
   treatment, with the caveat its own comment raises — the rotor's wake is
   prescribed, so its near field is approximate in a way the vanes'
   explicit wake legs are not.
3. **Two-way coupling.** Currently one-way: the airframe sees the disk,
   the disk does not see the airframe, so every number here is a lower
   bound on the interaction. Airframe (static frame) and rotor (rotating,
   `fc.p = Omega`) cannot share one `FreestreamConditions`, so closing it
   needs the partitioned outer fixed point `SolveRotorVaneCoupled` already
   establishes.
4. **The transition operating line.** Thrust and airspeed are currently
   independent arguments; a real transition walks a coupled schedule
   (thrust roughly balancing weight minus wing lift as α and V change).
   The sweep should follow that line, not the full rectangle.

## Verification plan

The previous two papers rest on closed-form checks rather than mesh
refinement, and this one has to do the same or it is just a plausible
trend.

- **Vortex cylinder against momentum theory.** The upstream induction must
  give zero far upstream, `vi` at the disk, and `2 vi` far downstream, with
  the classical one-half relation at the disk plane. Exact, and it pins the
  whole model.
- **Zero-thrust limit.** With the fan off, the separation boundary must
  reproduce the second paper's α ≈ 6° exactly. This is the strongest
  regression check available and costs nothing to run.
- **Thrust consistency.** The induced field must be self-consistent with
  the disk loading that produced it (the momentum balance
  `ComputeSlipstreamBands` already solves in the downstream direction).
- **Monotonicity.** Δα_sep must increase with disk loading and decrease
  with airspeed; a model that failed this would be wrong in a way a
  single-point comparison could hide.

## Validation

RANS is planned for the second JoA paper's attachment/separation
predictions. This paper wants the same treatment and, if only one case can
be run, arguably needs it more: the induction magnitude is the whole
result, and unlike the attachment line it has no exact-solution anchor
on the real geometry.

## Stated limitations to carry into the draft

- Potential-flow induction. The fuselage boundary layer that an aft fan
  genuinely ingests (a real BLI effect) is not represented, and it acts in
  the direction of reducing the benefit.
- Quasi-steady. Transition is a manoeuvre, not a sequence of trim points.
- The duct is straight-walled in the schema — no lip camber (see
  `Geometry::DuctGeometry.h`), so lip suction is only as good as that
  idealization.
- The separation criterion inherits the bubble-bursting caveat from the
  second paper: turbulent separation is an upper bound on usable
  incidence, not a stall prediction.

## Style files (`style/`)

Copied from [../journal-of-aircraft/style/](../journal-of-aircraft/style/)
— see that folder's README for provenance, including the caveat that
`aiaa-tc.cls` / `aiaa.bst` were extracted by a hand-written docstrip
equivalent and should be regenerated with `latex aiaa.ins` before a real
submission. Use `\documentclass[submit]{aiaa-tc}` for journal formatting.

## Open items

- [ ] Confirm the fan operating points that bracket a real JETTAIL
      transition (rpm schedule vs. airspeed) — the sweep is only as
      meaningful as the operating line it walks.
- [ ] Decide whether the duct's own lip suction is reported separately
      from the fan's induction; they act on the same stations and a reader
      will ask which is which.
- [ ] Author block: same three authors as the other two drafts, same
      outstanding departments/emails/member grades.
- [ ] Check the AIAA duplicate-submission position against the second JoA
      article — different question and different result, but the same
      airframe and the same solver, so the relationship should be stated
      in the introduction rather than left for a reviewer to notice.
