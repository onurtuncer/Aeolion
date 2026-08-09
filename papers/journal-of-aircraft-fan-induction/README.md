# Journal of Aircraft submission — aft-fan inflow induction

**Venue:** AIAA Journal of Aircraft
**Type:** full research article
**Status:** scoped — no draft yet. The method it depends on does not exist
in the solver yet; see "What has to be built" below.

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

**But the separation ONSET does not move.** It is α = 6° powered and
unpowered. The fan pushes the separation point aft *once the wing has
separated*; it does not delay the incidence at which separation begins.
That is a real result and it is not the one the paper was scoped around,
so it needs confronting rather than presenting as Δα_sep ≈ 0:

- The onset is set by the station that separates *first*, and the
  criterion is a threshold crossing (x/c < 0.90). A 1–4% chord shift
  moves the crossing within an α interval but not across one — the α grid
  is 2° and may simply be too coarse to resolve a sub-degree shift.
  **Refining α near onset is the first thing to try.**
- The onset station is inboard at |2y/b| ≈ 0.15, which is where the
  induction is strongest, so the mechanism is acting in the right place;
  the question is only whether it is worth a degree.
- Δα_sep may genuinely be small at these thrust settings, in which case
  the honest paper reports the *separation-extent* shift rather than an
  onset shift, and the headline becomes "the fan does not delay stall
  onset but substantially reduces separated area once stalled" — still a
  result, and a more surprising one.

Either way the deliverable needs restating from "Δα_sep" to something the
data can actually support. Do not write the abstract before the refined-α
sweep settles this.

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

Nothing here is a small edit; this is the honest cost.

1. **An upstream induction model.** Two levels, the cheap one verifying
   the expensive one:
   - *Level A* — semi-infinite vortex cylinder / actuator disk, which has
     a **closed-form** upstream axial induction. This is what makes the
     paper verifiable in the way the previous two are.
   - *Level B* — azimuthal mean of the solved blade lattice's own
     Biot–Savart field. `Solver::MeanInducedField` already does exactly
     this construction in the vane→rotor direction; the rotor→airframe
     direction needs the same treatment, with the caveat its own comment
     raises — the rotor's wake is prescribed, so its near field is
     approximate in a way the vanes' explicit wake legs are not.
2. **A coupled transition solve.** Airframe (static frame) and rotor
   (rotating frame, `fc.p = Omega`) cannot share one
   `FreestreamConditions`, so this needs the partitioned outer fixed point
   `SolveRotorVaneCoupled` already established for the rotor–vane problem.
3. **Tests**, in the repo's style — see the verification plan below.

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
