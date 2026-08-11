# Journal of Aircraft submission — Part I

**Venue:** AIAA Journal of Aircraft
**Type:** full research article, Part I of a two-part series
**Status:** drafting — split 2026-08-11: this part carries the method
(skin-flow topology, attachment lines) through its verification and the
attachment/sideslip application; the separation march, the separation
boundary, and the coefficient/derivative tables moved to Part II
(`papers/journal-of-aircraft-separation/`). Author block,
acknowledgments and external validation remain open in both parts.

**Corrected strip frames (2026-08-11).** The sweep originally measured
the section incidence in a camber-tilted panel frame, double-counting
the camber the section contour already carries and biasing α_n high by
~4.3° (the slope of the camber surface at the control point). The
driver now uses the true chord frame; offsets, α_n and every R̄ in the
prose changed accordingly, and CL/Cm were verified bit-identical.

**Figures/data pipeline.** Every number and figure in the application
section comes from the solver's own output — no screenshots, and no
hand-typed tables:

```
aeolion_attachment_sweep tests/Data/AeolionGeometryHandoff-1.8.0.json \
    papers/journal-of-aircraft/figures/attachment-sweep.json
cd papers/journal-of-aircraft/figures
python render-attachment-figures.py     # the four physics figures
cd .. && pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper
```

`aeolion_attachment_sweep` (app/AttachmentSweepExport.cpp) solves the
α × β matrix on the coupled system and exports coefficients, stability
derivatives, the per-strip attachment line, the separation survey and the
fuselage skin-flow topology. The intermediate `attachment-sweep.json` is
**not** tracked (the repo's `.gitignore` excludes `*.json`) — regenerate it
with the command above; the same convention the SciTech paper's
`lattice-solution.json` follows. The tables under `figures/tables/` are
generated and `\input`-ed by `paper.tex`, so a stale table is impossible:
rerun the scripts and the paper follows. The sweep takes ~20 minutes (55
attitudes × 13 solves for the central-difference derivatives).

## Scope

Stagnation points and attachment lines on the coupled wing–body solver,
for angle-of-attack and sideslip sweeps. The novel content is the method,
not an application of an existing VLM:

- the asymmetric treatment the two halves force: a source-panelled body's
  stagnation points are *found* by searching a real surface field, while
  a vortex lattice is a zero-thickness sheet on which the flow never
  stops and whose attachment line must come from a 2-D section solve
- skin-flow topology resolved in the panelling's own surface-index space,
  with critical points classified from a strain-rate tensor taken in a
  local orthonormal frame (not from contravariant components, whose chart
  factors do not converge away)
- the wing section posed in the plane **normal to the leading edge**,
  which is what makes sideslip's left/right asymmetry representable at
  all — and hence which wing is closer to leading-edge contamination
- everything an integral boundary-layer march needs to start: attachment
  location, `U_e(s)` from it, Thwaites' `θ₀`, the spreading metric `h(s)`,
  and Poll's attachment-line Reynolds number

The earlier planned scope — an applied wing-body-propeller case study of
the VBAT airframe (permeable base panels coupling propeller efflux,
BEMT-to-slipstream handoff) — is **not** this paper. The ducted-fan
vectored-thrust material is covered by the SciTech draft in
`papers/aiaa-scitech/`.

## Style files (`style/`)

Official AIAA LaTeX package, pulled unmodified from CTAN
(<https://ctan.org/tex-archive/macros/latex/contrib/aiaa>):

- `aiaa.dtx` / `aiaa.ins` — the authoritative docstrip source. Running
  `latex aiaa.ins` in that directory (any TeX distribution, e.g.
  MiKTeX) regenerates `aiaa-tc.cls` and `aiaa.bst` officially — **do
  this yourself before a real submission.**
- `aiaa-tc.cls`, `aiaa.bst` — already extracted here as a convenience,
  but via a hand-written docstrip-equivalent script (no local TeX
  engine was available to run the real `latex aiaa.ins` and diff
  against it). Treat these two files as unverified until regenerated
  or spot-checked against a compile.
- `template_basic.tex` / `template_advanced.tex` (+ rendered `.pdf`) —
  official starting-point documents. Use `\documentclass[submit]{aiaa-tc}`
  for journal formatting (12pt, double-spaced) as opposed to the
  default conference mode (10pt).
- `author_guide.tex/.pdf`, `aiaa.pdf` (class user manual),
  `bibtex_database.bib` — supporting reference material from the same
  package.

## Open items

- [x] **Section VI, Application** — written against the α × β matrix
      (α ∈ [−4, 16]°, β ∈ [−10, 10]°) on
      `tests/Data/AeolionGeometryHandoff-1.8.0.json`: attachment line,
      separation boundary, and a coefficient/derivative table restricted
      to the attached attitudes (see the pipeline note above).
- [x] Figures: attachment line vs. span across α; effective sweep and
      `Rbar` at β = ±10° against the 245/583 thresholds; fuselage
      attachment-node migration; traced fuselage surface streamlines
      coloured by Cp with the wing lattice; separation boundary.
- [ ] **`CDi` is not trustworthy here** and Sec. VI.E says so: the
      near-field force integration collects a spurious thrust from the
      discretized closed bodies, so CDi goes negative at zero lift and the
      apparent Oswald efficiency is 1.53 (impossible). A Trefftz-plane
      integration would fix it; see `TODO.md` at the repo root. CL and the
      moments are unaffected.
- [ ] Author block — co-authors (the SciTech draft carries three),
      departments, emails, AIAA member grades; acknowledgments/funding.
- [ ] Validation beyond closed-form verification. Sections V's checks are
      exact-solution comparisons (sphere, cylinder, matched asymptotics),
      which is strong for *verification* but is not *validation* — a
      wind-tunnel or higher-fidelity CFD attachment-line comparison would
      strengthen the article considerably.
- [ ] Regenerate `style/aiaa-tc.cls` and `style/aiaa.bst` with
      `latex aiaa.ins` before submission (see below).
