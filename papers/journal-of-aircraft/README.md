# Journal of Aircraft submission

**Venue:** AIAA Journal of Aircraft
**Type:** full research article
**Status:** drafting — `paper.tex` exists, Sections I–V written,
Section VI (Application) is a stub

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

- [ ] **Section VI, Application** — the substantive gap. Needs an α/β
      sweep on `tests/Data/AeolionGeometryHandoff-1.8.0.json`, which
      needs `panelbuilder` and therefore vcpkg. Planned content is a
      comment block in `paper.tex`; see `TODO.md` at the repo root.
- [ ] Figures: attachment line vs. span across α and β; `Rbar` against
      the 245/583 thresholds; traced fuselage surface streamlines
      coloured by Cp.
- [ ] Author block — co-authors (the SciTech draft carries three),
      departments, emails, AIAA member grades; acknowledgments/funding.
- [ ] Validation beyond closed-form verification. Sections V's checks are
      exact-solution comparisons (sphere, cylinder, matched asymptotics),
      which is strong for *verification* but is not *validation* — a
      wind-tunnel or higher-fidelity CFD attachment-line comparison would
      strengthen the article considerably.
- [ ] Regenerate `style/aiaa-tc.cls` and `style/aiaa.bst` with
      `latex aiaa.ins` before submission (see below).
