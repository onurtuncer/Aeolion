# Journal of Aircraft submission

**Venue:** AIAA Journal of Aircraft
**Type:** full research article
**Status:** planned

## Scope

Applied case study: coupled wing-body-propeller analysis of the VBAT
airframe using Aeolion. Core novel content is the numerical coupling, not
just an application of an existing VLM:

- permeable source panels coupling propeller efflux to the fuselage base
- wing circulation carried through the fuselage (source-panel + horseshoe
  vortex blocked linear system)
- BEMT-to-slipstream handoff feeding the body/control-vane interaction

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

## Open items before drafting

- [ ] Decide whether to frame this as analysis-only or as a stepping stone
      toward the CppAD-templated design-optimization path noted in the
      main README (affects scope and related-work framing)
- [ ] Validation case(s) beyond the unit-test physical bounds — need a
      wind-tunnel or higher-fidelity CFD comparison point for a full
      article (JOSS-level checks are not sufficient evidence here)
- [ ] Figures: lattice + source panel mesh, pressure/circulation
      distribution, propeller-airframe interaction results
