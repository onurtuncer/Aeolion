# AIAA SciTech Forum submission

**Venue:** AIAA SciTech Forum (conference proceedings)
**Type:** conference paper
**Status:** drafting — `paper.tex` + `references.bib` started; TODOs
marked inline (author block, paper number, intro literature review,
regenerated result numbers, acknowledgments)

**Figures:** `figures/lattice-views.pdf` is rendered from the solver's
own output, not screenshotted: `aeolion_lattice_export`
(app/LatticeFigureExport.cpp) solves the paper's reference case (the
TestVaneCascade coupled configuration; its thrust/pass count must match
Table 1's 6000 rpm row) and writes `figures/lattice-solution.json`,
then `figures/render-lattice-figures.py` (matplotlib) draws the vector
PDF. Rerun both after any solver change that moves the reference
numbers.

## Scope

The ducted-fan / vane slice, presented ahead of the full
[Journal of Propulsion and Power article](../journal-of-propulsion-and-power/README.md)
(AIAA explicitly permits extending one's own conference paper into a
journal article):

- duct paneled into the rotating-frame propeller solve
- downstream duct-jet vanes: panels, viscous solve, propulsive wrench
- vane wakes leaving along the local mean flow
- two-way rotor-vane coupling closed with a swirl momentum budget

Demonstration case: VBAT-class tail-sitter thrust vectoring. The point
of publishing this at a conference first is priority — this is the
freshest and most contested-space material — while journal-grade
validation is still being assembled.

Fallback venue if the SciTech cycle is missed: Aerospace Europe
Conference (Word template already kept in
[../ceas-aerospace-europe/style/](../ceas-aerospace-europe/style/)).

## Style files (`style/`)

Official AIAA LaTeX package from CTAN
(<https://ctan.org/tex-archive/macros/latex/contrib/aiaa>), copied from
[../journal-of-aircraft/style/](../journal-of-aircraft/style/) — see
that folder's README for provenance notes, including the caveat that
`aiaa-tc.cls` / `aiaa.bst` were extracted by a hand-written docstrip
equivalent and should be regenerated with `latex aiaa.ins` before a
real submission.

For a conference paper use `\documentclass{aiaa-tc}` in its default
conference mode (10pt, single-spaced) — not the `[submit]` journal
mode. Start from `template_basic.tex`. Check AIAA's current SciTech
author kit for the target year in case it diverges from the CTAN
package.

## Open items before drafting

- [ ] Check the SciTech abstract deadline for the target year (typically
      early June for the following January's forum)
- [ ] Pick the demonstration figure set: duct+blade+vane lattice, swirl
      budget closure, vectored-thrust wrench vs. vane deflection
- [ ] Decide the validation anchor for the conference version (physical
      bounds + momentum-budget closure may suffice at conference level;
      the journal version needs more — see JPP folder)
- [ ] Confirm scope split against the JPP article: this paper stakes the
      rotor-vane coupling claim; the journal version adds the viscous
      coupling levels and full validation
