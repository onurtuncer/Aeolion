# Journal of Propulsion and Power submission

**Venue:** AIAA Journal of Propulsion and Power
**Type:** full research article (extension of the
[SciTech conference paper](../aiaa-scitech/README.md))
**Status:** planned

## Scope

The ducted propulsor as a system — JPP rather than Journal of Aircraft
because the deliverable is a propulsive wrench and a swirl budget,
which is propulsion-integration territory:

- rotating-frame propeller VLM: blades on CST camber surfaces wrapped
  on radius cylinders, trailing legs along the local relative wind
  (the choice that keeps hover from diverging), thrust `= -Di`,
  shaft torque `= -Mx`
- leveled viscous coupling: Level 2 sectional-lift feedback, Level 3
  transpiration-coupled boundary-layer section solver
- duct paneled into the propeller solve (ducted-fan interaction)
- downstream duct-jet vanes and the two-way rotor-vane coupling with
  the swirl momentum budget

Scope split against the sibling papers: this paper never claims the
fuselage circulation carry-through (that is the
[Journal of Aircraft paper](../journal-of-aircraft/README.md)'s), and
it extends — not duplicates — the SciTech paper by adding the viscous
coupling levels and full validation.

## Style files (`style/`)

Official AIAA LaTeX package from CTAN
(<https://ctan.org/tex-archive/macros/latex/contrib/aiaa>), copied from
[../journal-of-aircraft/style/](../journal-of-aircraft/style/) — see
that folder's README for provenance notes, including the caveat that
`aiaa-tc.cls` / `aiaa.bst` were extracted by a hand-written docstrip
equivalent and should be regenerated with `latex aiaa.ins` before a
real submission.

Use `\documentclass[submit]{aiaa-tc}` for journal formatting (12pt,
double-spaced). Start from `template_basic.tex` or
`template_advanced.tex`.

## Open items before drafting

- [ ] Validation anchors — the schedule gate for the whole paper:
  - [ ] isolated propeller: published thrust/torque data (e.g. APC
        performance datasets or a UIUC prop-database case) against the
        rotating-frame VLM at Level 1/2/3
  - [ ] ducted fan with vanes: published ducted-fan/swirl-recovery-vane
        experimental data, VBAT program tunnel data, or a RANS
        comparison run for this paper
- [ ] Decide how the swirl momentum budget is presented: as a closure
      check on the solve, or as a predictive quantity validated in its
      own right
- [ ] Figures: duct+blade+vane lattice, viscous-level convergence,
      thrust/torque vs. advance ratio against the anchor data,
      vectored wrench vs. vane deflection
- [ ] Cite the JOSS paper's DOI for the software once minted
