# JOSS submission

**Venue:** Journal of Open Source Software
**Type:** software description paper (~250-1000 words + metadata), not a novel-methods paper
**Status:** planned

## Scope

Describe Aeolion itself: 3D VLM + BEMT + component-buildup drag toolkit,
JSON-contract geometry input, validation methodology (thin-wing theory,
Oswald efficiency, FOM<=1/eta<=1 physical bounds). Draws directly on the
README's "Validation" and "Known limitations" sections.

## Style files

JOSS has no LaTeX class — submissions are Markdown, built to PDF by
JOSS's own pandoc-based pipeline. `paper.md` (YAML front matter +
required section skeleton: Summary, Statement of need, State of the
field, Software design, Research impact statement, AI usage
disclosure, Acknowledgements, References) and `paper.bib` are here,
transcribed from JOSS's current documented format
(<https://joss.readthedocs.io/en/latest/paper.html>). Sections still
need actual content — placeholders are HTML comments.

## Checklist (JOSS review criteria)

- [ ] Open license (MIT — already in place)
- [ ] Installation instructions (README — already in place)
- [ ] Example usage (`solver_demo`, `aeolion_geometry` — already in place)
- [ ] Automated tests (ctest suite — already in place)
- [ ] Community guidelines (CONTRIBUTING / issue template — not yet)
- [ ] `paper.md` + `paper.bib` (JOSS's required format)
