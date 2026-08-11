# Journal of Aircraft submission — Part II (post-stall)

**Venue:** AIAA Journal of Aircraft
**Type:** full research article, Part II of a two-part series
**Status:** drafting — created 2026-08-11. Part I
(`papers/journal-of-aircraft/`) carries the aerodynamics *up to* the
separation boundary (attachment lines, the separation march, the
attached-flow coefficient tables); this part reports what lies beyond:
the post-stall map to α = 90° and β = 30°, from the sectional-feedback
coupling with the anchored post-stall section model
(`Solver/PostStallSection.h` — Kirchhoff attenuation from the computed
separation points, Viterna-Corrigan deep stall, Rayleigh centre of
pressure). Author block and member grades must stay synchronized with
Part I.

**Headline numbers (2026-08-11 data):** σ-collapse across β within 13%;
plate regime by α ≈ 55–65°; CN(90°) = 1.52 vs the Viterna ceiling 1.218
(+25%; the hand-tuned baseline sat +73%); x_cp walks 0.25c → 0.52c,
landing on Rayleigh's mid-chord; emergent CLmax = 1.76 at α = 18°
(upper bound — bubble bursting is not modelled, see Part I).

**Figures/data pipeline.** Two sweep runs feed the figures — the
anchored model and the Phase-0 baseline it replaces:

```
aeolion_poststall_sweep tests/Data/AeolionGeometryHandoff-1.8.0.json \
    papers/journal-of-aircraft-poststall/figures/poststall-map-anchored.json \
    25 0.05 0 1000 all anchored
aeolion_poststall_sweep tests/Data/AeolionGeometryHandoff-1.8.0.json \
    papers/journal-of-aircraft-poststall/figures/poststall-map-baseline.json \
    25 0.05 0 1000 all analytic
cd papers/journal-of-aircraft-poststall/figures
python render-poststall-figures.py
cd .. && pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper
```

The JSONs are not tracked (`.gitignore` excludes `*.json`); regenerate
with the commands above (~15 min each). The coupling knobs in the
command line matter: plain damped iteration (relaxation 0.05, Anderson
off) is what suppresses the checkerboard multistability the paper's
method section describes. Style files under `style/` are the same
hand-extracted `aiaa-tc.cls`/`aiaa.bst` as Part I's — regenerate with
`latex aiaa.ins` before a real submission.

Full study record: repo `TODO.md` §3c (Phase 0/1 of the post-separation
study).
