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

**Scope (agreed 2026-08-11, narrowed 2026-08-14):** Part II carries
the quasi-steady post-separation program — Phase 1 (this draft's
anchored map, plus the still-pending Ostowari–Naik validation and the
up/down hysteresis map; the Sheldahl–Klimas Re=3.6e5 comparison is in)
and Phase 2 (Maskew–Dvorak double-wake polars on the Hess–Smith
section solve). The unsteady cross-check (Phase 3: 2D LESP
discrete-vortex sections, 3D particle wake) was removed from the paper
on 2026-08-14 — that work lives in the repo (TODO.md §3d/3e), not in
this part. The declared-but-pending pieces are outlined in the paper's
"Validation and the remaining tiers" section with TODO markers; no
number in the results sections depends on them.

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

The section-level validation figure has its own one-liner:

```
aeolion_section_validation \
    papers/journal-of-aircraft-poststall/figures/section-validation.json
cd papers/journal-of-aircraft-poststall/figures
python render-validation-figure.py   # prints the prose numbers
```

`NACA_0015_SheldahlKlimas.dat` is the Sandia-distributed 360-degree
table (SAND80-2114 data, taken verbatim from
github.com/sandialabs/CACTUS `test/Airfoil_Section_Data/NACA_0015.dat`)
and IS tracked — it is the validation's frozen reference. The model
polar `section-validation.json` is tracked too, as the data of record
behind Fig. "section-validation" (it is small and the paper quotes its
numbers).

The JSONs are not tracked (`.gitignore` excludes `*.json`); regenerate
with the commands above (~15 min each). The coupling knobs in the
command line matter: plain damped iteration (relaxation 0.05, Anderson
off) is what suppresses the checkerboard multistability the paper's
method section describes. Style files under `style/` are the same
hand-extracted `aiaa-tc.cls`/`aiaa.bst` as Part I's — regenerate with
`latex aiaa.ins` before a real submission.

Full study record: repo `TODO.md` §3c (Phase 0/1 of the post-separation
study).
