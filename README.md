# Aeolion

A small, dependency-free (by default) aerodynamics toolchain built around a
vortex lattice method: 3D lifting-surface VLM, mesh-derived geometry from
STL, viscous drag estimation (both flat-plate buildup and a real 2D
panel+boundary-layer method), a hover-safe propeller BEMT solver, and a
raw-WebGL viewer. Built incrementally and validated against closed-form
theory at each stage rather than assumed correct — see "Validation" below.

## Layout

```
include/Aeolion/     header-only library (this is the actual toolkit)
  Vlm.h                3D VLM core: horseshoe vortices, LU-factored solve
                        (LAPACK/OpenBLAS), sideslip, body rates, moments,
                        stability derivatives, external-velocity-field hook
  Mesh.h               STL read/write (binary + ASCII, multi-solid)
  MeshSlice.h          slices a triangle-mesh lifting surface into VLM
                        panels (camber-line extraction from real geometry)
  DragEstimate.h       viscous CD0: flat-plate skin friction + form factor
                        + interference factor component buildup
  AirfoilPanel.h       2D panel method (Hess-Smith source+vortex)
  BoundaryLayer.h      Thwaites / Michel / Head / Squire-Young viscous BL
  Bemt.h               propeller BEMT (hover-safe: solves for induced
                        velocities directly, not induction factors) +
                        slipstream field for downstream control vanes
  VizExport.h          panel geometry + scalar fields -> JSON

src/                  driver programs (link against include/Aeolion)
  Main.cpp               parametric single-wing demo
  GeometryVlm.cpp        multi-surface (wing+htail+vtail+fuselage) CLI:
                          load STL per component, solve, derivatives,
                          drag buildup, JSON export for the viewer
  AirfoilPolar.cpp       standalone 2D airfoil polar tool (panel + BL)

tests/                 regression suite, wired into ctest
  fixtures/SyntheticWing.h      shared in-memory test-mesh generator
  TestVlmCore.cpp               VLM vs. thin-wing theory, Oswald efficiency
  TestPanelMethod.cpp           panel method vs. exact cylinder Cp, thin
                                 airfoil Cl_alpha
  TestBoundaryLayer.cpp         BL solver vs. Blasius + turbulent
                                 flat-plate correlations
  TestMeshSlice.cpp             mesh-derived panels vs. parametric panels
  TestBemt.cpp                  BEMT vs. hard physical bounds (FOM <= 1,
                                 efficiency <= 1)
  TestPropVane.cpp             propwash -> vane control authority
                                 integration test
  TestSolverBackend.cpp         LAPACK/OpenBLAS dense-solve sanity check

examples/
  vlm_viewer.html         raw-WebGL panel viewer (no Three.js, no VTK)
  spanwise_loading.png, polar.png   example plots
```

## Building

The dense linear solve calls LAPACK's dgetrf/dgetrs directly (Fortran ABI,
no LAPACKE header needed), backed by OpenBLAS or reference LAPACK/BLAS.
Install a provider first:

```
# Debian/Ubuntu:
apt-get install liblapacke-dev libopenblas-dev
```

Then:

```
mkdir build && cd build
cmake ..
make -j
ctest --output-on-failure
```

Without CMake, any single driver builds directly, e.g.:

```
g++ -std=c++17 -O2 -Iinclude -o geometry_vlm src/GeometryVlm.cpp \
    -llapack -lblas
```

## Quick usage

```
# single parametric wing, prints CL/CDi/derivatives
./vlm_demo

# a wing+tail aircraft from STL, with viscous drag and stability derivatives
./geometry_vlm --alpha 4 --vinf 25 --rho 1.225 --derivatives \
    --surface name=wing,file=wing.stl,span_axis=y,panels=24 \
    --surface name=htail,file=htail.stl,span_axis=y,panels=12 \
    --drag-component surface=wing,tc=0.12,sweep=8 \
    --json aircraft.json          # feed this into examples/vlm_viewer.html

# 2D airfoil polar via panel method + boundary layer (not the crude
# flat-plate estimate -- real transition prediction, separation warning)
./airfoil_polar --naca 2412 --re 500000 --cl 0.1,0.3,0.5,0.7,0.9
```

See the header comment at the top of `src/GeometryVlm.cpp` for the full
CLI reference (surfaces, sideslip/rates, drag components, CG/moments).

## Validation

Every module was checked against a closed-form or independently-known
result before being trusted, not just eyeballed for plausibility:

- **VLM core**: CL tracks thin-wing lifting-line theory within a few
  percent across aspect ratios 6-20; Oswald efficiency ~1.0 for a plain
  rectangular wing; panel-count convergence confirmed.
- **Panel method**: source-panel-only solve around a circle matches the
  exact potential-flow Cp(theta)=1-4sin^2(theta) to ~1e-7; symmetric NACA
  section gives exactly zero lift at alpha=0; Cl_alpha converges near 2*pi.
- **Boundary layer**: laminar Thwaites matches the exact Blasius solution
  to ~1%; turbulent Head's method converges to the known H~1.4 equilibrium
  shape factor and agrees with standard flat-plate correlations to ~15%
  (normal cross-method scatter for independent empirical closures).
- **Mesh-derived geometry**: STL-sliced panels match the parametric wing
  builder to <1% CL / <3% CDi across sweep+dihedral+twist combinations.
- **BEMT**: checked against hard physical bounds, not just plausibility --
  Figure of Merit and propulsive efficiency must both be <=1.0 (real
  thermodynamic constraints). This caught two real bugs during
  development: a torque-equation exponent error and a sign-flipped
  drag term in the thrust/torque resolution.
- **Solver**: LAPACK LU factorization (dgetrf/dgetrs), factorized once per
  geometry and reused across right-hand sides.

Bugs found and fixed along the way are left documented in the relevant
code comments rather than scrubbed from history -- several were only
caught by checking a *hard* physical bound (FOM<=1, efficiency<=1, exact
zero lift on a symmetric section) rather than a "looks about right" check,
which is worth keeping in mind if you extend this further.

## Known limitations (read before trusting results on a new geometry)

- **Thin, flat lattice**: no thickness effects beyond what `MeshSlice.h`
  captures via camber-line averaging. No boundary-layer/inviscid coupling
  in the 3D solver (that only exists in the standalone 2D `airfoil_polar`
  tool, and even there it's uncoupled/one-shot, same family as XFOIL minus
  its Newton iteration).
- **CDi is INDUCED drag only.** VLM is inviscid; total CD needs the
  `--drag-component` buildup (or, for higher fidelity per-section,
  `airfoil_polar`'s panel+BL cd) added on top.
- **Linear, attached-flow method.** No stall, no real separation physics
  (the BL solver flags a separation/laminar-bubble *warning*, doesn't
  model post-separation behavior). Don't trust it near stall.
- **Wake trails along the global x-axis**, not the true local freestream
  direction -- standard small-to-moderate-AoA VLM simplification.
- **BEMT is mid-fidelity blade element momentum theory**, not a
  substitute for measured prop data or a full rotor VLM (rotating lattice
  + helical wake) -- that would be a substantially larger, separate
  undertaking.
- **`MeshSlice.h` assumes one simply-connected thin-shell surface per
  STL file** -- export each lifting surface (wing, htail, vtail) as its
  own file; don't rely on automatic multi-body segmentation.
