# Aeolion

A small aerodynamics toolkit built around a 3D vortex lattice method:
lifting-surface VLM, mesh-derived geometry from STL, component-buildup
viscous drag estimation, and a hover-safe propeller BEMT solver. Built
incrementally and validated against closed-form theory at each stage
rather than assumed correct — see "Validation" below.

## Layout

```
include/Aeolion/     header-only library (this is the actual toolkit)
  Math/
    Vec3.h             3D double vector + dot/cross/axis-rotation helpers
  Types/               the plain data model (no algorithms, no LAPACK)
    Panel.h              one horseshoe-vortex panel
    WingParams.h         parametric single-wing planform
    FreestreamConditions.h   flight condition (freestream, rates, ref point)
    ReferenceGeometry.h  coefficient normalization constants
    StationResult.h      per-spanwise-station output
    SolveResult.h        full solve result (coeffs, forces, per-surface)
    StabilityDerivatives.h   central-difference derivative table
  VLM.h                3D VLM core: horseshoe vortices, LAPACK LU solve,
                        sideslip, body rates, moments, stability
                        derivatives, external-velocity-field hook
  Mesh.h               STL read/write (binary + ASCII, multi-solid)
  MeshSlice.h          slices a triangle-mesh lifting surface into VLM
                        panels (camber-line extraction from real geometry)
  DragEstimate.h       viscous CD0: flat-plate skin friction + form factor
                        + interference factor component buildup
  Bemt.h               propeller BEMT (hover-safe: solves for induced
                        velocities directly, not induction factors) +
                        slipstream field for downstream control vanes

src/                  driver programs (link against the aeolion library)
  main.cpp               parametric single-wing demo (vlm_demo)
  GeometryVLM.cpp        multi-surface (wing+htail+vtail+fuselage) CLI:
                          load STL per component, solve, derivatives,
                          drag buildup (geometry_vlm)

tests/                 regression suite, wired into ctest
  fixtures/SyntheticWing.h      shared in-memory test-mesh generator
  TestVLMCore.cpp               VLM vs. thin-wing theory, Oswald efficiency
  TestMeshSlice.cpp             mesh-derived panels vs. parametric panels
  TestBemt.cpp                  BEMT vs. hard physical bounds (FOM <= 1,
                                 efficiency <= 1)
  TestPropVane.cpp              propwash -> vane control authority
                                 integration test
  TestSolverBackend.cpp         LAPACK dense-solve sanity check

cmake/                build modules (CompilerWarnings, LapackBackend) + vcpkg triplets
```

## Building

Requires a **C++20** compiler and a LAPACK provider (OpenBLAS, reference
LAPACK, MKL, Accelerate, ...). The dense solve calls LAPACK's dgetrf/dgetrs
through their Fortran ABI directly — no LAPACKE C header needed.

```
# Debian/Ubuntu system LAPACK:
apt-get install liblapack-dev libblas-dev
```

Then:

```
mkdir build && cd build
cmake ..
make -j
ctest --output-on-failure
```

CMake presets (`CMakePresets.json`) wire up a vcpkg toolchain if you prefer
to have the dependency fetched for you.

Without CMake, a single driver builds directly:

```
g++ -std=c++20 -O2 -Iinclude -o geometry_vlm src/GeometryVLM.cpp -llapack -lblas
```

## Quick usage

```
# single parametric wing, prints CL/CDi/derivatives
./vlm_demo

# a wing+tail aircraft from STL, with viscous drag and stability derivatives
./geometry_vlm --alpha 4 --vinf 25 --rho 1.225 --derivatives \
    --surface name=wing,file=wing.stl,span_axis=y,panels=24 \
    --surface name=htail,file=htail.stl,span_axis=y,panels=12 \
    --drag-component surface=wing,tc=0.12,sweep=8
```

See the header comment at the top of `src/GeometryVLM.cpp` for the full
CLI reference (surfaces, sideslip/rates, drag components, CG/moments).

## Validation

Every module was checked against a closed-form or independently-known
result before being trusted, not just eyeballed for plausibility:

- **VLM core**: CL tracks thin-wing lifting-line theory within a few
  percent across aspect ratios 6-20; Oswald efficiency ~1.0 for a plain
  rectangular wing; panel-count convergence confirmed.
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
caught by checking a *hard* physical bound (FOM<=1, efficiency<=1) rather
than a "looks about right" check, which is worth keeping in mind if you
extend this further.

## Known limitations (read before trusting results on a new geometry)

- **Thin, flat lattice**: no thickness effects beyond what `MeshSlice.h`
  captures via camber-line averaging.
- **CDi is INDUCED drag only.** VLM is inviscid; total CD needs the
  `--drag-component` buildup added on top.
- **Linear, attached-flow method.** No stall, no separation physics.
  Don't trust it near stall.
- **Wake trails along the global x-axis**, not the true local freestream
  direction -- standard small-to-moderate-AoA VLM simplification.
- **BEMT is mid-fidelity blade element momentum theory**, not a
  substitute for measured prop data or a full rotor VLM (rotating lattice
  + helical wake) -- that would be a substantially larger, separate
  undertaking.
- **`MeshSlice.h` assumes one simply-connected thin-shell surface per
  STL file** -- export each lifting surface (wing, htail, vtail) as its
  own file; don't rely on automatic multi-body segmentation.
