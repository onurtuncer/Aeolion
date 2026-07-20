# Aeolion

A small aerodynamics toolkit built around a 3D vortex lattice method:
lifting-surface VLM, component-buildup viscous drag estimation, and a
hover-safe propeller BEMT solver, with geometry loaded from a JSON
contract. Built incrementally and validated against closed-form theory at
each stage rather than assumed correct — see "Validation" below.

## Layout

```
include/Aeolion/     header-only library (this is the actual toolkit)
                     one folder per namespace; each holds its module header
                     plus that module's plain data structs, one per file
  Math/                                                   (Aeolion::Math)
    Vec3.h             3D double vector + dot/cross/axis-rotation helpers
    Constants.h        shared numeric constants and angle conversions
  VLM/                                                    (Aeolion::VLM)
    VLM.h              3D VLM core: horseshoe vortices, LAPACK LU solve,
                        sideslip, body rates, moments, stability
                        derivatives, external-velocity-field hook
    Panel.h              one horseshoe-vortex panel
    WingParams.h         parametric single-wing planform
    FreestreamConditions.h   flight condition (freestream, rates, ref point)
    ReferenceGeometry.h  coefficient normalization constants
    StationResult.h      per-spanwise-station output
    SolveResult.h        full solve result (coeffs, forces, per-surface)
    StabilityDerivatives.h   central-difference derivative table
  BEMT/                                                   (Aeolion::BEMT)
    BEMT.h             propeller BEMT (hover-safe: solves for induced
                        velocities directly, not induction factors) +
                        slipstream field for downstream control vanes
  Geometry/                                           (Aeolion::Geometry)
    HandoffContract.h  strict parser for the aeolion_geometry.json handoff
    PlanformStation.h    one spanwise planform station
    AirfoilSection.h     CST section shape at a station
    ControlSurface.h     hinged surface + which body it binds to
    MeshTopology.h       requested lattice discretization
    PropulsionSpec.h     propeller blade geometry for a BEMT run
  DragEstimate/                                   (Aeolion::DragEstimate)
    DragEstimate.h     viscous CD0: flat-plate skin friction + form factor
                        + interference factor component buildup

src/                  driver programs (link against the aeolion library)
  main.cpp               parametric single-wing demo (vlm_demo)
  GeometryContractCLI.cpp   solve a wing loaded from a JSON contract (aeolion_geometry)

tests/                 regression suite, wired into ctest
  TestVLMCore.cpp               VLM vs. thin-wing theory, Oswald efficiency
  TestBEMT.cpp                  BEMT vs. hard physical bounds (FOM <= 1,
                                 efficiency <= 1)
  TestPropVane.cpp              propwash -> vane control authority
                                 integration test
  TestSolverBackend.cpp         LAPACK dense-solve sanity check
  TestHandoffContract.cpp       JSON handoff parsing, contract invariants,
                                 surface binding, trapezoid reduction

cmake/                build modules (CompilerWarnings, LapackBackend) + vcpkg triplets
```

## Building

Requires a **C++23** compiler and a LAPACK provider (OpenBLAS, reference
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

Without CMake, the parametric demo builds directly (the JSON contract CLI
also needs nlohmann/json on the include path):

```
g++ -std=c++23 -O2 -Iinclude -o vlm_demo src/main.cpp -llapack -lblas
```

## Quick usage

```
# single parametric wing, prints CL/CDi/derivatives
./vlm_demo

# solve a wing loaded from a JSON geometry contract
./aeolion_geometry geometry.json
```

## Validation

Every module was checked against a closed-form or independently-known
result before being trusted, not just eyeballed for plausibility:

- **VLM core**: CL tracks thin-wing lifting-line theory within a few
  percent across aspect ratios 6-20; Oswald efficiency ~1.0 for a plain
  rectangular wing; panel-count convergence confirmed.
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

## VBAT geometry contract

`aeolion_geometry <path/to/aeolion_geometry.json>` reads the versioned,
SI-unit parametric handoff generated by `vbat-uav-notebooks`, constructs the
lattice directly, and runs a baseline solve. The same delivery contains STEP
files for CAD traceability and geometric cross-checks. STEP/STL conversion is
not used for optimization because tessellation and section extraction would
cut the derivative chain.

The JSON fields under `planform` are intentionally solver-native design
variables. A future CppAD integration should template geometry construction,
the influence assembly, linear solve, and objective reduction on a scalar
type, while parsing JSON into ordinary `double` starting values at the API
boundary. Discretization must stay fixed during one derivative evaluation.

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
