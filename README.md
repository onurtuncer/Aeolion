# Aeolion

A small aerodynamics toolkit built around a 3D vortex lattice method:
lifting-surface VLM coupled to a fuselage source-panel mesh and
component-buildup viscous drag estimation, with geometry loaded from a
JSON contract. An interactive OpenGL viewer visualizes the lattice and
the live solve. Built incrementally and validated against closed-form
theory at each stage rather than assumed correct — see "Validation"
below. (The propeller BEMT solver that used to live here is a momentum
method, not a panel method, so it moved to its own project:
[onurtuncer/BEMT](https://github.com/onurtuncer/BEMT). This repo carries
no dependency on it.)

![aeolion_viewer rendering the wing + fuselage + duct lattice for a real geometry handoff, colored by circulation, with the solved coefficients shown live](doc/_static/viewer_airframe.png)

*`aeolion_viewer` rendering the coupled wing + fuselage + duct lattice built
by `PanelBuilder::LatticeBuilder` from a real geometry handoff (352 wing
panels colored by circulation `gamma`, 448 fuselage source panels, 160 duct
source panels -- all one coupled potential-flow solve). See
[doc/theory.rst](doc/theory.rst) for the math and
[doc/_static/viewer_airframe.png](doc/_static/viewer_airframe.png) full-size.*

## Layout

Every compiled/header-only module is a top-level component that owns its
own tests: `solver/tests/`, `panelbuilder/tests/`, and `tests/` for the
core header-only library itself. All of them are one `ctest` invocation
regardless of which folder a test lives in -- see "Building" below.

```
solver/                          header-only, part of `aeolion`
                                 kept as its own top-level component (like
                                 panelbuilder/ and viewer/) since
                                 every other module builds on top of it
  include/Aeolion/Solver/                              (Aeolion::Solver)
    Solver.h           3D VLM core: horseshoe vortices + fuselage source panels
                        in one blocked linear system, LAPACK LU solve, sideslip,
                        body rates, moments, stability derivatives,
                        external-velocity-field hook
    SourceInfluence.h    constant-strength quadrilateral source-panel kernel
                          (Hess & Smith in-plane, solid-angle normal component)
    WingParams.h         parametric single-wing planform
    FreestreamConditions.h   flight condition (freestream, rates, ref point)
    ReferenceGeometry.h  coefficient normalization constants
    StationResult.h      per-spanwise-station output
    SolveResult.h        full solve result (coeffs, forces, per-surface)
    StabilityDerivatives.h   central-difference derivative table
  tests/
    TestSolverCore.cpp            VLM vs. thin-wing theory, Oswald efficiency
    TestDenseSolve.cpp            LAPACK dense-solve sanity check
    TestSourcePanel.cpp           source-panel kernel vs. closed-form: far
                                   field, on-sheet jump, source-panelled sphere

panelbuilder/                    aeolion_panelbuilder (STATIC library)
                                 a component with compiled sources, not
                                 headers alone                    (Aeolion::PanelBuilder)
  include/Aeolion/PanelBuilder/
    PanelBuilder.h     LatticeBuilder: builds the solver's wing lattice, the
                        fuselage's source-panel mesh, AND the (disjoint)
                        duct's source-panel mesh from a parsed handoff --
                        CST camber surface, spanwise + chordwise
                        discretization, deflected control surfaces,
                        wing/body trim + carry-through, breakpoints on
                        every surface edge
  src/PanelBuilder.cpp
  tests/
    TestBodyPanels.cpp            fuselage panelling: closure, d'Alembert,
                                   Munk moment, frame convention
    TestDuctPanels.cpp            duct panelling: closure, enclosed annulus
                                   volume, d'Alembert, frame convention
    TestPanelBuilder.cpp          camber lifts at zero incidence, chordwise
                                   rows, control deflection, mesh lands on
                                   control surface edges
    TestAirframe.cpp              the whole chain coupled: real handoff ->
                                   wing + fuselage + duct lattice -> solve,
                                   including a control-surface deflection
                                   re-solved on the same cached builder
    TestPropellerLattice.cpp      rotating-frame blade VLM: thrust upstream,
                                   torque opposing rotation, blade-symmetry
                                   force cancellation, thrust falling with
                                   advance speed
    TestViscousCoupling.cpp       Level-2 sectional lift feedback: residual
                                   convergence, stall capping hover thrust,
                                   profile torque real and switchable,
                                   CST camber -> negative zero-lift angle

viewer/                          aeolion_viewer (exe, GL application)
                                 interactive OpenGL visualizer, not part of
                                 the library                    (not a namespace)
  include/Core/          Application (window + ImGui + solve loop), Window
                          (GLFW/glad bootstrap)
  include/Renderer/      Shader, Buffer, VertexArray, OrbitCamera (GL, no aero)
  include/Visualization/ LatticeRenderer (Solver::Panel + Lattice::SourcePanel
                          -> GPU meshes, colored by a chosen scalar field) +
                          PropellerRenderer (Geometry::Propeller -> twisted
                          blade surfaces, hub, disk; geometry only, no solve)
                          + ColorMap (viridis)

include/Aeolion/     header-only library (the rest of the toolkit)
                     one folder per namespace; each holds its module header
                     plus that module's plain data structs, one per file
  Math/                                                   (Aeolion::Math)
    Vec3.h             3D double vector + dot/cross/axis-rotation helpers
    Constants.h        shared numeric constants and angle conversions
  Lattice/                                             (Aeolion::Lattice)
    Panel.h            one horseshoe-vortex panel -- shared vocabulary
                        between the builder, the solver and the viewer,
                        so it belongs to none of them (re-exported into
                        Aeolion::Solver for source compatibility)
    SourcePanel.h      one constant-strength quadrilateral source panel --
                        the fuselage's atomic surface unit, the way Panel
                        is the wing's
  Geometry/                                           (Aeolion::Geometry)
    HandoffContract.h  strict parser for the aeolion_geometry.json handoff
    CstSurface.h       CST evaluation: camber mean line, its slope, and
                        camber-line arc length
    PlanformStation.h    one spanwise planform station
    AirfoilSection.h     CST section shape at a station
    ControlSurface.h     hinged surface + which body it binds to
    MeshTopology.h       requested lattice discretization
    PropulsionSpec.h     propeller blade geometry (the propulsion_bemt block)
    Propeller.h          the metric propeller consumers work with (blade
                          count, radii, chord/twist stations) -- built from
                          PropulsionSpec by HandoffContract.h's ToPropeller()
    BodyGeometry.h       fuselage body of revolution (axial x radius
                          stations)
    DuctGeometry.h        the duct: a single annular ring (chord, inner/outer
                           diameter, center), disjoint from the fuselage
    WingPlacement.h       where the wing sits relative to the body
  DragEstimate/                                   (Aeolion::DragEstimate)
    DragEstimate.h     viscous CD0: flat-plate skin friction + form factor
                        + interference factor component buildup

app/                  driver programs (link against the aeolion library)
  main.cpp               parametric single-wing demo (solver_demo)
  GeometryContractCLI.cpp   solve a wing loaded from a JSON contract (aeolion_geometry)

tests/                 the core header-only library's OWN tests, wired into
                       ctest the same way as every other module's tests/
  TestHandoffContract.cpp       JSON handoff parsing, contract invariants,
                                 surface binding, trapezoid reduction
  TestPropeller.cpp             contract-to-metric propeller conversion
                                 (r/R -> metres, refusals)
  Data/                  one real geometry handoff JSON per schema revision
                         tested against -- shared by every module's tests,
                         not just this folder's

cmake/                build modules (CompilerWarnings, LapackBackend,
                      AeolionTest -- the aeolion_add_test() every module's
                      tests/CMakeLists.txt calls) + vcpkg triplets

doc/                   Sphinx + Doxygen documentation (theory, API reference,
                       test-suite writeup) -- see "Documentation" below
```

## Building

Requires a **C++23** compiler and a LAPACK provider (OpenBLAS, reference
LAPACK, MKL, Accelerate, ...). The dense solve calls LAPACK's dgetrf/dgetrs
through their Fortran ABI directly — no LAPACKE C header needed. The
viewer additionally needs GLFW, glad, glm, and Dear ImGui (OpenGL 3.3
core); set `AEOLION_BUILD_VIEWER=OFF` to skip it and build only the
header-only library, `panelbuilder`, and the CLIs.

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
to have every dependency (LAPACK, nlohmann-json, and the viewer's GL stack)
fetched and built for you:

```
cmake --preset windows   # or: linux
cmake --build --preset windows
ctest --test-dir build/windows --output-on-failure
```

Without CMake, the parametric demo builds directly (the JSON contract CLI
also needs nlohmann/json on the include path):

```
g++ -std=c++23 -O2 -Iinclude -Isolver/include -o solver_demo app/main.cpp -llapack -lblas
```

## Quick usage

```
# single parametric wing, prints CL/CDi/derivatives
./solver_demo

# solve a wing loaded from a JSON geometry contract
./aeolion_geometry geometry.json
```

## Viewer

`aeolion_viewer` is an interactive OpenGL 3.3 visualizer: orbit camera,
panel-count/planform/freestream sliders, live CL/CDi/Cm/stability-derivative
readouts, and a scalar-field colormap (circulation, sectional cl, lift per
span) over the lattice. With no arguments it opens on the built-in
parametric-wing demo.

```
./aeolion_viewer                              # interactive, parametric wing demo
./aeolion_viewer --geometry geometry.json     # load a real handoff instead (wing + fuselage)
./aeolion_viewer --frames N                   # render N frames then exit (smoke test)
./aeolion_viewer --screen propeller           # open on the propeller screen
./aeolion_viewer --geometry geometry.json --screenshot out.ppm --frames 5
                                               # headless-ish capture: PPM screenshot, then exit
```

In handoff mode the Planform sliders are replaced by a read-only geometry
summary (design ID, panel counts, trim eta); the Freestream and Display
controls still apply. Screenshot output is a binary PPM (P6); convert with
any image tool (e.g. `Pillow`: `Image.open("out.ppm").save("out.png")`).

A menu-bar toggle (or `--screen propeller`) switches to the **propeller
screen**: the handoff's `propulsion_bemt` blade geometry (or a built-in
default prop) meshed into a rotating-frame vortex lattice
(`PanelBuilder::BuildPropellerLattice`) and solved by the same VLM core as
the airframe -- the rotation enters as the solver's roll rate about +x, so
every blade panel sees its true `Omega x r` onset flow. The blades sit on
their CST camber surface when the contract states blade sections (schema
>= 1.8.0), with chords wrapped on their radius cylinders, one Weissinger
row per radial strip, and trailing legs along the local relative wind
plus the momentum-theory hover inflow (a prescribed linearized helical
wake). The solved lattice renders colored by circulation / sectional cl /
lift per span, with RPM, axial-speed, and density controls and
dimensional thrust / torque / power / disk-loading readouts. On top of
the lattice sits **Level-2 viscous coupling** (`Solver::SolveViscousCoupled`,
on by default): per radial strip the lattice supplies `alpha_eff`, a 2-D
viscous section model supplies `cl/cd(alpha_eff, Re, Ma)` (today an
analytic polar -- finite lift slope about the CST camber line's
thin-airfoil zero-lift angle, smooth stall saturation, Re-scaled drag
polar; the interface is the seam a boundary-layer section solver will
implement), and the circulation relaxes until the sectional-lift
residual vanishes. That caps the stalled root loading and adds real
**profile torque**, reported separately from the induced part.
Remaining caveats: quasi-steady prescribed straight-line wake (no
roll-up or contraction), single chordwise row (integrated loads, not
pressure distributions), no thickness, incompressible (Ma carried but
unused).

## Validation

Every module was checked against a closed-form or independently-known
result before being trusted, not just eyeballed for plausibility:

- **Solver core**: CL tracks thin-wing lifting-line theory within a few
  percent across aspect ratios 6-20; Oswald efficiency ~1.0 for a plain
  rectangular wing; panel-count convergence confirmed.
- **Dense solve**: LAPACK LU factorization (dgetrf/dgetrs), factorized once
  per geometry and reused across right-hand sides.
- **Source-panel kernel**: matches the closed-form far-field point source,
  the analytic on-sheet velocity jump, and a source-panelled sphere's exact
  potential-flow surface speed.
- **Fuselage panelling**: a closed body of revolution encloses its stated
  volume with outward normals; feels zero net force at any incidence
  (d'Alembert) and a destabilizing pitching moment at incidence (the Munk
  moment) -- the dominant effect a fuselage has on airframe stability.
- **Wing-body coupling**: trimming the wing lattice where it enters the
  fuselage costs lift; the carry-through vortex restores essentially all
  of it (round-trips to the untrimmed full-span wing within 2%); the body
  then adds a modest lift increment on top (blockage + upwash).

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

## Documentation

`doc/` is a Sphinx + Doxygen/Breathe site: `doc/theory.rst` derives the
full linear system (influence-matrix entries, boundary condition, near-field
Kutta-Joukowski and pressure force integration, coefficient normalization)
and the drag-buildup math; `doc/tests.rst` walks what each module's
tests actually validate and why; `doc/api.rst` pulls the Doxygen
comments in every header into one API reference. Build it locally with
Doxygen + Graphviz and the Python packages in `doc/requirements.txt`:

```
pip install -r doc/requirements.txt
doxygen build/Doxyfile   # generated by doc/CMakeLists.txt, or substitute
                         # @CMAKE_SOURCE_DIR@/@DOXYGEN_OUTPUT_DIR@ by hand
sphinx-build -b html -Dbreathe_projects.Aeolion=<doxygen-xml-dir> doc build/sphinx
```

`.github/workflows/deploy-docs.yml` builds and publishes the same site to
GitHub Pages on every push to `main`.

## Known limitations (read before trusting results on a new geometry)

- **Wing lattice is thin**: panels sit on the CST camber mean line
  (`Geometry/CstSurface.h`), not on the true upper/lower surface -- no
  thickness effects beyond what camber-line placement captures. The
  fuselage, by contrast, is a real 3D source-panel mesh with thickness.
- **CDi is INDUCED drag only.** VLM is inviscid; total CD needs
  `Aeolion::DragEstimate::EstimateCD0()`'s component buildup added on top
  -- there is no CLI flag for this yet, it is a library call you make
  yourself (see `doc/theory.rst`).
- **Linear, attached-flow method.** No stall, no separation physics.
  Don't trust it near stall.
- **Wake trails along the global x-axis**, not the true local freestream
  direction -- standard small-to-moderate-AoA VLM simplification.
- **Wing/body placement needs schema >= 1.5.0's `planform.placement`.**
  Without it (older contracts, or a newer one that omits the optional
  block) the wing and the body both default to the origin, which puts the
  wing root at the fuselage nose rather than wherever it actually sits.
  The screenshot above has placement stated, so its wing sits at the
  contract's real anchor point instead.
