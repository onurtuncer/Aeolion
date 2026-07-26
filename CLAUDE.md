# Aeolion

## Naming conventions

Follow the naming conventions used by TheCherno's **Hazel** game engine.

| Element | Style | Example |
|---|---|---|
| Files | PascalCase, `.h` / `.cpp` | `VortexLatticeMethod.h`, `AirfoilPanel.cpp` |
| Classes / structs / enums | PascalCase | `Vec3`, `VertexBuffer` |
| Methods / free functions | PascalCase | `Normalized()`, `OnUpdate()` |
| Member fields | `m_` + PascalCase | `m_ProjectionMatrix` |
| Static members | `s_` + PascalCase | `s_Instance` |
| Globals | `g_` + PascalCase | `g_Context` |
| Locals / parameters | camelCase | `controlPoint`, `vertexCount` |
| Constants / enum values | PascalCase | `Failed` |
| Macros | `AE_`-prefixed SCREAMING_CASE | `AE_ASSERT`, `AE_CORE_ASSERT` |
| Namespace | PascalCase, `::` form | `namespace Aeolion::Solver {` |

The whole codebase (C++23) follows this convention. Headers live under
`include/Aeolion/`, **one folder per namespace** — `Math/`, `Lattice/`,
`BEMT/`, `Geometry/`, `DragEstimate/`. Nothing sits naked at the top of
`include/Aeolion/`. Each folder holds its module header (`Geometry/HandoffContract.h`, …)
alongside that module's plain data structs, one struct per file
(`Geometry/ControlSurface.h`, …). A new module means a new folder, not a
new top-level header. Namespaces are written with the `::` form —
`Aeolion::Solver`, `Aeolion::BEMT`, `Aeolion::DragEstimate`,
`Aeolion::Geometry`, `Aeolion::Lattice` — never the nested-brace form.

**Top-level components.** A module that is large enough to own its build
target, or that has compiled sources rather than headers alone, becomes a
sibling of `include/` with its own `include/`(`/src/`) split, the way
`viewer/` does. Each such module also owns its own `tests/` — there is no
single central test directory; `tests/` at the repo root is the core
header-only library's *own* test folder, not a shared one:

| Folder | Target | Why it is top-level |
|---|---|---|
| `solver/` | header-only, part of `aeolion` | the VLM core every other module builds on |
| `panelbuilder/` | `aeolion_panelbuilder` (STATIC) | a component with compiled sources |
| `bemt/` | `aeolion_bemt` (STATIC) | a component with compiled sources |
| `viewer/` | `aeolion_viewer` (exe) | GL application, not part of the library |

`Lattice/Panel.h` stays in `include/` despite being shared by all of them:
it is pure data vocabulary, so it belongs to no one module. It is
re-exported into `Aeolion::Solver` the same way `Math/Vec3.h` is, so
`Solver::Panel` keeps resolving. `BEMT/PropGeometry.h` follows the same
rule relative to `bemt/`: it is the blade data vocabulary a geometry
handoff builds (`Geometry::ToPropGeometry`), and building that data never
requires calling into the BEMT solver itself, so it stays header-only in
`include/Aeolion/BEMT/` rather than moving into `bemt/include/`.

Prefer STATIC over SHARED for Aeolion's own libraries — a DLL would need
`__declspec` plumbing on every exported symbol for no ABI-stability or
plugin benefit.

Established aerodynamics/math notation is preserved verbatim; the solver is
3D-only, and geometry is loaded from a JSON contract (no STL/mesh).

Established aerodynamics/math notation is preserved verbatim (`CL`, `CDi`,
`Cm`, `alphaDeg`, `Vinf`, `gamma`, `rho`, `Vec3.x/y/z`, and the
`StabilityDerivatives` symbols); only English-word identifiers are
PascalCased. Use `std::numbers::pi` for π — do not define a local `Pi`.
