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
| Namespace | PascalCase, `::` form | `namespace Aeolion::VLM {` |

The whole codebase (C++20) follows this convention. Headers live under
`include/Aeolion/`: the math primitive in `Math/` (`Aeolion::Math`), the
plain data structs one-per-file under `Types/`, and the modules
(`VLM.h`, `DragEstimate.h`, `BEMT.h`, `GeometryContract.h`). Namespaces
are written with the `::` form — `Aeolion::VLM`, `Aeolion::BEMT`,
`Aeolion::DragEstimate`, `Aeolion::Geometry` — never the nested-brace form.
The solver is 3D-only; geometry is loaded from a JSON contract (no STL/mesh).

Established aerodynamics/math notation is preserved verbatim (`CL`, `CDi`,
`Cm`, `alphaDeg`, `Vinf`, `gamma`, `rho`, `Vec3.x/y/z`, and the
`StabilityDerivatives` symbols); only English-word identifiers are
PascalCased. Use `std::numbers::pi` for π — do not define a local `Pi`.
