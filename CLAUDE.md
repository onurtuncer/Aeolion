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
| Constants / enum values | PascalCase | `Pi`, `Failed` |
| Macros | `AE_`-prefixed SCREAMING_CASE | `AE_ASSERT`, `AE_CORE_ASSERT` |
| Namespace | PascalCase | `Aeolion` |

The whole codebase follows this convention. Headers live under
`include/Aeolion/` in a single `Aeolion` namespace with per-module nested
namespaces (`Aeolion::Vlm`, `Aeolion::Airfoil`, `Aeolion::Bemt`,
`Aeolion::BoundaryLayer`, `Aeolion::DragEstimate`, `Aeolion::MeshIO`,
`Aeolion::MeshSlice`, `Aeolion::VizExport`). Established aerodynamics/math
notation is preserved verbatim (`CL`, `CDi`, `Cm`, `alphaDeg`, `Vinf`,
`gamma`, `rho`, `Vec3.x/y/z`, and the `StabilityDerivatives` symbols);
only English-word identifiers are PascalCased.
