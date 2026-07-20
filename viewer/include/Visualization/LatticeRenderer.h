// Visualization/LatticeRenderer.h
//
// Turns the solver's data model (Solver::Panel + Solver::SolveResult) into GPU
// meshes and draws them: panel quads colored by a chosen scalar field
// through the viridis colormap, plus a line overlay (panel outlines, bound
// vortex segments, control-point normals, reference grid, body axes).
//
// This is the boundary between aerodynamics and rendering -- nothing above
// it touches GL, nothing below it knows what a panel is.

#pragma once

#include "Renderer/Buffer.h"
#include "Renderer/Shader.h"
#include "Renderer/VertexArray.h"

#include "Aeolion/Solver/Panel.h"
#include "Aeolion/Solver/SolveResult.h"
#include "Aeolion/Solver/StationResult.h"

#include <glm/glm.hpp>

#include <cstddef>
#include <vector>

namespace Aeolion::Viewer {

enum class FieldMode {
    Gamma = 0,      // circulation per panel
    LocalCl,        // sectional lift coefficient
    LiftPerSpan,    // dL/dy [N/m]
};

struct LatticeDisplayOptions {
    FieldMode Field = FieldMode::Gamma;
    bool ShowWireframe = true;
    bool ShowNormals = false;
    bool ShowGrid = true;
    bool ShowAxes = true;
};

class LatticeRenderer {
public:
    LatticeRenderer();

    // Rebuild all vertex data from a solved lattice. `span` scales the
    // grid, axes, and normal-vector glyphs.
    void Update(const std::vector<Solver::Panel>& panels, const Solver::SolveResult& result,
                const LatticeDisplayOptions& options, double span);

    void Draw(const glm::mat4& viewProjection) const;

    // Field range mapped to the colormap ends by the last Update().
    [[nodiscard]] double FieldMin() const { return m_FieldMin; }
    [[nodiscard]] double FieldMax() const { return m_FieldMax; }

private:
    void AppendLine(std::vector<float>& lines, const glm::vec3& a, const glm::vec3& b,
                    const glm::vec3& color) const;

    Shader m_PanelShader;
    Shader m_LineShader;

    VertexArray m_PanelArray;
    VertexBuffer m_PanelBuffer;
    std::size_t m_PanelVertexCount = 0;

    VertexArray m_LineArray;
    VertexBuffer m_LineBuffer;
    std::size_t m_LineVertexCount = 0;

    double m_FieldMin = 0.0;
    double m_FieldMax = 0.0;
};

} // namespace Aeolion::Viewer
