// Visualization/PropellerRenderer.h
//
// Turns a BEMT propeller (geometry + solved radial distribution) into GPU
// meshes and draws them: the blades as twisted, tapered surfaces colored by
// a chosen blade-element scalar field through the viridis colormap, plus the
// hub, the actuator-disk tip circle, and the rotation axis.
//
// The propeller "screen" counterpart to LatticeRenderer: same aero-to-GL
// boundary, same two-shader (shaded triangles + colored lines) approach, but
// posed on a BEMT::Result instead of a Solver::SolveResult.

#pragma once

#include "Renderer/Buffer.h"
#include "Renderer/Shader.h"
#include "Renderer/VertexArray.h"

#include "Aeolion/BEMT/BEMT.h"
#include "Aeolion/BEMT/PropGeometry.h"

#include <glm/glm.hpp>

#include <cstddef>
#include <vector>

namespace Aeolion::Viewer {

// Which per-station BEMT quantity colors the blade surface.
enum class PropFieldMode {
    Alpha = 0,    // blade-element angle of attack [deg]
    LoadingdTdr,  // thrust grading dT/dr [N/m]
    InflowPhi,    // inflow angle phi [deg]
    SectionCl,    // section lift coefficient
};

struct PropellerDisplayOptions {
    PropFieldMode Field = PropFieldMode::Alpha;
    bool ShowBlades = true;
    bool ShowHub = true;
    bool ShowDisk = true;      // tip circle, hub circle, and rotation axis
    bool ShowWireframe = true; // blade panel outlines
};

// The propeller is drawn about the +x axis (thrust along +x, disk in the
// y-z plane), matching the solver's body axes (x aft, y right, z up).
class PropellerRenderer {
public:
    PropellerRenderer();

    // Rebuild all vertex data from a propeller and its solved radial
    // distribution. `result.Stations` is expected to align 1:1 with
    // `prop.Stations` (as BEMT::Solve produces); any shortfall just leaves
    // those blade strips at the field minimum.
    void Update(const BEMT::PropGeometry& prop, const BEMT::Result& result,
                const PropellerDisplayOptions& options);

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
