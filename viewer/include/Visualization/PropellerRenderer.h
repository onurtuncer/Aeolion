// Visualization/PropellerRenderer.h
//
// Turns a Geometry::Propeller into GPU meshes and draws them: the blades as
// twisted, tapered surfaces colored by a chosen per-station geometric field
// through the viridis colormap, plus the hub, the actuator-disk tip circle,
// and the rotation axis.
//
// The propeller "screen" counterpart to LatticeRenderer: same aero-to-GL
// boundary, same two-shader (shaded triangles + colored lines) approach --
// but posed on geometry alone. The propeller's aerodynamic solve moved out
// with BEMT (now a separate project this repo does not depend on); when a
// calculation method lands here, its per-station results become additional
// PropFieldMode entries and an extra Update() input, nothing more.

#pragma once

#include "Renderer/Buffer.h"
#include "Renderer/Shader.h"
#include "Renderer/VertexArray.h"

#include "Aeolion/Geometry/Propeller.h"

#include <glm/glm.hpp>

#include <cstddef>
#include <vector>

namespace Aeolion::Viewer {

// Which per-station geometric quantity colors the blade surface.
enum class PropFieldMode {
    Twist = 0, // local geometric pitch [deg]
    Chord,     // local chord [m]
};

struct PropellerDisplayOptions {
    PropFieldMode Field = PropFieldMode::Twist;
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

    // Rebuild all vertex data from a propeller's geometry.
    void Update(const Geometry::Propeller& prop, const PropellerDisplayOptions& options);

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
