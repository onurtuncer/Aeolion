// Visualization/LatticeRenderer.cpp

#include "Visualization/LatticeRenderer.h"
#include "Visualization/ColorMap.h"

#include <glad/glad.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace Aeolion::Viewer {

namespace {

// Panel quads sit exactly on the lattice surface; offsetting the fill pass
// lets the coplanar outline/vortex lines win the depth test cleanly.
constexpr float PanelPolygonOffset = 1.0f;

constexpr float QuarterChord = 0.25f;
constexpr float ThreeQuarterChord = 0.75f;

// Overlay scale factors, all relative to the wing span.
constexpr float NormalGlyphSpanFactor = 0.04f;
constexpr float GridHalfExtentSpanFactor = 0.8f;
constexpr float GridDropSpanFactor = 0.3f;   // grid plane sits this far below the lattice
constexpr int   GridLineCount = 17;          // lines per direction (odd -> line through center)
constexpr float AxisLengthSpanFactor = 0.3f;

const glm::vec3 OutlineColor(0.42f, 0.44f, 0.48f);
const glm::vec3 BoundVortexColor(0.95f, 0.95f, 0.98f);
const glm::vec3 NormalGlyphColor(1.0f, 0.62f, 0.25f);
const glm::vec3 GridColor(0.20f, 0.21f, 0.24f);
const glm::vec3 AxisXColor(0.85f, 0.30f, 0.30f);
const glm::vec3 AxisYColor(0.35f, 0.80f, 0.35f);
const glm::vec3 AxisZColor(0.35f, 0.55f, 0.95f);

const char* PanelVertexSource = R"(#version 330 core
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec3 aColor;
uniform mat4 u_ViewProjection;
out vec3 vNormal;
out vec3 vColor;
void main() {
    vNormal = aNormal;
    vColor = aColor;
    gl_Position = u_ViewProjection * vec4(aPosition, 1.0);
})";

const char* PanelFragmentSource = R"(#version 330 core
in vec3 vNormal;
in vec3 vColor;
uniform vec3 u_LightDirection;
out vec4 FragColor;
void main() {
    // Two-sided Lambert-ish shading; the field color stays dominant.
    float incidence = abs(dot(normalize(vNormal), normalize(u_LightDirection)));
    float shade = 0.55 + 0.45 * incidence;
    FragColor = vec4(vColor * shade, 1.0);
})";

const char* LineVertexSource = R"(#version 330 core
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aColor;
uniform mat4 u_ViewProjection;
out vec3 vColor;
void main() {
    vColor = aColor;
    gl_Position = u_ViewProjection * vec4(aPosition, 1.0);
})";

const char* LineFragmentSource = R"(#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
    FragColor = vec4(vColor, 1.0);
})";

[[nodiscard]] glm::vec3 ToGlm(const Math::Vec3& v) {
    return {static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z)};
}

} // namespace

LatticeRenderer::LatticeRenderer()
    : m_PanelShader(PanelVertexSource, PanelFragmentSource),
      m_LineShader(LineVertexSource, LineFragmentSource) {
    m_PanelArray.SetLayout(m_PanelBuffer, {3, 3, 3}); // position, normal, color
    m_LineArray.SetLayout(m_LineBuffer, {3, 3});      // position, color
}

void LatticeRenderer::AppendLine(std::vector<float>& lines, const glm::vec3& a, const glm::vec3& b,
                                 const glm::vec3& color) const {
    lines.insert(lines.end(), {a.x, a.y, a.z, color.r, color.g, color.b,
                               b.x, b.y, b.z, color.r, color.g, color.b});
}

void LatticeRenderer::Update(const std::vector<Solver::Panel>& panels, const Solver::SolveResult& result,
                             const LatticeDisplayOptions& options, double span) {
    const std::size_t n = panels.size();

    // --- scalar field per panel (indexed like `panels`, not by station order)
    // Circulation is per panel, but cl and lift-per-span are per SPANWISE
    // STATION -- one value for a whole chordwise stack. Those are broadcast
    // back over every panel sharing the station's strip so a multi-row
    // lattice colors as coherent spanwise bands rather than leaving all but
    // the stack's first panel black.
    std::vector<double> field(n, 0.0);
    auto stripKey = [&](std::size_t panelIndex) {
        return panels[panelIndex].StripIndex >= 0 ? panels[panelIndex].StripIndex
                                                  : -static_cast<int>(panelIndex) - 1;
    };
    auto broadcastByStrip = [&](auto valueOf) {
        std::map<int, double> byStrip;
        for (const auto& station : result.Stations) {
            if (station.PanelIndex < 0 || static_cast<std::size_t>(station.PanelIndex) >= n) continue;
            byStrip[stripKey(static_cast<std::size_t>(station.PanelIndex))] = valueOf(station);
        }
        for (std::size_t i = 0; i < n; ++i) {
            const auto found = byStrip.find(stripKey(i));
            if (found != byStrip.end()) field[i] = found->second;
        }
    };

    switch (options.Field) {
    case FieldMode::Gamma:
        for (std::size_t i = 0; i < n && i < result.gamma.size(); ++i) field[i] = result.gamma[i];
        break;
    case FieldMode::LocalCl:
        broadcastByStrip([](const Solver::StationResult& s) { return s.cl_local; });
        break;
    case FieldMode::LiftPerSpan:
        broadcastByStrip([](const Solver::StationResult& s) { return s.LiftPerSpan; });
        break;
    }

    m_FieldMin = std::numeric_limits<double>::max();
    m_FieldMax = std::numeric_limits<double>::lowest();
    for (double v : field) {
        m_FieldMin = std::min(m_FieldMin, v);
        m_FieldMax = std::max(m_FieldMax, v);
    }
    if (n == 0) m_FieldMin = m_FieldMax = 0.0;
    double range = m_FieldMax - m_FieldMin;
    if (range < 1e-12) range = 1.0; // flat field -> single colormap value

    // --- panel quads --------------------------------------------------------
    // A panel stores its bound (quarter-chord) segment A-B, area, and width;
    // the chordwise footprint is reconstructed from the local chord
    // c = Area/width along the control-point direction, which already
    // carries the geometric twist: quarter chord ahead, three quarters aft.
    std::vector<float> panelVerts;
    panelVerts.reserve(n * 6 * 9);
    std::vector<float> lineVerts;

    for (std::size_t i = 0; i < n; ++i) {
        const Solver::Panel& p = panels[i];
        // The drawn quad is the panel's planform footprint, so it is built
        // from the projected area -- Area is the curved surface and would
        // over-length the chord (see Solver::Panel).
        double chord = (p.SpanwiseWidth > 0.0) ? p.PlanformArea / p.SpanwiseWidth : 0.0;

        glm::vec3 a = ToGlm(p.A);
        glm::vec3 b = ToGlm(p.B);
        glm::vec3 chordDir = ToGlm((p.ControlPoint - (p.A + p.B) * 0.5).Normalized());
        if (glm::dot(chordDir, chordDir) < 0.5f) chordDir = glm::vec3(1, 0, 0); // degenerate: fall back to +x

        float ahead = static_cast<float>(QuarterChord * chord);
        float aft = static_cast<float>(ThreeQuarterChord * chord);

        glm::vec3 leadingA = a - chordDir * ahead;
        glm::vec3 leadingB = b - chordDir * ahead;
        glm::vec3 trailingA = a + chordDir * aft;
        glm::vec3 trailingB = b + chordDir * aft;

        glm::vec3 normal = ToGlm(p.Normal);
        float t = static_cast<float>((field[i] - m_FieldMin) / range);
        glm::vec3 color = Viridis(t);

        auto pushVertex = [&](const glm::vec3& pos) {
            panelVerts.insert(panelVerts.end(),
                              {pos.x, pos.y, pos.z, normal.x, normal.y, normal.z,
                               color.r, color.g, color.b});
        };
        pushVertex(leadingA); pushVertex(leadingB); pushVertex(trailingB);
        pushVertex(leadingA); pushVertex(trailingB); pushVertex(trailingA);

        if (options.ShowWireframe) {
            AppendLine(lineVerts, leadingA, leadingB, OutlineColor);
            AppendLine(lineVerts, leadingB, trailingB, OutlineColor);
            AppendLine(lineVerts, trailingB, trailingA, OutlineColor);
            AppendLine(lineVerts, trailingA, leadingA, OutlineColor);
            AppendLine(lineVerts, a, b, BoundVortexColor);
        }
        if (options.ShowNormals) {
            glm::vec3 cp = ToGlm(p.ControlPoint);
            float glyphLength = static_cast<float>(span) * NormalGlyphSpanFactor;
            AppendLine(lineVerts, cp, cp + normal * glyphLength, NormalGlyphColor);
        }
    }

    // --- reference grid and body axes --------------------------------------
    float s = static_cast<float>(span);
    if (options.ShowGrid && s > 0.0f) {
        float half = s * GridHalfExtentSpanFactor;
        float z = -s * GridDropSpanFactor;
        for (int i = 0; i < GridLineCount; ++i) {
            float u = -half + 2.0f * half * static_cast<float>(i) / (GridLineCount - 1);
            AppendLine(lineVerts, {u, -half, z}, {u, half, z}, GridColor);
            AppendLine(lineVerts, {-half, u, z}, {half, u, z}, GridColor);
        }
    }
    if (options.ShowAxes && s > 0.0f) {
        float len = s * AxisLengthSpanFactor;
        AppendLine(lineVerts, {0, 0, 0}, {len, 0, 0}, AxisXColor); // +x aft
        AppendLine(lineVerts, {0, 0, 0}, {0, len, 0}, AxisYColor); // +y right
        AppendLine(lineVerts, {0, 0, 0}, {0, 0, len}, AxisZColor); // +z up
    }

    m_PanelVertexCount = panelVerts.size() / 9;
    m_LineVertexCount = lineVerts.size() / 6;
    m_PanelBuffer.SetData(panelVerts.data(), panelVerts.size() * sizeof(float));
    m_LineBuffer.SetData(lineVerts.data(), lineVerts.size() * sizeof(float));
}

void LatticeRenderer::Draw(const glm::mat4& viewProjection) const {
    if (m_PanelVertexCount > 0) {
        m_PanelShader.Bind();
        m_PanelShader.SetMat4("u_ViewProjection", viewProjection);
        m_PanelShader.SetVec3("u_LightDirection", glm::normalize(glm::vec3(-0.4f, 0.3f, -1.0f)));
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(PanelPolygonOffset, PanelPolygonOffset);
        m_PanelArray.Bind();
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(m_PanelVertexCount));
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
    if (m_LineVertexCount > 0) {
        m_LineShader.Bind();
        m_LineShader.SetMat4("u_ViewProjection", viewProjection);
        m_LineArray.Bind();
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_LineVertexCount));
    }
}

} // namespace Aeolion::Viewer
