// Visualization/PropellerRenderer.cpp

#include "Visualization/PropellerRenderer.h"
#include "Visualization/ColorMap.h"

#include "Aeolion/Math/Constants.h"

#include <glad/glad.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace Aeolion::Viewer {

namespace {

// Blade/hub fills sit on the same surfaces as their outlines; offsetting the
// fill pass lets the coplanar lines win the depth test cleanly.
constexpr float PanelPolygonOffset = 1.0f;

constexpr int HubSegments = 24;   // radial facets of the hub cylinder
constexpr int DiskSegments = 96;  // segments of the tip/hub circle line loops
constexpr float HubHalfLengthFactor = 0.6f; // hub cylinder half-length = this * hub radius
constexpr float AxisOverhangFactor = 1.6f;  // rotation axis extends this * radius each way

const glm::vec3 HubColor(0.55f, 0.55f, 0.58f);
const glm::vec3 BladeOutlineColor(0.12f, 0.12f, 0.14f);
const glm::vec3 HubOutlineColor(0.30f, 0.30f, 0.33f);
const glm::vec3 DiskColor(0.35f, 0.40f, 0.50f);
const glm::vec3 AxisColor(0.90f, 0.75f, 0.30f);

// Shaders match LatticeRenderer's so the two screens shade identically: a
// two-sided Lambert panel program and a flat colored-line program.
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

// The per-station scalar the caller chose to color blades by.
[[nodiscard]] double FieldValue(const BEMT::StationResult& s, PropFieldMode field) {
    switch (field) {
    case PropFieldMode::Alpha:       return s.alphaDeg;
    case PropFieldMode::LoadingdTdr: return s.dT_dr;
    case PropFieldMode::InflowPhi:   return s.phiDeg;
    case PropFieldMode::SectionCl:   return s.cl;
    }
    return 0.0;
}

} // namespace

PropellerRenderer::PropellerRenderer()
    : m_PanelShader(PanelVertexSource, PanelFragmentSource),
      m_LineShader(LineVertexSource, LineFragmentSource) {
    m_PanelArray.SetLayout(m_PanelBuffer, {3, 3, 3}); // position, normal, color
    m_LineArray.SetLayout(m_LineBuffer, {3, 3});      // position, color
}

void PropellerRenderer::AppendLine(std::vector<float>& lines, const glm::vec3& a, const glm::vec3& b,
                                   const glm::vec3& color) const {
    lines.insert(lines.end(), {a.x, a.y, a.z, color.r, color.g, color.b,
                               b.x, b.y, b.z, color.r, color.g, color.b});
}

void PropellerRenderer::Update(const BEMT::PropGeometry& prop, const BEMT::Result& result,
                               const PropellerDisplayOptions& options) {
    const std::vector<BEMT::BladeStation>& stations = prop.Stations;
    const std::size_t ns = stations.size();

    // --- per-station scalar field, mapped onto the colormap -----------------
    std::vector<double> field(ns, 0.0);
    for (std::size_t i = 0; i < ns && i < result.Stations.size(); ++i)
        field[i] = FieldValue(result.Stations[i], options.Field);

    m_FieldMin = std::numeric_limits<double>::max();
    m_FieldMax = std::numeric_limits<double>::lowest();
    for (double v : field) {
        m_FieldMin = std::min(m_FieldMin, v);
        m_FieldMax = std::max(m_FieldMax, v);
    }
    if (ns == 0) m_FieldMin = m_FieldMax = 0.0;
    double range = m_FieldMax - m_FieldMin;
    if (range < 1e-12) range = 1.0; // flat field -> single colormap value

    std::vector<float> panelVerts;
    std::vector<float> lineVerts;

    const glm::vec3 axial(1.0f, 0.0f, 0.0f); // thrust axis, +x

    auto pushVertex = [&](const glm::vec3& pos, const glm::vec3& normal, const glm::vec3& color) {
        panelVerts.insert(panelVerts.end(),
                          {pos.x, pos.y, pos.z, normal.x, normal.y, normal.z,
                           color.r, color.g, color.b});
    };

    // --- blades -------------------------------------------------------------
    // Each blade is the same chord/twist distribution swept out at its own
    // azimuth. A blade element at radius r sits on the radial line rhat; its
    // chord lies in the plane of the axial (+x) and tangential directions,
    // tilted out of the rotor plane by the local geometric twist.
    if (options.ShowBlades && ns >= 2) {
        const int blades = std::max(prop.NBlades, 1);
        for (int b = 0; b < blades; ++b) {
            const float psi = 2.0f * std::numbers::pi_v<float> * static_cast<float>(b) /
                              static_cast<float>(blades);
            const glm::vec3 rhat(0.0f, std::cos(psi), std::sin(psi));       // radial, in y-z plane
            const glm::vec3 that(0.0f, -std::sin(psi), std::cos(psi));      // tangential = axial x rhat

            // Leading/trailing edge point of every station's chord segment.
            std::vector<glm::vec3> leading(ns), trailing(ns);
            for (std::size_t i = 0; i < ns; ++i) {
                const float beta = static_cast<float>(Math::DegToRad(stations[i].TwistDeg));
                // Twist tilts the chord from the rotor plane (tangential) toward
                // the axial direction; the chord is centered on the radial line.
                const glm::vec3 chordDir = that * std::cos(beta) + axial * std::sin(beta);
                const glm::vec3 center = rhat * static_cast<float>(stations[i].r);
                const float half = 0.5f * static_cast<float>(stations[i].Chord);
                leading[i] = center - chordDir * half;
                trailing[i] = center + chordDir * half;
            }

            for (std::size_t i = 0; i + 1 < ns; ++i) {
                const glm::vec3& l0 = leading[i];
                const glm::vec3& t0 = trailing[i];
                const glm::vec3& l1 = leading[i + 1];
                const glm::vec3& t1 = trailing[i + 1];

                glm::vec3 normal = glm::cross(t0 - l0, l1 - l0);
                float nlen = glm::length(normal);
                normal = (nlen > 1e-8f) ? normal / nlen : glm::vec3(1.0f, 0.0f, 0.0f);

                const float ta = static_cast<float>((field[i] - m_FieldMin) / range);
                const glm::vec3 color = Viridis(ta);

                pushVertex(l0, normal, color); pushVertex(t0, normal, color); pushVertex(t1, normal, color);
                pushVertex(l0, normal, color); pushVertex(t1, normal, color); pushVertex(l1, normal, color);

                if (options.ShowWireframe) {
                    AppendLine(lineVerts, l0, t0, BladeOutlineColor);
                    AppendLine(lineVerts, t0, t1, BladeOutlineColor);
                    AppendLine(lineVerts, t1, l1, BladeOutlineColor);
                    AppendLine(lineVerts, l1, l0, BladeOutlineColor);
                }
            }
        }
    }

    // --- hub cylinder -------------------------------------------------------
    if (options.ShowHub && prop.HubRadius > 0.0) {
        const float hubR = static_cast<float>(prop.HubRadius);
        const float halfLen = hubR * HubHalfLengthFactor;
        for (int i = 0; i < HubSegments; ++i) {
            const float a0 = 2.0f * std::numbers::pi_v<float> * static_cast<float>(i) / HubSegments;
            const float a1 = 2.0f * std::numbers::pi_v<float> * static_cast<float>(i + 1) / HubSegments;
            const glm::vec3 r0(0.0f, std::cos(a0), std::sin(a0));
            const glm::vec3 r1(0.0f, std::cos(a1), std::sin(a1));
            const glm::vec3 front = axial * halfLen;
            const glm::vec3 back = -axial * halfLen;

            const glm::vec3 p00 = r0 * hubR + back;
            const glm::vec3 p01 = r0 * hubR + front;
            const glm::vec3 p10 = r1 * hubR + back;
            const glm::vec3 p11 = r1 * hubR + front;

            pushVertex(p00, r0, HubColor); pushVertex(p01, r0, HubColor); pushVertex(p11, r1, HubColor);
            pushVertex(p00, r0, HubColor); pushVertex(p11, r1, HubColor); pushVertex(p10, r1, HubColor);

            // End caps (front and back), fanned to the axis point.
            pushVertex(front, axial, HubColor); pushVertex(p01, axial, HubColor); pushVertex(p11, axial, HubColor);
            pushVertex(back, -axial, HubColor); pushVertex(p10, -axial, HubColor); pushVertex(p00, -axial, HubColor);
        }
    }

    // --- actuator-disk tip circle, hub circle, and rotation axis ------------
    if (options.ShowDisk && prop.Radius > 0.0) {
        const float tipR = static_cast<float>(prop.Radius);
        const float hubR = static_cast<float>(prop.HubRadius);
        glm::vec3 prevTip(0.0f), prevHub(0.0f);
        for (int i = 0; i <= DiskSegments; ++i) {
            const float a = 2.0f * std::numbers::pi_v<float> * static_cast<float>(i) / DiskSegments;
            const glm::vec3 rhat(0.0f, std::cos(a), std::sin(a));
            const glm::vec3 tip = rhat * tipR;
            const glm::vec3 hub = rhat * hubR;
            if (i > 0) {
                AppendLine(lineVerts, prevTip, tip, DiskColor);
                if (hubR > 0.0f) AppendLine(lineVerts, prevHub, hub, HubOutlineColor);
            }
            prevTip = tip;
            prevHub = hub;
        }
        const float overhang = tipR * AxisOverhangFactor;
        AppendLine(lineVerts, -axial * overhang, axial * overhang, AxisColor);
    }

    m_PanelVertexCount = panelVerts.size() / 9;
    m_LineVertexCount = lineVerts.size() / 6;
    m_PanelBuffer.SetData(panelVerts.data(), panelVerts.size() * sizeof(float));
    m_LineBuffer.SetData(lineVerts.data(), lineVerts.size() * sizeof(float));
}

void PropellerRenderer::Draw(const glm::mat4& viewProjection) const {
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
