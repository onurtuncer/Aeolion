// Renderer/OrbitCamera.h
//
// Turntable camera orbiting a target point, z-up to match the solver's
// body axes (x aft, y right, z up). Rotate/pan/zoom are driven with mouse
// deltas from the frame loop; ViewProjection() feeds the shaders.

#pragma once

#include <glm/glm.hpp>

namespace Aeolion::Viewer {

class OrbitCamera {
public:
    void SetViewport(int width, int height);
    void SetTarget(const glm::vec3& target) { m_Target = target; }
    void SetDistance(float distance);

    void Rotate(float dYawDeg, float dPitchDeg);
    void Pan(float dxPixels, float dyPixels, int viewportHeight);
    void Zoom(float scrollSteps);

    [[nodiscard]] glm::vec3 Position() const;
    [[nodiscard]] glm::mat4 ViewProjection() const;

private:
    glm::vec3 m_Target{0.0f};
    float m_YawDeg = -135.0f;
    float m_PitchDeg = 25.0f;
    float m_Distance = 5.0f;
    float m_AspectRatio = 16.0f / 9.0f;
};

} // namespace Aeolion::Viewer
