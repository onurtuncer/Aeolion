// Renderer/OrbitCamera.cpp

#include "Renderer/OrbitCamera.h"

#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <cmath>

namespace Aeolion::Viewer {

namespace {
constexpr float MinDistance = 0.05f;
constexpr float MaxDistance = 500.0f;
constexpr float PitchLimitDeg = 89.0f;      // keep away from the up-vector singularity
constexpr float ZoomStepFactor = 0.9f;      // multiplicative zoom per wheel notch
constexpr float FovYDeg = 45.0f;
constexpr float PanWorldPerScreen = 2.0f;   // world units panned across one viewport height at unit distance (~2*tan(fov/2))
} // namespace

void OrbitCamera::SetViewport(int width, int height) {
    if (width > 0 && height > 0)
        m_AspectRatio = static_cast<float>(width) / static_cast<float>(height);
}

void OrbitCamera::SetDistance(float distance) {
    m_Distance = std::clamp(distance, MinDistance, MaxDistance);
}

void OrbitCamera::Rotate(float dYawDeg, float dPitchDeg) {
    m_YawDeg += dYawDeg;
    m_PitchDeg = std::clamp(m_PitchDeg + dPitchDeg, -PitchLimitDeg, PitchLimitDeg);
}

void OrbitCamera::Pan(float dxPixels, float dyPixels, int viewportHeight) {
    if (viewportHeight <= 0) return;
    float worldPerPixel = PanWorldPerScreen * m_Distance / static_cast<float>(viewportHeight);

    glm::vec3 forward = glm::normalize(m_Target - Position());
    glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0, 0, 1)));
    glm::vec3 up = glm::cross(right, forward);

    m_Target = m_Target + (right * -dxPixels + up * dyPixels) * worldPerPixel;
}

void OrbitCamera::Zoom(float scrollSteps) {
    SetDistance(m_Distance * std::pow(ZoomStepFactor, scrollSteps));
}

glm::vec3 OrbitCamera::Position() const {
    float yaw = glm::radians(m_YawDeg);
    float pitch = glm::radians(m_PitchDeg);
    glm::vec3 offset(std::cos(pitch) * std::cos(yaw),
                     std::cos(pitch) * std::sin(yaw),
                     std::sin(pitch));
    return m_Target + offset * m_Distance;
}

glm::mat4 OrbitCamera::ViewProjection() const {
    glm::mat4 view = glm::lookAt(Position(), m_Target, glm::vec3(0, 0, 1));
    // Near/far track the orbit distance so depth precision follows the scene scale.
    float nearPlane = std::max(0.001f, m_Distance * 0.01f);
    float farPlane = m_Distance * 100.0f;
    glm::mat4 proj = glm::perspective(glm::radians(FovYDeg), m_AspectRatio, nearPlane, farPlane);
    return proj * view;
}

} // namespace Aeolion::Viewer
