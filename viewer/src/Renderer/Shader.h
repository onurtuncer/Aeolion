// Renderer/Shader.h
//
// GLSL program wrapper: compile/link from in-memory source strings (the
// viewer embeds its shaders as string literals -- no runtime asset paths)
// plus the handful of uniform setters the visualization needs.

#pragma once

#include <glm/glm.hpp>

#include <cstdint>
#include <string>

namespace Aeolion::Viewer {

class Shader {
public:
    Shader(const std::string& vertexSource, const std::string& fragmentSource);
    ~Shader();

    Shader(const Shader&) = delete;
    Shader& operator=(const Shader&) = delete;

    void Bind() const;

    void SetMat4(const std::string& name, const glm::mat4& value) const;
    void SetVec3(const std::string& name, const glm::vec3& value) const;
    void SetFloat(const std::string& name, float value) const;

private:
    [[nodiscard]] int UniformLocation(const std::string& name) const;

    std::uint32_t m_ProgramId = 0;
};

} // namespace Aeolion::Viewer
