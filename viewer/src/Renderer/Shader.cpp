// Renderer/Shader.cpp

#include "Renderer/Shader.h"

#include <glad/glad.h>
#include <glm/gtc/type_ptr.hpp>

#include <stdexcept>
#include <vector>

namespace Aeolion::Viewer {

namespace {

[[nodiscard]] GLuint CompileStage(GLenum stage, const std::string& source) {
    GLuint shader = glCreateShader(stage);
    const char* src = source.c_str();
    glShaderSource(shader, 1, &src, nullptr);
    glCompileShader(shader);

    GLint ok = GL_FALSE;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        GLint length = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
        std::vector<char> log(static_cast<std::size_t>(length) + 1, '\0');
        glGetShaderInfoLog(shader, length, nullptr, log.data());
        glDeleteShader(shader);
        throw std::runtime_error(std::string("shader compile failed: ") + log.data());
    }
    return shader;
}

} // namespace

Shader::Shader(const std::string& vertexSource, const std::string& fragmentSource) {
    GLuint vs = CompileStage(GL_VERTEX_SHADER, vertexSource);
    GLuint fs = CompileStage(GL_FRAGMENT_SHADER, fragmentSource);

    m_ProgramId = glCreateProgram();
    glAttachShader(m_ProgramId, vs);
    glAttachShader(m_ProgramId, fs);
    glLinkProgram(m_ProgramId);
    glDeleteShader(vs);
    glDeleteShader(fs);

    GLint ok = GL_FALSE;
    glGetProgramiv(m_ProgramId, GL_LINK_STATUS, &ok);
    if (!ok) {
        GLint length = 0;
        glGetProgramiv(m_ProgramId, GL_INFO_LOG_LENGTH, &length);
        std::vector<char> log(static_cast<std::size_t>(length) + 1, '\0');
        glGetProgramInfoLog(m_ProgramId, length, nullptr, log.data());
        glDeleteProgram(m_ProgramId);
        throw std::runtime_error(std::string("shader link failed: ") + log.data());
    }
}

Shader::~Shader() { glDeleteProgram(m_ProgramId); }

void Shader::Bind() const { glUseProgram(m_ProgramId); }

int Shader::UniformLocation(const std::string& name) const {
    return glGetUniformLocation(m_ProgramId, name.c_str());
}

void Shader::SetMat4(const std::string& name, const glm::mat4& value) const {
    glUniformMatrix4fv(UniformLocation(name), 1, GL_FALSE, glm::value_ptr(value));
}

void Shader::SetVec3(const std::string& name, const glm::vec3& value) const {
    glUniform3fv(UniformLocation(name), 1, glm::value_ptr(value));
}

void Shader::SetFloat(const std::string& name, float value) const {
    glUniform1f(UniformLocation(name), value);
}

} // namespace Aeolion::Viewer
