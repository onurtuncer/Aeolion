// Renderer/Buffer.cpp

#include "Renderer/Buffer.h"

#include <glad/glad.h>

namespace Aeolion::Viewer {

VertexBuffer::VertexBuffer() { glGenBuffers(1, &m_BufferId); }
VertexBuffer::~VertexBuffer() { glDeleteBuffers(1, &m_BufferId); }

void VertexBuffer::Bind() const { glBindBuffer(GL_ARRAY_BUFFER, m_BufferId); }

void VertexBuffer::SetData(const void* data, std::size_t sizeBytes) {
    Bind();
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sizeBytes), data, GL_DYNAMIC_DRAW);
}

} // namespace Aeolion::Viewer
