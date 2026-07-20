// Renderer/VertexArray.cpp

#include "Renderer/VertexArray.h"
#include "Renderer/Buffer.h"

#include <glad/glad.h>

namespace Aeolion::Viewer {

VertexArray::VertexArray() { glGenVertexArrays(1, &m_ArrayId); }
VertexArray::~VertexArray() { glDeleteVertexArrays(1, &m_ArrayId); }

void VertexArray::Bind() const { glBindVertexArray(m_ArrayId); }

void VertexArray::SetLayout(const VertexBuffer& buffer, std::initializer_list<int> componentCounts) {
    Bind();
    buffer.Bind();

    int floatsPerVertex = 0;
    for (int count : componentCounts) floatsPerVertex += count;
    GLsizei stride = floatsPerVertex * static_cast<GLsizei>(sizeof(float));

    GLuint location = 0;
    std::size_t offsetFloats = 0;
    for (int count : componentCounts) {
        glEnableVertexAttribArray(location);
        glVertexAttribPointer(location, count, GL_FLOAT, GL_FALSE, stride,
                              reinterpret_cast<const void*>(offsetFloats * sizeof(float)));
        ++location;
        offsetFloats += static_cast<std::size_t>(count);
    }
}

} // namespace Aeolion::Viewer
