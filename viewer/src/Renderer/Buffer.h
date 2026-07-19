// Renderer/Buffer.h
//
// GL_ARRAY_BUFFER wrapper. The visualization meshes are rebuilt whenever
// the solve changes, so buffers are GL_DYNAMIC_DRAW and re-uploaded whole
// via SetData (orphaning realloc -- fine at lattice sizes).

#pragma once

#include <cstddef>
#include <cstdint>

namespace Aeolion::Viewer {

class VertexBuffer {
public:
    VertexBuffer();
    ~VertexBuffer();

    VertexBuffer(const VertexBuffer&) = delete;
    VertexBuffer& operator=(const VertexBuffer&) = delete;

    void Bind() const;
    void SetData(const void* data, std::size_t sizeBytes);

private:
    std::uint32_t m_BufferId = 0;
};

} // namespace Aeolion::Viewer
