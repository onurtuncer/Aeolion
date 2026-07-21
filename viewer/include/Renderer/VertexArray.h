// Renderer/VertexArray.h
//
// VAO wrapper with a deliberately small layout API: every attribute is
// float-typed, declared as a list of component counts in location order
// (e.g. {3, 3, 3} = position, normal, color). Stride and offsets are
// derived from that list.

#pragma once

#include <cstdint>
#include <initializer_list>

namespace Aeolion::Viewer {

class VertexBuffer;

class VertexArray {
public:
    VertexArray();
    ~VertexArray();

    VertexArray(const VertexArray&) = delete;
    VertexArray& operator=(const VertexArray&) = delete;

    void Bind() const;

    // Binds `buffer` into this VAO and declares its interleaved all-float
    // layout; attribute locations are assigned in list order.
    void SetLayout(const VertexBuffer& buffer, std::initializer_list<int> componentCounts);

private:
    std::uint32_t m_ArrayId = 0;
};

} // namespace Aeolion::Viewer
