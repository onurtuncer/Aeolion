// Core/Window.h
//
// GLFW window + OpenGL context bootstrap. Owns the one window the viewer
// runs in: creates the GL 3.3 core context, loads GL entry points through
// glad, and accumulates scroll input for the frame loop to consume.

#pragma once

#include <string>

struct GLFWwindow;

namespace Aeolion::Viewer {

class Window {
public:
    Window(const std::string& title, int width, int height);
    ~Window();

    Window(const Window&) = delete;
    Window& operator=(const Window&) = delete;

    [[nodiscard]] bool ShouldClose() const;
    void PollEvents();
    void SwapBuffers();

    [[nodiscard]] GLFWwindow* Handle() const { return m_Window; }
    [[nodiscard]] int FramebufferWidth() const;
    [[nodiscard]] int FramebufferHeight() const;

    // Wheel movement accumulated since the last call; calling it clears it.
    [[nodiscard]] double TakeScrollDelta();

private:
    GLFWwindow* m_Window = nullptr;
    double m_ScrollDelta = 0.0;
};

} // namespace Aeolion::Viewer
