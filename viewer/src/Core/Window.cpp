// Core/Window.cpp

#include "Core/Window.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <stdexcept>

namespace Aeolion::Viewer {

Window::Window(const std::string& title, int width, int height) {
    if (!glfwInit())
        throw std::runtime_error("glfwInit failed");

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    m_Window = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
    if (!m_Window) {
        glfwTerminate();
        throw std::runtime_error("glfwCreateWindow failed");
    }

    glfwMakeContextCurrent(m_Window);
    glfwSwapInterval(1); // vsync

    if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
        glfwDestroyWindow(m_Window);
        glfwTerminate();
        throw std::runtime_error("gladLoadGLLoader failed");
    }

    // Installed BEFORE ImGui's backend init so ImGui chains to it.
    glfwSetWindowUserPointer(m_Window, this);
    glfwSetScrollCallback(m_Window, [](GLFWwindow* window, double /*dx*/, double dy) {
        auto* self = static_cast<Window*>(glfwGetWindowUserPointer(window));
        self->m_ScrollDelta += dy;
    });
}

Window::~Window() {
    glfwDestroyWindow(m_Window);
    glfwTerminate();
}

bool Window::ShouldClose() const { return glfwWindowShouldClose(m_Window); }
void Window::PollEvents() { glfwPollEvents(); }
void Window::SwapBuffers() { glfwSwapBuffers(m_Window); }

int Window::FramebufferWidth() const {
    int w = 0, h = 0;
    glfwGetFramebufferSize(m_Window, &w, &h);
    return w;
}

int Window::FramebufferHeight() const {
    int w = 0, h = 0;
    glfwGetFramebufferSize(m_Window, &w, &h);
    return h;
}

double Window::TakeScrollDelta() {
    double d = m_ScrollDelta;
    m_ScrollDelta = 0.0;
    return d;
}

} // namespace Aeolion::Viewer
