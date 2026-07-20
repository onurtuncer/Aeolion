# Aeolion viewer: interactive OpenGL visualizer for the Aeolion solver.
#
# Hazel-style layering, trimmed to what a single-backend visualizer needs:
#   Core/          Application, Window (GLFW + glad bootstrap)
#   Renderer/      Shader, Buffer, VertexArray, OrbitCamera (GL, no aero)
#   Visualization/ LatticeRenderer + ColorMap (solver types -> GPU meshes)
#
# The aero library stays GL-free; GLFW/glad/ImGui/glm are viewer-only deps.

find_package(glfw3 CONFIG REQUIRED)
find_package(glad CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)
find_package(OpenGL REQUIRED)

add_executable(aeolion_viewer
    src/Main.cpp
    src/Core/Application.cpp
    src/Core/Window.cpp
    src/Renderer/Buffer.cpp
    src/Renderer/OrbitCamera.cpp
    src/Renderer/Shader.cpp
    src/Renderer/VertexArray.cpp
    src/Visualization/LatticeRenderer.cpp
)

target_include_directories(aeolion_viewer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(aeolion_viewer PRIVATE
    aeolion
    glfw
    glad::glad
    glm::glm
    imgui::imgui
    OpenGL::GL
)
