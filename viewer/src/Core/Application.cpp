// Core/Application.cpp

#include "Core/Application.h"
#include "Visualization/ColorMap.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <algorithm>
#include <cfloat>
#include <numbers>
#include <string>
#include <vector>

namespace Aeolion::Viewer {

namespace {

constexpr int DefaultWindowWidth = 1600;
constexpr int DefaultWindowHeight = 900;
constexpr float RotateDegPerPixel = 0.25f;
constexpr float CameraFrameDistanceFactor = 1.4f; // orbit distance = this * span when framing

// ImGui sliders for the solver's double-typed parameters.
bool SliderD(const char* label, double* value, double min, double max, const char* format = "%.3f") {
    return ImGui::SliderScalar(label, ImGuiDataType_Double, value, &min, &max, format);
}

void ReadoutRow(const char* label, double value, const char* format = "%.5f") {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::TextUnformatted(label);
    ImGui::TableSetColumnIndex(1);
    ImGui::Text(format, value);
}

// Horizontal viridis gradient with the current field range at its ends.
void ColorBar(double minValue, double maxValue) {
    ImDrawList* drawList = ImGui::GetWindowDrawList();
    ImVec2 origin = ImGui::GetCursorScreenPos();
    float width = ImGui::GetContentRegionAvail().x;
    const float height = 14.0f;
    const int strips = 48;
    for (int i = 0; i < strips; ++i) {
        float t0 = static_cast<float>(i) / strips;
        float t1 = static_cast<float>(i + 1) / strips;
        glm::vec3 c = Viridis((t0 + t1) * 0.5f);
        drawList->AddRectFilled(ImVec2(origin.x + t0 * width, origin.y),
                                ImVec2(origin.x + t1 * width, origin.y + height),
                                IM_COL32(static_cast<int>(c.r * 255), static_cast<int>(c.g * 255),
                                         static_cast<int>(c.b * 255), 255));
    }
    ImGui::Dummy(ImVec2(width, height + 2.0f));
    ImGui::Text("%.4g", minValue);
    ImGui::SameLine(std::max(0.0f, width - 60.0f));
    ImGui::Text("%.4g", maxValue);
}

} // namespace

Application::Application()
    : m_Window("Aeolion Viewer", DefaultWindowWidth, DefaultWindowHeight) {
    // Sensible default study: a moderately swept, tapered, washed-out wing.
    m_Wing.Span = 4.0;
    m_Wing.RootChord = 0.6;
    m_Wing.TipChord = 0.35;
    m_Wing.SweepQuarterChordDeg = 15.0;
    m_Wing.DihedralDeg = 3.0;
    m_Wing.TwistTipDeg = -2.0;
    m_Wing.NPanelsSemiSpan = 16;
    m_Wing.CosineSpacing = true;

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(m_Window.Handle(), true);
    ImGui_ImplOpenGL3_Init("#version 330");

    FrameView();
}

Application::~Application() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void Application::FrameView() {
    // Center on the mid-chord of the root section; the lattice grows aft (+x).
    m_Camera.SetTarget(glm::vec3(static_cast<float>(m_Wing.RootChord) * 0.25f, 0.0f, 0.0f));
    m_Camera.SetDistance(static_cast<float>(m_Wing.Span) * CameraFrameDistanceFactor);
}

void Application::Resolve() {
    m_Panels = VLM::BuildWing(m_Wing);
    for (auto& panel : m_Panels) panel.Surface = "wing";

    double area = 0.0;
    for (const auto& panel : m_Panels) area += panel.Area;

    VLM::ReferenceGeometry ref;
    ref.Area = area;
    ref.Span = m_Wing.Span;
    ref.Chord = area / m_Wing.Span;

    double trail = (m_Wing.TrailLength > 0.0)
                       ? m_Wing.TrailLength
                       : VLM::DefaultTrailSpanFactor * m_Wing.Span;

    m_Result = VLM::Solve(m_Panels, m_Freestream, ref, trail);
    m_Derivatives.reset(); // stale for the new geometry/condition; recompute on demand
}

void Application::HandleCameraInput() {
    GLFWwindow* window = m_Window.Handle();
    double mouseX = 0.0, mouseY = 0.0;
    glfwGetCursorPos(window, &mouseX, &mouseY);
    float dx = static_cast<float>(mouseX - m_LastMouseX);
    float dy = static_cast<float>(mouseY - m_LastMouseY);
    m_LastMouseX = mouseX;
    m_LastMouseY = mouseY;

    double scroll = m_Window.TakeScrollDelta();

    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse) return;

    bool rotating = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    bool panning = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS ||
                   glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;

    if (rotating)
        m_Camera.Rotate(dx * RotateDegPerPixel, -dy * RotateDegPerPixel);
    else if (panning)
        m_Camera.Pan(dx, dy, m_Window.FramebufferHeight());
    if (scroll != 0.0)
        m_Camera.Zoom(static_cast<float>(scroll));
}

void Application::DrawUI() {
    // --- controls -----------------------------------------------------------
    ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(360, 560), ImGuiCond_FirstUseEver);
    ImGui::Begin("Controls");

    if (ImGui::CollapsingHeader("Planform", ImGuiTreeNodeFlags_DefaultOpen)) {
        m_Dirty |= SliderD("Span [m]", &m_Wing.Span, 0.5, 12.0);
        m_Dirty |= SliderD("Root chord [m]", &m_Wing.RootChord, 0.05, 2.0);
        m_Dirty |= SliderD("Tip chord [m]", &m_Wing.TipChord, 0.05, 2.0);
        m_Dirty |= SliderD("Sweep c/4 [deg]", &m_Wing.SweepQuarterChordDeg, -40.0, 40.0, "%.1f");
        m_Dirty |= SliderD("Dihedral [deg]", &m_Wing.DihedralDeg, -15.0, 15.0, "%.1f");
        m_Dirty |= SliderD("Tip twist [deg]", &m_Wing.TwistTipDeg, -10.0, 10.0, "%.1f");
        m_Dirty |= ImGui::SliderInt("Panels / semi-span", &m_Wing.NPanelsSemiSpan, 4, 60);
        m_Dirty |= ImGui::Checkbox("Cosine spacing", &m_Wing.CosineSpacing);
    }

    if (ImGui::CollapsingHeader("Freestream", ImGuiTreeNodeFlags_DefaultOpen)) {
        m_Dirty |= SliderD("Alpha [deg]", &m_Freestream.alphaDeg, -15.0, 20.0, "%.2f");
        m_Dirty |= SliderD("Beta [deg]", &m_Freestream.betaDeg, -15.0, 15.0, "%.2f");
        m_Dirty |= SliderD("Vinf [m/s]", &m_Freestream.Vinf, 1.0, 100.0, "%.1f");
        m_Dirty |= SliderD("rho [kg/m^3]", &m_Freestream.rho, 0.4, 1.4);
    }

    if (ImGui::CollapsingHeader("Display", ImGuiTreeNodeFlags_DefaultOpen)) {
        const char* fieldNames[] = {"Circulation gamma", "Sectional cl", "Lift per span"};
        int field = static_cast<int>(m_Display.Field);
        if (ImGui::Combo("Field", &field, fieldNames, IM_ARRAYSIZE(fieldNames))) {
            m_Display.Field = static_cast<FieldMode>(field);
            m_DisplayDirty = true;
        }
        m_DisplayDirty |= ImGui::Checkbox("Wireframe", &m_Display.ShowWireframe);
        m_DisplayDirty |= ImGui::Checkbox("Normals", &m_Display.ShowNormals);
        m_DisplayDirty |= ImGui::Checkbox("Grid", &m_Display.ShowGrid);
        m_DisplayDirty |= ImGui::Checkbox("Axes", &m_Display.ShowAxes);
        if (ImGui::Button("Frame view")) FrameView();
    }
    ImGui::End();

    // --- results ------------------------------------------------------------
    ImGui::SetNextWindowPos(ImVec2(static_cast<float>(m_Window.FramebufferWidth()) - 370, 10),
                            ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(360, 560), ImGuiCond_FirstUseEver);
    ImGui::Begin("Results");

    if (ImGui::BeginTable("coefficients", 2, ImGuiTableFlags_SizingStretchProp)) {
        ReadoutRow("CL", m_Result.CL);
        ReadoutRow("CDi", m_Result.CDi);
        ReadoutRow("CY", m_Result.CY);
        ReadoutRow("Cm", m_Result.Cm);
        ReadoutRow("Croll", m_Result.Croll);
        ReadoutRow("Cn", m_Result.Cn);
        ReadoutRow("L [N]", m_Result.L, "%.2f");
        ReadoutRow("Di [N]", m_Result.Di, "%.3f");
        double aspectRatio = (m_Result.ReferenceArea > 0.0)
                                 ? m_Wing.Span * m_Wing.Span / m_Result.ReferenceArea
                                 : 0.0;
        ReadoutRow("AR", aspectRatio, "%.2f");
        if (m_Result.CDi > 1e-9) {
            double e = m_Result.CL * m_Result.CL /
                       (std::numbers::pi * aspectRatio * m_Result.CDi);
            ReadoutRow("Span efficiency e", e, "%.3f");
        }
        ImGui::EndTable();
    }

    ImGui::SeparatorText("Field range");
    ColorBar(m_Lattice.FieldMin(), m_Lattice.FieldMax());

    ImGui::SeparatorText("Spanwise cl");
    if (!m_Result.Stations.empty()) {
        std::vector<float> cl(m_Result.Stations.size());
        for (std::size_t i = 0; i < cl.size(); ++i)
            cl[i] = static_cast<float>(m_Result.Stations[i].cl_local);
        ImGui::PlotLines("##cl_span", cl.data(), static_cast<int>(cl.size()), 0,
                         "cl(y), tip-to-tip", FLT_MAX, FLT_MAX, ImVec2(-1, 90));
    }

    ImGui::SeparatorText("Stability derivatives");
    if (ImGui::Button("Compute derivatives")) {
        double area = 0.0;
        for (const auto& panel : m_Panels) area += panel.Area;
        VLM::ReferenceGeometry ref;
        ref.Area = area;
        ref.Span = m_Wing.Span;
        ref.Chord = area / m_Wing.Span;
        double trail = (m_Wing.TrailLength > 0.0)
                           ? m_Wing.TrailLength
                           : VLM::DefaultTrailSpanFactor * m_Wing.Span;
        m_Derivatives = VLM::ComputeDerivatives(m_Panels, m_Freestream, ref, trail);
    }
    if (m_Derivatives) {
        if (ImGui::BeginTable("derivatives", 2, ImGuiTableFlags_SizingStretchProp)) {
            ReadoutRow("CL_alpha [/rad]", m_Derivatives->CL_alpha, "%.4f");
            ReadoutRow("Cm_alpha [/rad]", m_Derivatives->Cm_alpha, "%.4f");
            ReadoutRow("CY_beta [/rad]", m_Derivatives->CY_beta, "%.4f");
            ReadoutRow("Croll_beta [/rad]", m_Derivatives->Croll_beta, "%.4f");
            ReadoutRow("Cn_beta [/rad]", m_Derivatives->Cn_beta, "%.4f");
            ReadoutRow("CL_q (nd)", m_Derivatives->CL_q_nd, "%.4f");
            ReadoutRow("Cm_q (nd)", m_Derivatives->Cm_q_nd, "%.4f");
            ReadoutRow("Croll_p (nd)", m_Derivatives->Croll_p_nd, "%.4f");
            ReadoutRow("Cn_r (nd)", m_Derivatives->Cn_r_nd, "%.4f");
            ImGui::EndTable();
        }
    }
    ImGui::End();
}

void Application::Run(int maxFrames) {
    int frame = 0;
    while (!m_Window.ShouldClose()) {
        if (maxFrames >= 0 && frame++ >= maxFrames) break;
        m_Window.PollEvents();

        HandleCameraInput();

        if (m_Dirty) {
            Resolve();
            m_Dirty = false;
            m_DisplayDirty = true;
        }
        if (m_DisplayDirty) {
            m_Lattice.Update(m_Panels, m_Result, m_Display, m_Wing.Span);
            m_DisplayDirty = false;
        }

        int width = m_Window.FramebufferWidth();
        int height = m_Window.FramebufferHeight();
        m_Camera.SetViewport(width, height);

        glViewport(0, 0, width, height);
        glClearColor(0.09f, 0.10f, 0.13f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        m_Lattice.Draw(m_Camera.ViewProjection());

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        DrawUI();
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        m_Window.SwapBuffers();
    }
}

} // namespace Aeolion::Viewer
