// Core/Application.cpp

#include "Core/Application.h"
#include "Visualization/ColorMap.h"

#include "Aeolion/PanelBuilder/PanelBuilder.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <algorithm>
#include <cfloat>
#include <cstdio>
#include <fstream>
#include <numbers>
#include <string>
#include <vector>

namespace Aeolion::Viewer {

namespace {

constexpr int DefaultWindowWidth = 1600;
constexpr int DefaultWindowHeight = 900;
constexpr float RotateDegPerPixel = 0.25f;
constexpr float CameraFrameDistanceFactor = 1.4f; // orbit distance = this * span when framing

// Deflection slider range used when a control surface states no
// deflection_limits_deg (schema < 1.4.0, where MinDeg == MaxDeg == 0 means
// "not stated" rather than "cannot move" -- see ControlSurface.h).
constexpr double DefaultDeflectionLimitDeg = 30.0;

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

Application::Application(const std::string& geometryPath, const std::string& screenshotPath)
    : m_Window("Aeolion Viewer", DefaultWindowWidth, DefaultWindowHeight),
      m_ScreenshotPath(screenshotPath) {
    // Sensible default study: a moderately swept, tapered, washed-out wing.
    // Ignored once LoadHandoff() below switches to handoff mode.
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

    if (!geometryPath.empty()) LoadHandoff(geometryPath);

    FrameView();
}

Application::~Application() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void Application::LoadHandoff(const std::string& geometryPath) {
    m_Contract = Geometry::LoadHandoff(geometryPath);
    m_LatticeBuilder.emplace(m_Contract);
    m_Panels = m_LatticeBuilder->Build();
    m_BodyPanels = m_LatticeBuilder->BuildBody();
    m_DuctPanels = m_LatticeBuilder->BuildDuct();
    m_SourcePanels = m_BodyPanels;
    m_SourcePanels.insert(m_SourcePanels.end(), m_DuctPanels.begin(), m_DuctPanels.end());
    m_ReferenceArea = m_LatticeBuilder->GrossPlanformArea();
    m_ReferenceSpan = m_Contract.Span;
    m_TrimEta = m_LatticeBuilder->TrimEta();
    m_UseHandoff = true;

    // Fresh contract, fresh commands: indices from a previous handoff would
    // otherwise dangle against this one's (possibly shorter, differently
    // ordered) ControlSurfaces list.
    m_RightDeflectionDeg.assign(m_Contract.ControlSurfaces.size(), 0.0);
    m_LeftDeflectionDeg.assign(m_Contract.ControlSurfaces.size(), 0.0);

    // Bounding sphere of the whole airframe (wing + fuselage + duct), for
    // framing -- the wing's own span alone would leave the fuselage's
    // nose/tail (or a duct sitting aft of it) cut off whenever the wing is
    // placed near one end of the body.
    Solver::Vec3 lo(1e30, 1e30, 1e30), hi(-1e30, -1e30, -1e30);
    auto accumulate = [&](const Solver::Vec3& p) {
        lo = Solver::Vec3(std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z));
        hi = Solver::Vec3(std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z));
    };
    for (const auto& panel : m_Panels) { accumulate(panel.A); accumulate(panel.B); }
    for (const auto& panel : m_SourcePanels)
        for (const auto& corner : panel.Corners) accumulate(corner);

    m_SceneCenter = (lo + hi) * 0.5;
    Solver::Vec3 extent = hi - lo;
    m_SceneRadius = std::max({extent.x, extent.y, extent.z, 0.1});

    Resolve();
}

void Application::FrameView() {
    if (m_UseHandoff) {
        m_Camera.SetTarget(glm::vec3(static_cast<float>(m_SceneCenter.x), static_cast<float>(m_SceneCenter.y),
                                      static_cast<float>(m_SceneCenter.z)));
        m_Camera.SetDistance(static_cast<float>(m_SceneRadius) * CameraFrameDistanceFactor);
        return;
    }
    // Center on the mid-chord of the root section; the lattice grows aft (+x).
    m_Camera.SetTarget(glm::vec3(static_cast<float>(m_Wing.RootChord) * 0.25f, 0.0f, 0.0f));
    m_Camera.SetDistance(static_cast<float>(m_Wing.Span) * CameraFrameDistanceFactor);
}

void Application::Resolve() {
    Solver::ReferenceGeometry ref;
    double trail;

    if (m_UseHandoff) {
        // The planform itself is fixed once loaded (no Planform sliders in
        // this mode), but a control-surface deflection reshapes the wing
        // panels, so those are rebuilt from the cached LatticeBuilder every
        // Resolve() -- cheap relative to the solve, since it skips the
        // spanwise march (done once in LoadHandoff()) and only redoes the
        // per-Build() chordwise/hinge work.
        m_LatticeBuilder->ClearDeflections();
        for (std::size_t i = 0; i < m_Contract.ControlSurfaces.size(); ++i) {
            if (m_Contract.ControlSurfaces[i].Binding != Geometry::ControlSurfaceBinding::Wing) continue;
            const double rightDeg = m_RightDeflectionDeg[i];
            const double leftDeg = m_LeftDeflectionDeg[i];
            if (rightDeg != 0.0 || leftDeg != 0.0)
                m_LatticeBuilder->Deflect({i, rightDeg, leftDeg});
        }
        m_Panels = m_LatticeBuilder->Build();

        ref.Area = m_ReferenceArea;
        ref.Span = m_ReferenceSpan;
        ref.Chord = m_ReferenceArea / m_ReferenceSpan;
        trail = Solver::DefaultTrailSpanFactor * m_ReferenceSpan;

        const Solver::PanelSystem system{m_Panels, m_SourcePanels};
        const Solver::PreparedSystem prepared = Solver::Prepare(system, trail);
        m_Result = Solver::SolveWithSystem(prepared, m_Freestream, ref);
    } else {
        m_Panels = Solver::BuildWing(m_Wing);
        for (auto& panel : m_Panels) panel.Surface = "wing";

        double area = 0.0;
        for (const auto& panel : m_Panels) area += panel.Area;
        ref.Area = area;
        ref.Span = m_Wing.Span;
        ref.Chord = area / m_Wing.Span;
        trail = (m_Wing.TrailLength > 0.0) ? m_Wing.TrailLength
                                            : Solver::DefaultTrailSpanFactor * m_Wing.Span;

        m_Result = Solver::Solve(m_Panels, m_Freestream, ref, trail);
    }
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

    if (m_UseHandoff) {
        if (ImGui::CollapsingHeader("Geometry (handoff)", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Text("design_id: %.16s...", m_Contract.DesignId.c_str());
            ImGui::Text("schema: %s", m_Contract.SchemaVersion.c_str());
            ImGui::Text("wing panels: %zu", m_Panels.size());
            ImGui::Text("body panels: %zu", m_BodyPanels.size());
            ImGui::Text("duct panels: %zu", m_DuctPanels.size());
            ImGui::Text("trimmed at eta: %.3f", m_TrimEta);
        }
        if (ImGui::CollapsingHeader("Control surfaces", ImGuiTreeNodeFlags_DefaultOpen)) {
            bool any = false;
            for (std::size_t i = 0; i < m_Contract.ControlSurfaces.size(); ++i) {
                const Geometry::ControlSurface& surface = m_Contract.ControlSurfaces[i];
                // Duct-jet vanes never enter the wing lattice (their eta is a
                // duct-exit radius fraction, not a semi-span fraction -- see
                // ControlSurface.h) and PanelBuilder builds no vane panels for
                // them, so there is nothing here to deflect.
                if (surface.Binding != Geometry::ControlSurfaceBinding::Wing) continue;
                any = true;

                // Hard limits where the contract states them (schema >= 1.4.0);
                // a generous default otherwise, since MinDeg == MaxDeg == 0
                // means "not stated", not "cannot move".
                const bool hasLimits = surface.Limits.MaxDeg > surface.Limits.MinDeg;
                const double lo = hasLimits ? surface.Limits.MinDeg : -DefaultDeflectionLimitDeg;
                const double hi = hasLimits ? surface.Limits.MaxDeg : DefaultDeflectionLimitDeg;

                ImGui::PushID(static_cast<int>(i));
                ImGui::Text("%s (eta %.2f-%.2f)", surface.Name.c_str(), surface.EtaStart, surface.EtaEnd);
                // Same Right/Left addressing as PanelBuilder::ControlDeflection:
                // equal angles give a flap/elevator, opposite an aileron, and
                // nothing stops any other combination.
                m_Dirty |= SliderD("Right [deg]", &m_RightDeflectionDeg[i], lo, hi, "%.1f");
                m_Dirty |= SliderD("Left [deg]", &m_LeftDeflectionDeg[i], lo, hi, "%.1f");
                ImGui::PopID();
            }
            if (!any) ImGui::TextDisabled("(no wing-bound control surfaces)");
        }
    } else if (ImGui::CollapsingHeader("Planform", ImGuiTreeNodeFlags_DefaultOpen)) {
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
        double refSpan = m_UseHandoff ? m_ReferenceSpan : m_Wing.Span;
        double aspectRatio = (m_Result.ReferenceArea > 0.0)
                                 ? refSpan * refSpan / m_Result.ReferenceArea
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

    // Solver::ComputeDerivatives() takes a wing-only panel list (no
    // PanelSystem/body-coupling overload exists), so it would silently
    // ignore the fuselage in handoff mode -- hidden there rather than
    // offering a derivative table that quietly excludes the body.
    if (!m_UseHandoff) {
        ImGui::SeparatorText("Stability derivatives");
        if (ImGui::Button("Compute derivatives")) {
            double area = 0.0;
            for (const auto& panel : m_Panels) area += panel.Area;
            Solver::ReferenceGeometry ref;
            ref.Area = area;
            ref.Span = m_Wing.Span;
            ref.Chord = area / m_Wing.Span;
            double trail = (m_Wing.TrailLength > 0.0)
                               ? m_Wing.TrailLength
                               : Solver::DefaultTrailSpanFactor * m_Wing.Span;
            m_Derivatives = Solver::ComputeDerivatives(m_Panels, m_Freestream, ref, trail);
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
            double span = m_UseHandoff ? m_SceneRadius : m_Wing.Span;
            m_Lattice.Update(m_Panels, m_Result, m_Display, span, m_SourcePanels);
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

        if (!m_ScreenshotPath.empty() && maxFrames >= 0 && frame == maxFrames) CaptureScreenshot();
    }
}

void Application::CaptureScreenshot() const {
    int width = m_Window.FramebufferWidth();
    int height = m_Window.FramebufferHeight();
    std::vector<unsigned char> pixels(static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 3);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

    std::ofstream file(m_ScreenshotPath, std::ios::binary);
    if (!file) {
        std::fprintf(stderr, "aeolion_viewer: could not open '%s' for writing\n", m_ScreenshotPath.c_str());
        return;
    }
    file << "P6\n" << width << ' ' << height << "\n255\n";
    // glReadPixels rows run bottom-to-top; flip to the top-down row order a
    // PPM reader expects.
    for (int row = height - 1; row >= 0; --row) {
        const auto* rowStart = pixels.data() + static_cast<std::size_t>(row) * static_cast<std::size_t>(width) * 3;
        file.write(reinterpret_cast<const char*>(rowStart), static_cast<std::streamsize>(width) * 3);
    }
}

} // namespace Aeolion::Viewer
