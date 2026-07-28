// Core/Application.cpp

#include "Core/Application.h"
#include "Visualization/ColorMap.h"

#include "Aeolion/Logger/Log.h"
#include "Aeolion/PanelBuilder/PanelBuilder.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <algorithm>
#include <cfloat>
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
constexpr float PropFrameDistanceFactor = 3.2f;   // orbit distance = this * tip radius on the propeller screen

// Built-in default propeller (used until a contract states propulsion_bemt):
// a small 2-blade fixed-pitch prop with a conventional taper + washout.
constexpr int DefaultPropBlades = 2;
constexpr double DefaultPropRadius = 0.15;   // [m]
constexpr double DefaultPropHubRadius = 0.02; // [m]
constexpr int DefaultPropStations = 12;
constexpr double DefaultPropGeometricPitch = 0.12; // [m], sets the twist(r) = atan(P / 2*pi*r) washout
constexpr double DefaultPropRootChord = 0.032;     // [m]
constexpr double DefaultPropTaper = 0.55;          // chord = root * (1 - taper * r/R)

// Default shroud proportions when the contract states no duct: a snug
// tip clearance, a modest ring thickness, and a chord that wraps the
// rotor plane -- enough to show the ducted-fan interaction without
// pretending to be a designed duct.
constexpr double DefaultShroudTipClearance = 1.04; // inner radius = this * tip radius
constexpr double DefaultShroudOuterFactor = 1.20;  // outer radius = this * tip radius
constexpr double DefaultShroudChordFactor = 0.5;   // chord = this * tip radius

// Hover means zero axial inflow, but the solver's wind-axis force
// projections need a nonzero translational freestream to be defined
// (drag/lift directions normalize Vinf) -- solve at this floor instead.
// The dimensional forces are dominated by the Omega*r rotation term, so
// the floor's contribution is negligible.
constexpr double MinPropAirspeed = 0.01; // [m/s]

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

Application::Application(const std::string& geometryPath, const std::string& screenshotPath,
                         Screen initialScreen)
    : m_Window("Aeolion Viewer", DefaultWindowWidth, DefaultWindowHeight),
      m_Screen(initialScreen),
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

    BuildDefaultPropeller();
    // The propeller screen's lattice is the subject itself -- the airframe
    // screen's ground grid and body axes would just clutter a 0.1 m prop.
    m_PropLatticeDisplay.ShowGrid = false;
    m_PropLatticeDisplay.ShowAxes = false;
    // The flat geometry ribbon coincides with the solved lattice surface, so
    // it stays off by default (z-fighting); hub/disk/axis decoration stays.
    m_PropDisplay.ShowBlades = false;

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(m_Window.Handle(), true);
    ImGui_ImplOpenGL3_Init("#version 330");

    if (!geometryPath.empty()) {
        LoadHandoff(geometryPath);
        // Prefer the contract's own propeller over the built-in default when
        // the handoff carries a propulsion_bemt block.
        if (!m_Contract.Propulsion.BladeStations.empty()) AdoptContractPropeller();
    }

    if (m_Screen == Screen::Propeller) FrameProp();
    else FrameView();
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

void Application::BuildDefaultPropeller() {
    m_Prop = Geometry::Propeller{};
    m_Prop.BladeCount = DefaultPropBlades;
    m_Prop.Radius = DefaultPropRadius;
    m_Prop.HubRadius = DefaultPropHubRadius;
    m_Prop.Stations.clear();
    m_Prop.Stations.reserve(DefaultPropStations);
    for (int i = 0; i < DefaultPropStations; ++i) {
        const double frac = static_cast<double>(i) / (DefaultPropStations - 1); // 0..1
        const double r = m_Prop.HubRadius + frac * (m_Prop.Radius - m_Prop.HubRadius);
        const double rOverR = r / m_Prop.Radius;
        const double chord = DefaultPropRootChord * (1.0 - DefaultPropTaper * rOverR);
        // Ideal fixed-pitch twist: the chordline of every element subtends the
        // same helix, so twist(r) = atan(pitch / (2*pi*r)) -- steep at the
        // root, shallow at the tip.
        const double twistDeg = Math::RadToDeg(std::atan2(DefaultPropGeometricPitch, 2.0 * std::numbers::pi * r));
        m_Prop.Stations.push_back({r, chord, twistDeg});
    }
    m_PropFromContract = false;
    m_PropDirty = true;
}

void Application::AdoptContractPropeller() {
    try {
        m_Prop = Geometry::ToPropeller(m_Contract.Propulsion);
        m_PropFromContract = true;
        m_PropDirty = true;
    } catch (const std::exception& e) {
        // Older contracts (schema < 1.4.0) state no blade count, so
        // ToPropeller throws; keep the built-in default rather than fail.
        AE_WARN("propeller: contract propulsion_bemt unusable ({}); keeping the default propeller", e.what());
    }
}

void Application::ResolveProp() {
    const double omega = m_PropRpm * 2.0 * std::numbers::pi / Math::SecondsPerMinute;
    m_PropPanels = PanelBuilder::BuildPropellerLattice(m_Prop, m_PropSpeed, omega);

    // The duct shroud: the contract's duct ring when it states one, a
    // default-proportioned shroud otherwise -- centered on the rotor plane
    // (BuildPropellerDuct's own convention).
    m_PropDuctPanels.clear();
    if (m_PropDuctEnabled) {
        if (m_UseHandoff && m_Contract.Duct.IsStated) {
            m_PropDuctPanels = PanelBuilder::BuildPropellerDuct(
                m_Contract.Duct.InnerDiameter * 0.5, m_Contract.Duct.OuterDiameter * 0.5,
                m_Contract.Duct.Chord);
        } else {
            m_PropDuctPanels = PanelBuilder::BuildPropellerDuct(
                DefaultShroudTipClearance * m_Prop.Radius, DefaultShroudOuterFactor * m_Prop.Radius,
                DefaultShroudChordFactor * m_Prop.Radius);
        }
    }

    // The rotation IS the flight condition: Omega enters as the solver's
    // roll rate about +x, whose kinematic-velocity term gives every blade
    // panel its true Omega x r tangential onset flow about the hub
    // (RefPoint at the origin, where the lattice is built). Thrust and
    // shaft torque then fall out of the body-axis force/moment sums:
    // thrust = -Di (upstream, -x), torque = -Mx (opposing +Omega).
    Solver::FreestreamConditions fc;
    fc.Vinf = std::max(m_PropSpeed, MinPropAirspeed);
    fc.alphaDeg = 0.0;
    fc.betaDeg = 0.0;
    fc.rho = m_PropDensity;
    fc.p = omega;
    fc.RefPoint = Solver::Vec3(0.0, 0.0, 0.0);

    // Reference quantities only normalize coefficients, which are not
    // meaningful at the hover floor speed anyway -- the readouts below are
    // dimensional. Still, state something sensible.
    Solver::ReferenceGeometry ref;
    double area = 0.0;
    for (const auto& panel : m_PropPanels) area += panel.PlanformArea;
    ref.Area = std::max(area, 1e-9);
    ref.Span = 2.0 * m_Prop.Radius;
    ref.Chord = ref.Area / std::max(ref.Span, 1e-9);

    const double trail = Solver::DefaultTrailSpanFactor * 2.0 * m_Prop.Radius;
    if (m_PropUseViscous) {
        // Level-2 sectional lift feedback: the lattice supplies alpha_eff,
        // the analytic viscous section model supplies cl/cd, and the
        // circulation relaxes until they agree. Base has SolveResult's
        // shape, so everything downstream reads it unchanged.
        m_PropStrips = PanelBuilder::BuildPropellerStrips(m_Prop);
        const Solver::SectionModel model =
            (m_PropSectionModel == 1) ? PanelBuilder::MakePropellerSectionModel(m_Prop)
                                      : Solver::SectionModel(Solver::AnalyticSectionModel{});
        m_PropCoupled = Solver::SolveViscousCoupled(m_PropPanels, m_PropStrips, fc, ref, trail, model,
                                                    {}, m_PropDuctPanels);
        m_PropSolve = m_PropCoupled.Base;
    } else {
        if (m_PropDuctPanels.empty()) {
            m_PropSolve = Solver::Solve(m_PropPanels, fc, ref, trail);
        } else {
            const Solver::PanelSystem system{m_PropPanels, m_PropDuctPanels};
            m_PropSolve = Solver::SolveWithSystem(Solver::Prepare(system, trail), fc, ref);
        }
        m_PropCoupled = {};
    }

    m_PropLattice.Update(m_PropPanels, m_PropSolve, m_PropLatticeDisplay, 2.0 * m_Prop.Radius,
                         m_PropDuctPanels);
    m_Propeller.Update(m_Prop, m_PropDisplay);
}

void Application::FrameProp() {
    m_Camera.SetTarget(glm::vec3(0.0f, 0.0f, 0.0f));
    m_Camera.SetDistance(static_cast<float>(m_Prop.Radius) * PropFrameDistanceFactor);
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

void Application::DrawScreenSelector() {
    // A menu-bar toggle between the airframe and propeller screens. Switching
    // reframes the camera, since the two subjects differ in scale by ~30x.
    if (ImGui::BeginMainMenuBar()) {
        ImGui::TextUnformatted("Screen:");
        if (ImGui::MenuItem("Airframe", nullptr, m_Screen == Screen::Airframe) &&
            m_Screen != Screen::Airframe) {
            m_Screen = Screen::Airframe;
            FrameView();
        }
        if (ImGui::MenuItem("Propeller", nullptr, m_Screen == Screen::Propeller) &&
            m_Screen != Screen::Propeller) {
            m_Screen = Screen::Propeller;
            FrameProp();
        }
        ImGui::EndMainMenuBar();
    }
}

void Application::DrawUI() {
    DrawScreenSelector();
    if (m_Screen == Screen::Propeller) {
        DrawPropellerUI();
        return;
    }

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

void Application::DrawPropellerUI() {
    // --- controls -----------------------------------------------------------
    ImGui::SetNextWindowPos(ImVec2(10, 28), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(360, 460), ImGuiCond_FirstUseEver);
    ImGui::Begin("Propeller");

    if (ImGui::CollapsingHeader("Geometry", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Text("source: %s", m_PropFromContract ? "handoff contract" : "built-in default");
        ImGui::Text("radius: %.4f m   hub: %.4f m", m_Prop.Radius, m_Prop.HubRadius);
        ImGui::Text("stations: %zu", m_Prop.Stations.size());
        int blades = m_Prop.BladeCount;
        if (ImGui::SliderInt("Blade count", &blades, 2, 6)) {
            m_Prop.BladeCount = blades;
            m_PropDirty = true;
        }
    }

    if (ImGui::CollapsingHeader("Operating point", ImGuiTreeNodeFlags_DefaultOpen)) {
        m_PropDirty |= SliderD("RPM", &m_PropRpm, 500.0, 15000.0, "%.0f");
        m_PropDirty |= SliderD("Airspeed [m/s]", &m_PropSpeed, 0.0, 40.0, "%.1f");
        m_PropDirty |= SliderD("rho [kg/m^3]", &m_PropDensity, 0.4, 1.4);
        m_PropDirty |= ImGui::Checkbox("Viscous coupling", &m_PropUseViscous);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Sectional lift feedback: a 2-D viscous section model\n"
                              "iterated against the lattice's induced flow.");
        if (m_PropUseViscous) {
            const char* modelNames[] = {"Analytic polar", "Boundary layer (transpiration)"};
            m_PropDirty |= ImGui::Combo("Section model", &m_PropSectionModel, modelNames,
                                        IM_ARRAYSIZE(modelNames));
        }
        m_PropDirty |= ImGui::Checkbox("Duct", &m_PropDuctEnabled);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Shroud the propeller with the contract's duct ring (or a\n"
                              "default-proportioned one) and solve them coupled: the blades\n"
                              "see the duct's induced flow, the duct sees the propwash.");
    }

    if (ImGui::CollapsingHeader("Lattice", ImGuiTreeNodeFlags_DefaultOpen)) {
        // One Weissinger row per radial strip -- see BuildPropellerLattice's
        // header for why chordwise stacking is off the table for now.
        ImGui::Text("panels: %zu (1 chordwise row per strip)", m_PropPanels.size());
    }

    if (ImGui::CollapsingHeader("Display", ImGuiTreeNodeFlags_DefaultOpen)) {
        const char* fieldNames[] = {"Circulation gamma", "Sectional cl", "Lift per span"};
        int field = static_cast<int>(m_PropLatticeDisplay.Field);
        if (ImGui::Combo("Field", &field, fieldNames, IM_ARRAYSIZE(fieldNames))) {
            m_PropLatticeDisplay.Field = static_cast<FieldMode>(field);
            m_PropDirty = true;
        }
        m_PropDirty |= ImGui::Checkbox("Wireframe", &m_PropLatticeDisplay.ShowWireframe);
        m_PropDirty |= ImGui::Checkbox("Normals", &m_PropLatticeDisplay.ShowNormals);
        m_PropDirty |= ImGui::Checkbox("Hub", &m_PropDisplay.ShowHub);
        m_PropDirty |= ImGui::Checkbox("Disk & axis", &m_PropDisplay.ShowDisk);
        m_PropDirty |= ImGui::Checkbox("Geometry ribbon", &m_PropDisplay.ShowBlades);
        if (ImGui::Button("Frame view")) FrameProp();
    }
    ImGui::End();

    // --- results ------------------------------------------------------------
    ImGui::SetNextWindowPos(ImVec2(static_cast<float>(m_Window.FramebufferWidth()) - 370, 28),
                            ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(360, 560), ImGuiCond_FirstUseEver);
    ImGui::Begin("Propeller results");

    const double diskArea = std::numbers::pi * m_Prop.Radius * m_Prop.Radius;
    // Blade planform area by the trapezoid rule over the chord distribution;
    // solidity is the classic sigma = N * bladeArea / diskArea.
    double bladeArea = 0.0;
    for (std::size_t i = 0; i + 1 < m_Prop.Stations.size(); ++i) {
        const double dr = m_Prop.Stations[i + 1].r - m_Prop.Stations[i].r;
        bladeArea += 0.5 * (m_Prop.Stations[i].Chord + m_Prop.Stations[i + 1].Chord) * dr;
    }

    // The rotating-lattice solve's body-axis sums, restated as propeller
    // quantities (see BuildPropellerLattice's conventions): thrust points
    // upstream (-x), shaft torque opposes +Omega. With viscous coupling the
    // torque is real (induced + profile); the bare inviscid lattice's is
    // induced only.
    const double omega = m_PropRpm * 2.0 * std::numbers::pi / Math::SecondsPerMinute;
    const double thrust = -m_PropSolve.Di;
    const double torque = -m_PropSolve.Mx;
    const double power = torque * omega;
    if (ImGui::BeginTable("prop_solve", 2, ImGuiTableFlags_SizingStretchProp)) {
        ReadoutRow("Thrust [N]", thrust, "%.3f");
        if (m_PropUseViscous) {
            ReadoutRow("Torque [N.m]", torque, "%.4f");
            ReadoutRow("  induced", -m_PropCoupled.InducedMoment.x, "%.4f");
            ReadoutRow("  profile", -m_PropCoupled.ProfileMoment.x, "%.4f");
            ReadoutRow("Power [W]", power, "%.1f");
        } else {
            ReadoutRow("Torque [N.m] (induced)", torque, "%.4f");
            ReadoutRow("Power [W] (induced)", power, "%.1f");
        }
        if (!m_PropDuctPanels.empty()) {
            // The duct's share of the axial force, from the per-surface
            // pressure bookkeeping (drag = force along +x; thrust is -x).
            double ductDrag = 0.0;
            for (const auto& [surface, drag] : m_PropSolve.DragBySurface)
                if (surface.rfind("duct", 0) == 0) ductDrag += drag;
            ReadoutRow("Duct thrust [N]", -ductDrag, "%.3f");
        }
        if (diskArea > 0.0) ReadoutRow("Disk loading [Pa]", thrust / diskArea, "%.1f");
        if (power > 1e-9 && m_PropSpeed > 0.0)
            ReadoutRow("Propulsive eff.", thrust * m_PropSpeed / power, "%.3f");
        ImGui::EndTable();
    }
    if (m_PropUseViscous) {
        ImGui::Text("coupling: %d iteration(s), residual %.2e", m_PropCoupled.Iterations,
                    m_PropCoupled.MaxResidual);
        if (!m_PropCoupled.Converged)
            ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.2f, 1.0f), "coupling did not converge");
    }

    ImGui::SeparatorText("Geometry");
    if (ImGui::BeginTable("prop_geom", 2, ImGuiTableFlags_SizingStretchProp)) {
        ReadoutRow("Blades", static_cast<double>(m_Prop.BladeCount), "%.0f");
        ReadoutRow("Tip radius [m]", m_Prop.Radius, "%.4f");
        ReadoutRow("Hub radius [m]", m_Prop.HubRadius, "%.4f");
        ReadoutRow("Disk area [m^2]", diskArea, "%.4f");
        if (diskArea > 0.0)
            ReadoutRow("Solidity", m_Prop.BladeCount * bladeArea / diskArea, "%.4f");
        ImGui::EndTable();
    }

    ImGui::SeparatorText("Field range");
    ColorBar(m_PropLattice.FieldMin(), m_PropLattice.FieldMax());

    ImGui::SeparatorText("Radial distribution (hub -> tip)");
    if (!m_Prop.Stations.empty()) {
        const std::size_t n = m_Prop.Stations.size();
        std::vector<float> chord(n), twist(n);
        for (std::size_t i = 0; i < n; ++i) {
            chord[i] = static_cast<float>(m_Prop.Stations[i].Chord);
            twist[i] = static_cast<float>(m_Prop.Stations[i].TwistDeg);
        }
        ImGui::PlotLines("##chord", chord.data(), static_cast<int>(n), 0, "chord [m]",
                         FLT_MAX, FLT_MAX, ImVec2(-1, 80));
        ImGui::PlotLines("##twist", twist.data(), static_cast<int>(n), 0, "twist [deg]",
                         FLT_MAX, FLT_MAX, ImVec2(-1, 80));
    }
    // One blade's converged section states -- the coupling's own output.
    if (m_PropUseViscous && m_Prop.Stations.size() > 1) {
        const std::size_t per = m_Prop.Stations.size() - 1;
        if (m_PropCoupled.Strips.size() >= per) {
            std::vector<float> alphaEff(per), cl(per);
            for (std::size_t i = 0; i < per; ++i) {
                alphaEff[i] = static_cast<float>(m_PropCoupled.Strips[i].alphaEffDeg);
                cl[i] = static_cast<float>(m_PropCoupled.Strips[i].cl);
            }
            ImGui::PlotLines("##alphaeff", alphaEff.data(), static_cast<int>(per), 0,
                             "alpha_eff [deg]", FLT_MAX, FLT_MAX, ImVec2(-1, 80));
            ImGui::PlotLines("##clvisc", cl.data(), static_cast<int>(per), 0, "section cl",
                             FLT_MAX, FLT_MAX, ImVec2(-1, 80));
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
        if (m_PropDirty) {
            ResolveProp();
            m_PropDirty = false;
        }

        int width = m_Window.FramebufferWidth();
        int height = m_Window.FramebufferHeight();
        m_Camera.SetViewport(width, height);

        glViewport(0, 0, width, height);
        glClearColor(0.09f, 0.10f, 0.13f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        if (m_Screen == Screen::Propeller) {
            m_PropLattice.Draw(m_Camera.ViewProjection());
            m_Propeller.Draw(m_Camera.ViewProjection()); // hub, disk circle, axis (and optional ribbon)
        } else {
            m_Lattice.Draw(m_Camera.ViewProjection());
        }

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
        AE_ERROR("aeolion_viewer: could not open '{}' for writing", m_ScreenshotPath);
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
