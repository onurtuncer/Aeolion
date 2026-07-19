// Core/Application.h
//
// The viewer application: owns the window, ImGui lifetime, orbit camera,
// the solver state (WingParams + FreestreamConditions), and the lattice
// renderer. Each frame: consume input, re-solve if a parameter changed,
// draw the scene, draw the UI.

#pragma once

#include "Core/Window.h"
#include "Renderer/OrbitCamera.h"
#include "Visualization/LatticeRenderer.h"

#include "Aeolion/VLM/VLM.h"

#include <optional>
#include <vector>

namespace Aeolion::Viewer {

class Application {
public:
    Application();
    ~Application();

    Application(const Application&) = delete;
    Application& operator=(const Application&) = delete;

    // maxFrames >= 0 renders that many frames then exits (smoke-test mode).
    void Run(int maxFrames = -1);

private:
    void HandleCameraInput();
    void Resolve();
    void FrameView();
    void DrawUI();

    Window m_Window;
    OrbitCamera m_Camera;
    LatticeRenderer m_Lattice;

    VLM::WingParams m_Wing;
    VLM::FreestreamConditions m_Freestream;
    LatticeDisplayOptions m_Display;

    std::vector<VLM::Panel> m_Panels;
    VLM::SolveResult m_Result;
    std::optional<VLM::StabilityDerivatives> m_Derivatives;

    bool m_Dirty = true;         // geometry/condition changed -> re-solve + remesh
    bool m_DisplayDirty = false; // display options changed -> remesh only

    double m_LastMouseX = 0.0, m_LastMouseY = 0.0;
};

} // namespace Aeolion::Viewer
