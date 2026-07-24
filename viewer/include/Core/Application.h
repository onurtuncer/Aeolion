// Core/Application.h
//
// The viewer application: owns the window, ImGui lifetime, orbit camera,
// the solver state (WingParams + FreestreamConditions, or a loaded
// Geometry::HandoffContract's wing + fuselage lattice), and the lattice
// renderer. Each frame: consume input, re-solve if a parameter changed,
// draw the scene, draw the UI.

#pragma once

#include "Core/Window.h"
#include "Renderer/OrbitCamera.h"
#include "Visualization/LatticeRenderer.h"

#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/Solver/Solver.h"

#include <optional>
#include <string>
#include <vector>

namespace Aeolion::Viewer {

class Application {
public:
    // geometryPath: if non-empty, loads an aeolion_geometry.json handoff
    // (wing lattice + fuselage source panels via PanelBuilder::LatticeBuilder)
    // instead of the parametric single-wing demo. screenshotPath: if
    // non-empty, captures the framebuffer to a binary PPM file on the last
    // rendered frame -- meant to be paired with a small, finite maxFrames.
    explicit Application(const std::string& geometryPath = "", const std::string& screenshotPath = "");
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
    void LoadHandoff(const std::string& geometryPath);
    void CaptureScreenshot() const;

    Window m_Window;
    OrbitCamera m_Camera;
    LatticeRenderer m_Lattice;

    Solver::WingParams m_Wing;
    Solver::FreestreamConditions m_Freestream;
    LatticeDisplayOptions m_Display;

    // Handoff-geometry mode: the lattice is built once from the contract
    // (it does not respond to the Planform sliders, which drive m_Wing
    // instead) and cached here; only the freestream re-solves per frame.
    bool m_UseHandoff = false;
    Geometry::HandoffContract m_Contract;
    std::vector<Lattice::SourcePanel> m_BodyPanels;
    std::vector<Lattice::SourcePanel> m_DuctPanels;
    // Fuselage + duct combined: what the solve and the renderer actually
    // consume. Cached once alongside m_BodyPanels/m_DuctPanels (kept
    // separate for the per-surface panel-count readout in DrawUI()) rather
    // than reassembled every Resolve(), same as the rest of handoff mode's
    // "build once, only the freestream re-solves per frame" caching.
    std::vector<Lattice::SourcePanel> m_SourcePanels;
    double m_ReferenceArea = 0.0;
    double m_ReferenceSpan = 0.0;
    double m_TrimEta = 0.0;
    Solver::Vec3 m_SceneCenter{0.0, 0.0, 0.0};
    double m_SceneRadius = 1.0;

    std::string m_ScreenshotPath;

    std::vector<Solver::Panel> m_Panels;
    Solver::SolveResult m_Result;
    std::optional<Solver::StabilityDerivatives> m_Derivatives;

    bool m_Dirty = true;         // geometry/condition changed -> re-solve + remesh
    bool m_DisplayDirty = false; // display options changed -> remesh only

    double m_LastMouseX = 0.0, m_LastMouseY = 0.0;
};

} // namespace Aeolion::Viewer
