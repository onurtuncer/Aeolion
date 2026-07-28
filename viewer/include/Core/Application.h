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
#include "Visualization/PropellerRenderer.h"

#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/Solver.h"

#include <optional>
#include <string>
#include <vector>

namespace Aeolion::Viewer {

// Which subject the viewer is showing. Airframe = the wing/fuselage lattice
// (LatticeRenderer + VLM); Propeller = the isolated propeller geometry
// (PropellerRenderer -- no aerodynamic solve behind it since BEMT moved to
// its own project; a future calculation method plugs in behind this same
// screen). The two carry independent geometry, camera framing, and UI
// panels; only the window, GL context, and orbit camera are shared.
enum class Screen {
    Airframe = 0,
    Propeller,
};

class Application {
public:
    // geometryPath: if non-empty, loads an aeolion_geometry.json handoff
    // (wing lattice + fuselage source panels via PanelBuilder::LatticeBuilder)
    // instead of the parametric single-wing demo. screenshotPath: if
    // non-empty, captures the framebuffer to a binary PPM file on the last
    // rendered frame -- meant to be paired with a small, finite maxFrames.
    explicit Application(const std::string& geometryPath = "", const std::string& screenshotPath = "",
                         Screen initialScreen = Screen::Airframe);
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
    void DrawScreenSelector();
    void LoadHandoff(const std::string& geometryPath);
    void CaptureScreenshot() const;

    // --- propeller screen ---------------------------------------------------
    void BuildDefaultPropeller();               // fallback prop when no contract states one
    void AdoptContractPropeller();              // use the loaded handoff's propulsion_bemt block
    void ResolveProp();                         // mesh the blades + rotating-frame VLM solve + remesh
    void FrameProp();                           // frame the camera on the (small) propeller
    void DrawPropellerUI();

    Window m_Window;
    OrbitCamera m_Camera;
    Screen m_Screen = Screen::Airframe;
    LatticeRenderer m_Lattice;
    PropellerRenderer m_Propeller;

    Solver::WingParams m_Wing;
    Solver::FreestreamConditions m_Freestream;
    LatticeDisplayOptions m_Display;

    // Handoff-geometry mode: the geometry comes from the contract (it does
    // not respond to the Planform sliders, which drive m_Wing instead)
    // rather than being rebuilt from scratch. The builder itself is kept
    // (not just its Build() output) so a control-surface deflection can
    // re-run Build() cheaply without re-deriving the planform march.
    bool m_UseHandoff = false;
    Geometry::HandoffContract m_Contract;
    std::optional<PanelBuilder::LatticeBuilder> m_LatticeBuilder;
    std::vector<Lattice::SourcePanel> m_BodyPanels;
    std::vector<Lattice::SourcePanel> m_DuctPanels;
    // Fuselage + duct combined: what the solve and the renderer actually
    // consume. Cached once alongside m_BodyPanels/m_DuctPanels (kept
    // separate for the per-surface panel-count readout in DrawUI()) rather
    // than reassembled every Resolve(), same as the rest of handoff mode's
    // "build once, only the freestream re-solves per frame" caching. Neither
    // fuselage nor duct move with a control deflection, so both stay fixed
    // for the contract's lifetime, unlike m_Panels.
    std::vector<Lattice::SourcePanel> m_SourcePanels;
    double m_ReferenceArea = 0.0;
    double m_ReferenceSpan = 0.0;
    double m_TrimEta = 0.0;
    Solver::Vec3 m_SceneCenter{0.0, 0.0, 0.0};
    double m_SceneRadius = 1.0;

    // Commanded deflection per HandoffContract::ControlSurfaces index, in
    // degrees; only meaningful where Binding == Wing (see DrawUI() and
    // Resolve()). Sized to m_Contract.ControlSurfaces on every LoadHandoff(),
    // and indexed the same way LatticeBuilder::ControlDeflection is, so a
    // slider maps directly onto Deflect() without a lookup.
    std::vector<double> m_RightDeflectionDeg;
    std::vector<double> m_LeftDeflectionDeg;

    std::string m_ScreenshotPath;

    std::vector<Solver::Panel> m_Panels;
    Solver::SolveResult m_Result;
    std::optional<Solver::StabilityDerivatives> m_Derivatives;

    bool m_Dirty = true;         // geometry/condition changed -> re-solve + remesh
    bool m_DisplayDirty = false; // display options changed -> remesh only

    // --- propeller screen state ---------------------------------------------
    // The propeller is independent of the airframe: its own geometry (from
    // the handoff's propulsion_bemt block via Geometry::ToPropeller, or a
    // built-in default), its own operating point, and its own solve. The
    // blades are meshed by PanelBuilder::BuildPropellerLattice and run
    // through the SAME Solver::Solve as the airframe, spinning via the
    // body-rate term (fc.p = Omega about +x) -- Aeolion's panel-native
    // propeller method. m_PropDirty covers mesh + solve + remesh together;
    // the lattice is small enough not to split them.
    Geometry::Propeller m_Prop;
    std::vector<Solver::Panel> m_PropPanels;
    Solver::SolveResult m_PropSolve;
    LatticeRenderer m_PropLattice;              // the solved blade lattice, field-colored
    LatticeDisplayOptions m_PropLatticeDisplay; // grid/axes off; see constructor
    PropellerDisplayOptions m_PropDisplay;      // hub/disk/axis decoration (geometry ribbon off)
    bool m_PropFromContract = false; // true if m_Prop came from a loaded contract, not the default
    double m_PropRpm = 6000.0;       // shaft speed; Omega = rpm * 2*pi/60 about +x
    double m_PropSpeed = 0.0;        // axial inflow [m/s]; 0 = hover (solved at a tiny floor speed)
    double m_PropDensity = 1.225;    // [kg/m^3], sea-level ISA
    bool m_PropDirty = true;         // geometry/operating point/display changed -> remesh + re-solve

    double m_LastMouseX = 0.0, m_LastMouseY = 0.0;
};

} // namespace Aeolion::Viewer
