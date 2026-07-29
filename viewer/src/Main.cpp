// Main.cpp
//
// Entry point.
//   aeolion_viewer [--frames N] [--geometry FILE] [--screenshot FILE]
//                  [--screen airframe|propeller]
// --frames renders N frames and exits, used as a headless-ish smoke test in
// CI/builds. --geometry loads an aeolion_geometry.json handoff (wing +
// fuselage, and the propeller geometry if the contract carries a
// propulsion_bemt block) instead of the parametric single-wing demo.
// --screenshot captures the last rendered frame to a binary PPM file; meant
// to be paired with a small --frames so the window closes itself once the
// capture is written (a couple of frames lets ImGui settle its first-use
// layout). --screen picks which view opens first: the airframe lattice
// (VLM) or the propeller geometry.

#include "Core/Application.h"

#include "Aeolion/Logger/Log.h"

#include <charconv>
#include <string>
#include <string_view>

int main(int argc, char** argv) {
    Aeolion::Logger::Log::Init();

    int maxFrames = -1;
    std::string geometryPath;
    std::string screenshotPath;
    double orbitYawDeg = 0.0, orbitPitchDeg = 0.0;
    Aeolion::Viewer::Screen screen = Aeolion::Viewer::Screen::Airframe;
    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "--frames" && i + 1 < argc) {
            std::string_view value = argv[++i];
            std::from_chars(value.data(), value.data() + value.size(), maxFrames);
        } else if (arg == "--geometry" && i + 1 < argc) {
            geometryPath = argv[++i];
        } else if (arg == "--screenshot" && i + 1 < argc) {
            screenshotPath = argv[++i];
        } else if (arg == "--screen" && i + 1 < argc) {
            std::string_view value = argv[++i];
            if (value == "propeller") screen = Aeolion::Viewer::Screen::Propeller;
            else if (value == "airframe") screen = Aeolion::Viewer::Screen::Airframe;
            else AE_WARN("aeolion_viewer: unknown --screen '{}', staying on airframe", value);
        } else if (arg == "--orbit" && i + 2 < argc) {
            // Yaw and pitch offsets [deg] from the framed pose, for
            // scripted screenshots from chosen viewpoints.
            std::string_view yawText = argv[++i];
            std::string_view pitchText = argv[++i];
            std::from_chars(yawText.data(), yawText.data() + yawText.size(), orbitYawDeg);
            std::from_chars(pitchText.data(), pitchText.data() + pitchText.size(), orbitPitchDeg);
        }
    }

    try {
        Aeolion::Viewer::Application app(geometryPath, screenshotPath, screen);
        if (orbitYawDeg != 0.0 || orbitPitchDeg != 0.0) app.OrbitBy(orbitYawDeg, orbitPitchDeg);
        app.Run(maxFrames);
    } catch (const std::exception& e) {
        AE_CRITICAL("aeolion_viewer: fatal: {}", e.what());
        return 1;
    }
    return 0;
}
