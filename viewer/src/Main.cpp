// Main.cpp
//
// Entry point. `aeolion_viewer [--frames N]` -- the optional flag renders
// N frames and exits, used as a headless-ish smoke test in CI/builds.

#include "Core/Application.h"

#include <charconv>
#include <cstdio>
#include <string_view>

int main(int argc, char** argv) {
    int maxFrames = -1;
    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "--frames" && i + 1 < argc) {
            std::string_view value = argv[++i];
            std::from_chars(value.data(), value.data() + value.size(), maxFrames);
        }
    }

    try {
        Aeolion::Viewer::Application app;
        app.Run(maxFrames);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "aeolion_viewer: fatal: %s\n", e.what());
        return 1;
    }
    return 0;
}
