// VLM/StationResult.h
//
// Per-spanwise-station output of a solve (one entry per panel, then sorted
// by span coordinate).
#pragma once
#include <string>

namespace Aeolion::VLM {

struct StationResult {
    double y = 0.0;
    double gamma = 0.0;
    double LiftPerSpan = 0.0; // dL/dy  [N/m]
    double cl_local = 0.0;    // local (sectional) lift coefficient
    std::string Surface;
    int PanelIndex = -1;      // index into the original (unsorted) panel/gamma array
};

} // namespace Aeolion::VLM
