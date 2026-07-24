// Solver/StationResult.h
//
// Per-spanwise-station output of a solve (one entry per panel, then sorted
// by span coordinate).
#pragma once
#include <string>

namespace Aeolion::Solver {

/** Per-spanwise-station output of a solve. */
struct StationResult {
    double y = 0.0;
    double gamma = 0.0;
    double LiftPerSpan = 0.0; ///< dL/dy  [N/m].
    double cl_local = 0.0;    ///< Local (sectional) lift coefficient.
    std::string Surface;
    int PanelIndex = -1;      ///< Index into the original (unsorted) panel/gamma array.
};

} // namespace Aeolion::Solver
