// Solver/SolveResult.h
//
// Full result of a VLM solve: circulation distribution, per-station
// breakdown, integrated force/moment coefficients and their dimensional
// counterparts, and per-surface bookkeeping.
#pragma once
#include <vector>
#include <string>
#include <map>
#include "Aeolion/Solver/StationResult.h"

namespace Aeolion::Solver {

/** Full result of a VLM solve. */
struct SolveResult {
    std::vector<double> gamma;
    std::vector<StationResult> Stations;
    double CL = 0.0;
    /**
     * INDUCED drag only -- VLM is a potential-flow method and cannot predict
     * viscous/profile drag. Total CD = CDi + your own CD0 estimate.
     */
    double CDi = 0.0;
    double CY = 0.0;
    double Croll = 0.0; ///< Rolling moment coefficient about RefPoint (positive: right wingtip up).
    double Cm = 0.0;    ///< Pitching moment coefficient about RefPoint (positive: nose up).
    double Cn = 0.0;    ///< Yawing moment coefficient about RefPoint (positive: nose right, about +z).
    double L = 0.0, Di = 0.0, Y = 0.0;
    double Mx = 0.0, My = 0.0, Mz = 0.0; ///< Moments about RefPoint [N*m], same axis convention as Croll/Cm/Cn.
    double ReferenceArea = 0.0, ReferenceChord = 0.0, ReferenceSpan = 0.0;
    std::map<std::string, double> LiftBySurface;  ///< [N], per source surface tag.
    std::map<std::string, double> DragBySurface;  ///< [N]
    std::map<std::string, double> AreaBySurface;  ///< [m^2]
};

} // namespace Aeolion::Solver
