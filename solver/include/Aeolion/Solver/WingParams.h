// Solver/WingParams.h
//
// Parametric single-wing planform description consumed by BuildWing().
#pragma once

namespace Aeolion::Solver {

struct WingParams {
    double Span = 0.0;            // full span, tip-to-tip [m]
    double RootChord = 0.0;       // [m]
    double TipChord = 0.0;        // [m]
    double SweepQuarterChordDeg = 0.0; // sweep of the quarter-chord line
    double DihedralDeg = 0.0;
    double TwistTipDeg = 0.0;     // linear washout: root=0, tip=TwistTipDeg (negative = washout)
    int NPanelsSemiSpan = 0;      // panels per semi-span (total panels = 2x this)
    bool CosineSpacing = false;   // cluster panels near tips
    double TrailLength = 0.0;     // 0 => auto (50x span)
};

} // namespace Aeolion::Solver
