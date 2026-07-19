// main.cpp — example driver for the vortex lattice method in VLM.h
//
// Build:
//   g++ -std=c++20 -O2 -Iinclude -o vlm_demo src/main.cpp
// Run:
//   ./vlm_demo
//
// Edit the demo constants below (or wire up argv) to study your own planform.

#include "Aeolion/VLM.h"
#include "Aeolion/Math/Constants.h"
#include <numbers>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace Aeolion;

namespace {
// --- demo wing planform ---------------------------------------------------
constexpr double DemoSpan                 = 2.0;   // m, tip-to-tip
constexpr double DemoRootChord            = 0.30;  // m
constexpr double DemoTipChord             = 0.15;  // m (taper ratio 0.5)
constexpr double DemoSweepQuarterChordDeg = 10.0;
constexpr double DemoDihedralDeg          = 3.0;
constexpr double DemoTwistTipDeg          = -3.0;  // 3 deg washout at tip
constexpr int    DemoPanelsSemiSpan       = 20;

// --- demo flight condition ------------------------------------------------
constexpr double DemoVinf       = 20.0;  // m/s
constexpr double DemoAlphaDeg   = 5.0;
constexpr double DemoAirDensity = 1.225; // kg/m^3

// thin-wing lift-slope theory: CL = 2*pi*alpha / (1 + 2/AR)
constexpr double LiftingLineArCorrection = Math::Two; // the "2" in 1 + 2/AR
}

int main() {
    VLM::WingParams wing;
    wing.Span = DemoSpan;
    wing.RootChord = DemoRootChord;
    wing.TipChord = DemoTipChord;
    wing.SweepQuarterChordDeg = DemoSweepQuarterChordDeg;
    wing.DihedralDeg = DemoDihedralDeg;
    wing.TwistTipDeg = DemoTwistTipDeg;
    wing.NPanelsSemiSpan = DemoPanelsSemiSpan;
    wing.CosineSpacing = true;

    VLM::FreestreamConditions fc;
    fc.Vinf = DemoVinf;
    fc.alphaDeg = DemoAlphaDeg;
    fc.rho = DemoAirDensity;

    VLM::SolveResult res = VLM::Solve(wing, fc);

    double AR = (wing.Span * wing.Span) / res.ReferenceArea;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Reference area S   : " << res.ReferenceArea << " m^2\n";
    std::cout << "Aspect ratio AR    : " << AR << "\n";
    std::cout << "Panels (semi x2)   : " << res.gamma.size() << "\n";
    std::cout << "-----------------------------------------\n";
    std::cout << "CL                 : " << res.CL << "\n";
    std::cout << "CDi (induced)      : " << res.CDi << "\n";
    std::cout << "CY (side)          : " << res.CY << "\n";
    std::cout << "L/D (induced only) : " << (res.CDi > VLM::GeometryEps ? res.CL / res.CDi : 0.0) << "\n";

    // Sanity check against linear thin-wing theory for an unswept,
    // untwisted flat wing: CL ~= 2*pi*alpha / (1 + 2/AR)  [alpha in rad]
    double alphaRad = Math::DegToRad(fc.alphaDeg);
    double CL_flat_theory = Math::Two * std::numbers::pi * alphaRad / (1.0 + LiftingLineArCorrection / AR);
    std::cout << "-----------------------------------------\n";
    std::cout << "Flat-wing thin-airfoil-theory CL estimate (unswept/untwisted\n"
                 "reference, not directly comparable to the swept/twisted case\n"
                 "above): " << CL_flat_theory << "\n";

    // Write spanwise loading to CSV for plotting.
    std::ofstream csv("spanwise_loading.csv");
    csv << "y_m,gamma_m2_s,lift_per_span_N_m,cl_local\n";
    csv << std::setprecision(8);
    for (const auto& s : res.Stations) {
        csv << s.y << "," << s.gamma << "," << s.LiftPerSpan << "," << s.cl_local << "\n";
    }
    csv.close();
    std::cout << "\nSpanwise loading written to spanwise_loading.csv\n";

    // Stability derivatives about the quarter-chord of the mean aerodynamic
    // chord (a common default reference point for a wing-alone case; set
    // fc.RefPoint to your actual CG for a real aircraft).
    std::vector<VLM::Panel> panels = VLM::BuildWing(wing);
    for (auto& p : panels) p.Surface = "wing";
    VLM::ReferenceGeometry ref;
    ref.Area = res.ReferenceArea; ref.Span = wing.Span; ref.Chord = res.ReferenceArea / wing.Span;
    double trail = VLM::DefaultTrailSpanFactor * wing.Span;
    fc.RefPoint = VLM::Vec3(ref.Chord * Math::QuarterChord, 0, 0);
    VLM::StabilityDerivatives d = VLM::ComputeDerivatives(panels, fc, ref, trail);

    std::cout << std::setprecision(4);
    std::cout << "\n=== Stability derivatives (wing alone, about x=" << fc.RefPoint.x << ") ===\n";
    std::cout << "CL_alpha  = " << d.CL_alpha << " /rad\n";
    std::cout << "Cm_alpha  = " << d.Cm_alpha << " /rad  (" << (d.Cm_alpha < 0 ? "stable" : "UNSTABLE") << ")\n";
    std::cout << "CY_beta   = " << d.CY_beta << " /rad\n";
    std::cout << "Croll_beta= " << d.Croll_beta << " /rad\n";

    return 0;
}
