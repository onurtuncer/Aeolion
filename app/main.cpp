// main.cpp — example driver for the vortex lattice method in Solver.h
//
// Build:
//   g++ -std=c++23 -O2 -Iinclude -o solver_demo src/main.cpp
// Run:
//   ./solver_demo
//
// Edit the demo constants below (or wire up argv) to study your own planform.

#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Logger/Log.h"
#include <numbers>
#include <fstream>
#include <format>

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
    Logger::Log::Init();

    Solver::WingParams wing;
    wing.Span = DemoSpan;
    wing.RootChord = DemoRootChord;
    wing.TipChord = DemoTipChord;
    wing.SweepQuarterChordDeg = DemoSweepQuarterChordDeg;
    wing.DihedralDeg = DemoDihedralDeg;
    wing.TwistTipDeg = DemoTwistTipDeg;
    wing.NPanelsSemiSpan = DemoPanelsSemiSpan;
    wing.CosineSpacing = true;

    Solver::FreestreamConditions fc;
    fc.Vinf = DemoVinf;
    fc.alphaDeg = DemoAlphaDeg;
    fc.rho = DemoAirDensity;

    Solver::SolveResult res = Solver::Solve(wing, fc);

    double AR = (wing.Span * wing.Span) / res.ReferenceArea;

    AE_INFO("Reference area S   : {:.5f} m^2", res.ReferenceArea);
    AE_INFO("Aspect ratio AR    : {:.5f}", AR);
    AE_INFO("Panels (semi x2)   : {}", res.gamma.size());
    AE_INFO("-----------------------------------------");
    AE_INFO("CL                 : {:.5f}", res.CL);
    AE_INFO("CDi (induced)      : {:.5f}", res.CDi);
    AE_INFO("CY (side)          : {:.5f}", res.CY);
    AE_INFO("L/D (induced only) : {:.5f}", res.CDi > Solver::GeometryEps ? res.CL / res.CDi : 0.0);

    // Sanity check against linear thin-wing theory for an unswept,
    // untwisted flat wing: CL ~= 2*pi*alpha / (1 + 2/AR)  [alpha in rad]
    double alphaRad = Math::DegToRad(fc.alphaDeg);
    double CL_flat_theory = Math::Two * std::numbers::pi * alphaRad / (1.0 + LiftingLineArCorrection / AR);
    AE_INFO("-----------------------------------------");
    AE_INFO("Flat-wing thin-airfoil-theory CL estimate (unswept/untwisted\n"
            "reference, not directly comparable to the swept/twisted case\n"
            "above): {:.5f}", CL_flat_theory);

    // Write spanwise loading to CSV for plotting.
    std::ofstream csv("spanwise_loading.csv");
    csv << "y_m,gamma_m2_s,lift_per_span_N_m,cl_local\n";
    for (const auto& s : res.Stations) {
        csv << std::format("{:.8g},{:.8g},{:.8g},{:.8g}\n", s.y, s.gamma, s.LiftPerSpan, s.cl_local);
    }
    AE_INFO("Spanwise loading written to spanwise_loading.csv");

    // Stability derivatives about the quarter-chord of the mean aerodynamic
    // chord (a common default reference point for a wing-alone case; set
    // fc.RefPoint to your actual CG for a real aircraft).
    std::vector<Solver::Panel> panels = Solver::BuildWing(wing);
    for (auto& p : panels) p.Surface = "wing";
    Solver::ReferenceGeometry ref;
    ref.Area = res.ReferenceArea; ref.Span = wing.Span; ref.Chord = res.ReferenceArea / wing.Span;
    double trail = Solver::DefaultTrailSpanFactor * wing.Span;
    fc.RefPoint = Solver::Vec3(ref.Chord * Math::QuarterChord, 0, 0);
    Solver::StabilityDerivatives d = Solver::ComputeDerivatives(panels, fc, ref, trail);

    AE_INFO("=== Stability derivatives (wing alone, about x={:.4f}) ===", fc.RefPoint.x);
    AE_INFO("CL_alpha  = {:.4f} /rad", d.CL_alpha);
    AE_INFO("Cm_alpha  = {:.4f} /rad  ({})", d.Cm_alpha, d.Cm_alpha < 0 ? "stable" : "UNSTABLE");
    AE_INFO("CY_beta   = {:.4f} /rad", d.CY_beta);
    AE_INFO("Croll_beta= {:.4f} /rad", d.Croll_beta);

    return 0;
}
