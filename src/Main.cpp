// Main.cpp — example driver for the vortex lattice method in Vlm.h
//
// Build:
//   g++ -std=c++17 -O2 -Iinclude -o vlm_demo src/Main.cpp
// Run:
//   ./vlm_demo
//
// Edit the WingParams / FreestreamConditions below (or wire up argv) to
// study your own planform.

#include "Aeolion/Vlm.h"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace Aeolion;

int main() {
    Vlm::WingParams wing;
    wing.Span = 2.0;
    wing.RootChord = 0.30;
    wing.TipChord = 0.15;          // taper ratio 0.5
    wing.SweepQuarterChordDeg = 10.0;
    wing.DihedralDeg = 3.0;
    wing.TwistTipDeg = -3.0;       // 3 deg washout at tip
    wing.NPanelsSemiSpan = 20;
    wing.CosineSpacing = true;

    Vlm::FreestreamConditions fc;
    fc.Vinf = 20.0;   // m/s
    fc.alphaDeg = 5.0;
    fc.rho = 1.225;

    Vlm::SolveResult res = Vlm::Solve(wing, fc);

    double AR = (wing.Span * wing.Span) / res.ReferenceArea;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Reference area S   : " << res.ReferenceArea << " m^2\n";
    std::cout << "Aspect ratio AR    : " << AR << "\n";
    std::cout << "Panels (semi x2)   : " << res.gamma.size() << "\n";
    std::cout << "-----------------------------------------\n";
    std::cout << "CL                 : " << res.CL << "\n";
    std::cout << "CDi (induced)      : " << res.CDi << "\n";
    std::cout << "CY (side)          : " << res.CY << "\n";
    std::cout << "L/D (induced only) : " << (res.CDi > 1e-9 ? res.CL / res.CDi : 0.0) << "\n";

    // Sanity check against linear thin-wing theory for an unswept,
    // untwisted flat wing: CL ~= 2*pi*alpha / (1 + 2/AR)  [alpha in rad]
    double alphaRad = fc.alphaDeg * Vlm::Pi / 180.0;
    double CL_flat_theory = 2.0 * Vlm::Pi * alphaRad / (1.0 + 2.0 / AR);
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
    std::vector<Vlm::Panel> panels = Vlm::BuildWing(wing);
    for (auto& p : panels) p.Surface = "wing";
    Vlm::ReferenceGeometry ref;
    ref.Area = res.ReferenceArea; ref.Span = wing.Span; ref.Chord = res.ReferenceArea / wing.Span;
    double trail = 50.0 * wing.Span;
    fc.RefPoint = Vlm::Vec3(ref.Chord * 0.25, 0, 0);
    Vlm::StabilityDerivatives d = Vlm::ComputeDerivatives(panels, fc, ref, trail);

    std::cout << std::setprecision(4);
    std::cout << "\n=== Stability derivatives (wing alone, about x=" << fc.RefPoint.x << ") ===\n";
    std::cout << "CL_alpha  = " << d.CL_alpha << " /rad\n";
    std::cout << "Cm_alpha  = " << d.Cm_alpha << " /rad  (" << (d.Cm_alpha < 0 ? "stable" : "UNSTABLE") << ")\n";
    std::cout << "CY_beta   = " << d.CY_beta << " /rad\n";
    std::cout << "Croll_beta= " << d.Croll_beta << " /rad\n";

    return 0;
}
