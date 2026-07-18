// DragEstimate.h
//
// Viscous ("profile"/"parasite") drag estimate via the classic component
// buildup method (see e.g. Raymer, "Aircraft Design: A Conceptual
// Approach", ch. 12; or Hoerner, "Fluid-Dynamic Drag"). This is NOT
// something a vortex lattice method can give you -- VLM is a potential-flow
// (inviscid) method, so VLM::SolveResult::CDi is induced drag only. Total
// drag is:
//
//     CD_total = CDi (from VLM::Solve)  +  CD0 (from this file)
//
// CD0 is a mid-fidelity CONCEPTUAL-DESIGN-LEVEL estimate, not a substitute
// for wind-tunnel or CFD data -- typical accuracy is +/-10-20% for a
// clean configuration. It gets noticeably less reliable with flow
// separation, flap deflection, interference-heavy junctions, or anything
// transonic.
//
// Method, per wetted component c:
//
//     CD0 = sum_c( Cf_c * FF_c * Q_c * Swet_c ) / Sref  *  (1 + miscFraction)
//
//   Cf_c   : flat-plate skin friction coefficient at the component's local
//            Reynolds number (turbulent, laminar, or a simple mixed blend)
//   FF_c   : form factor -- how much extra (mostly pressure) drag the
//            component's real shape adds over an idealized flat plate,
//            from thickness ratio + sweep (airfoil-like components) or
//            fineness ratio (body-like components, e.g. a fuselage)
//   Q_c    : interference factor -- extra drag from how the component
//            junctions with the rest of the airframe (1.0 = isolated)
//   Swet_c : wetted area of the component
//   miscFraction : lumped allowance for excrescences, leakage, gaps,
//            antennas etc. (2-6% is a common conceptual-design allowance;
//            0 if you'd rather account for these yourself)
//
// No external dependencies beyond the standard library, Aeolion/Math
// (Vec3/Cross) and Aeolion/Mesh (PartMesh, to compute wetted area from a
// triangle mesh).

#pragma once
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <numbers>
#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Mesh.h"

namespace Aeolion::DragEstimate {

using VLM::Vec3;
using VLM::Cross;

struct AirProperties {
    double rho = 1.225;      // [kg/m^3] -- sea-level ISA by default
    double mu = 1.789e-5;    // [Pa*s]   -- sea-level ISA dynamic viscosity
};

inline double ReynoldsNumber(double Vinf, double charLength, const AirProperties& air) {
    return air.rho * Vinf * charLength / air.mu;
}

// Turbulent flat-plate skin friction, Prandtl-Schlichting correlation, with
// an optional Frankl-Voishel compressibility correction.
inline double CfTurbulent(double Re, double mach = 0.0) {
    if (Re < 1.0) return 0.0;
    double cf = 0.455 / std::pow(std::log10(Re), 2.58);
    if (mach > 1e-6) cf /= std::pow(1.0 + 0.144 * mach * mach, 0.65);
    return cf;
}

// Laminar flat-plate skin friction (Blasius).
inline double CfLaminar(double Re) {
    if (Re < 1.0) return 0.0;
    return 1.328 / std::sqrt(Re);
}

// Simple engineering blend between laminar and turbulent Cf, weighted by
// the fraction of the run assumed laminar before transition. This is an
// approximation (not the rigorous laminar-run-subtracted method) but is
// adequate for conceptual estimates. laminarFraction=0 (fully turbulent)
// is the conservative default -- most real surfaces trip well forward of
// the ideal transition point due to roughness, insects, rivets, etc.
// Genuinely smooth/glider-polished laminar-flow airfoils might justify
// 0.3-0.5; a rough or heavily-riveted metal surface should use 0.
inline double CfMixed(double Re, double laminarFraction, double mach = 0.0) {
    laminarFraction = std::clamp(laminarFraction, 0.0, 1.0);
    if (laminarFraction <= 0.0) return CfTurbulent(Re, mach);
    return laminarFraction * CfLaminar(Re) + (1.0 - laminarFraction) * CfTurbulent(Re, mach);
}

// Airfoil-like (wing/tail) form factor: thickness + sweep effect on
// pressure drag, beyond flat-plate friction. xcMaxThickness is the
// chordwise location of maximum thickness (~0.3 for most conventional
// airfoils, ~0.4-0.5 for laminar-flow sections); sweepDeg is the sweep of
// the max-thickness line (quarter-chord sweep is a fine approximation).
inline double FormFactorAirfoil(double thicknessRatio, double xcMaxThickness, double sweepDeg) {
    double tc = thicknessRatio;
    double xm = std::max(xcMaxThickness, 1e-3);
    double sweep = sweepDeg * std::numbers::pi / 180.0;
    double base = 1.0 + (0.6 / xm) * tc + 100.0 * std::pow(tc, 4);
    return base * std::pow(std::cos(sweep), 0.28);
}

// Body-like (fuselage, pod, boom) form factor from fineness ratio
// f = length / equivalent_diameter.
inline double FormFactorBody(double fineness) {
    double f = std::max(fineness, 1e-3);
    return 1.0 + 60.0 / (f * f * f) + f / 400.0;
}

// Sum of triangle areas -- for a thin-shell surface mesh (upper+lower skin)
// this IS the wetted area directly, no factor-of-2 needed since both skins
// are already separate triangles in the mesh.
inline double WettedArea(const MeshIO::PartMesh& part) {
    double A = 0.0;
    for (const auto& t : part.Tris) {
        Vec3 e1 = t.v1 - t.v0, e2 = t.v2 - t.v0;
        A += 0.5 * Cross(e1, e2).Norm();
    }
    return A;
}

struct ComponentSpec {
    std::string Name;
    double Swet = 0.0;           // wetted area [m^2]
    double CharLength = 0.0;     // Reynolds-number reference length [m] (chord, or fuselage length)
    bool IsBody = false;         // true: use FormFactorBody (needs fineness); false: FormFactorAirfoil (needs thickness/sweep)
    double ThicknessRatio = 0.12;
    double xcMaxThickness = 0.3;
    double SweepDeg = 0.0;
    double Fineness = 6.0;
    double Q = 1.0;              // interference factor (1.0 = isolated; 1.03-1.08 typical for a wing/tail-fuselage junction)
    double LaminarFraction = 0.0;
};

struct ComponentResult {
    std::string Name;
    double Re = 0.0, Cf = 0.0, FF = 0.0, Q = 0.0, Swet = 0.0;
    double CD0_contribution = 0.0; // already divided by Sref
};

inline ComponentResult EstimateComponent(const ComponentSpec& c, double Vinf, double Sref,
                                          const AirProperties& air, double mach = 0.0) {
    ComponentResult r;
    r.Name = c.Name;
    r.Re = ReynoldsNumber(Vinf, c.CharLength, air);
    r.Cf = CfMixed(r.Re, c.LaminarFraction, mach);
    r.FF = c.IsBody ? FormFactorBody(c.Fineness) : FormFactorAirfoil(c.ThicknessRatio, c.xcMaxThickness, c.SweepDeg);
    r.Q = c.Q;
    r.Swet = c.Swet;
    r.CD0_contribution = (Sref > 1e-9) ? (r.Cf * r.FF * r.Q * r.Swet / Sref) : 0.0;
    return r;
}

struct BuildupResult {
    std::vector<ComponentResult> Components;
    double CD0_clean = 0.0;   // sum of component contributions, before misc allowance
    double MiscFraction = 0.0;
    double CD0 = 0.0;         // CD0_clean * (1 + MiscFraction) -- use this one
};

inline BuildupResult EstimateCD0(const std::vector<ComponentSpec>& specs, double Vinf, double Sref,
                                  const AirProperties& air, double miscFraction = 0.03, double mach = 0.0) {
    BuildupResult res;
    res.MiscFraction = miscFraction;
    for (const auto& c : specs) {
        auto cr = EstimateComponent(c, Vinf, Sref, air, mach);
        res.Components.push_back(cr);
        res.CD0_clean += cr.CD0_contribution;
    }
    res.CD0 = res.CD0_clean * (1.0 + miscFraction);
    return res;
}

} // namespace Aeolion::DragEstimate
