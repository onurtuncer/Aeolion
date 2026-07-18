// AirfoilPolar.cpp — 2D viscous airfoil polar tool
//
// Combines AirfoilPanel.h (inviscid surface velocity distribution) with
// BoundaryLayer.h (Thwaites/Michel/Head/Squire-Young) to estimate
// section cd for a given cl and Reynolds number -- a much higher-fidelity
// viscous drag estimate than the flat-plate + form-factor buildup in
// DragEstimate.h, because it uses the REAL pressure-gradient-driven
// velocity distribution around the actual airfoil shape, including
// transition location and a separation warning.
//
// Build:
//   g++ -std=c++17 -O2 -o airfoil_polar AirfoilPolar.cpp
//
// Usage:
//   ./airfoil_polar --naca 2412 --re 500000 --cl 0.2,0.4,0.6,0.8
//   ./airfoil_polar --naca 0012 --re 1000000 --alpha -2,0,2,4,6,8
//   ./airfoil_polar --dat myfoil.dat --re 300000 --cl 0.3
//
// --dat file format: plain list of "x y" pairs (one airfoil-fraction
// coordinate per line, x in [0,1]), ordered lower-surface TE->LE then
// upper-surface LE->TE (Selig-style .dat convention) -- same convention
// Naca4() produces internally.
//
// NOTE on fidelity: this is an UNCOUPLED viscous-inviscid analysis (the
// boundary layer is solved once on the fixed inviscid Ue(s), it does not
// feed back to displace/correct the inviscid solution). That's the same
// core method XFOIL uses, minus its Newton-coupled iteration -- good for
// attached or mildly-loaded flow, less reliable near stall or with
// significant separation (a separation warning is printed when the
// relevant criterion is crossed, but the drag number past that point
// should not be trusted).

#include "Aeolion/AirfoilPanel.h"
#include "Aeolion/BoundaryLayer.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <algorithm>

using namespace Aeolion;
using namespace Aeolion::Airfoil;

std::vector<Vec2> ReadDatFile(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::vector<Vec2> pts;
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream iss(line);
        double x, y;
        if (iss >> x >> y) pts.push_back({x, y});
    }
    if (pts.size() < 10) throw std::runtime_error(".dat file has too few points: " + path);
    return pts;
}

// Given panel solve results (dimensionless, Vinf magnitude=1) and the
// physical flow conditions, find the stagnation point, split into
// upper/lower arc-length branches, and run the BL analysis.
struct PolarPoint {
    double alphaDeg = 0, cl = 0, cd = 0;
    double sTransUpper = -1, sTransLower = -1;
    bool Separated = false;
};

PolarPoint EvaluateAtAlpha(const std::vector<Panel>& panels, double alphaRad,
                            double Vinf, double nu, double chord) {
    SolveResult r = Solve(panels, alphaRad);
    int N = static_cast<int>(panels.size());

    // Cumulative arc length along the loop from panel 0's start point.
    std::vector<double> sCum(N + 1, 0.0);
    for (int i = 0; i < N; ++i) sCum[i + 1] = sCum[i] + panels[i].Length;

    // Find stagnation point: sign change of Ue (r.Ue, dimensionless, signed
    // along each panel's own tangent direction). Search the whole loop
    // (usually near the lower/upper boundary, but don't assume).
    int stagIdx = -1;
    for (int i = 0; i + 1 < N; ++i) {
        if (r.Ue[i] <= 0.0 && r.Ue[i + 1] > 0.0) { stagIdx = i; break; }
    }
    if (stagIdx < 0) {
        // fall back: also check wraparound
        if (r.Ue[N - 1] <= 0.0 && r.Ue[0] > 0.0) stagIdx = N - 1;
    }
    if (stagIdx < 0) throw std::runtime_error("EvaluateAtAlpha: could not locate a stagnation point");

    // Lower branch: walk BACKWARD from stagIdx down to panel 0 (this
    // traverses stagnation -> lower TE, physically the flow direction on
    // the lower surface).
    std::vector<double> sLower, UeLower;
    {
        double acc = 0.0;
        sLower.push_back(0.0);
        UeLower.push_back(std::fabs(r.Ue[stagIdx]) * Vinf);
        for (int i = stagIdx; i >= 0; --i) {
            acc += panels[i].Length;
            sLower.push_back(acc);
            UeLower.push_back(std::fabs(r.Ue[i]) * Vinf);
        }
    }
    // Upper branch: walk FORWARD from stagIdx+1 to panel N-1 (stagnation
    // -> upper TE).
    std::vector<double> sUpper, UeUpper;
    {
        double acc = 0.0;
        sUpper.push_back(0.0);
        UeUpper.push_back(std::fabs(r.Ue[stagIdx + 1]) * Vinf);
        for (int i = stagIdx + 1; i < N; ++i) {
            acc += panels[i].Length;
            sUpper.push_back(acc);
            UeUpper.push_back(std::fabs(r.Ue[i]) * Vinf);
        }
    }

    BoundaryLayer::AirfoilBLResult blr = BoundaryLayer::AnalyzeAirfoil(sUpper, UeUpper, sLower, UeLower, Vinf, nu, chord);

    PolarPoint pp;
    pp.alphaDeg = alphaRad * 180.0 / Pi;
    pp.cl = r.Cl / chord; // Solve()'s Cl is for the geometry as-built; Naca4() already scales by chord, so
                          // if chord!=1 was used consistently in Naca4()+panel build, r.Cl is already correct
                          // per-chord -- guard divide only matters if caller passes mismatched chord.
    pp.cd = blr.cd;
    pp.sTransUpper = blr.Upper.sTransition;
    pp.sTransLower = blr.Lower.sTransition;
    pp.Separated = blr.Upper.Separated || blr.Lower.Separated;
    return pp;
}

int main(int argc, char** argv) {
    std::string nacaCode, datFile;
    double Re = 500000.0;
    double chord = 1.0;
    std::vector<double> targetCl, targetAlpha;
    int panelsPerSide = 140;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : throw std::runtime_error("missing value after " + a); };
        auto parseList = [&](const std::string& s) {
            std::vector<double> v; std::istringstream iss(s); std::string tok;
            while (std::getline(iss, tok, ',')) v.push_back(std::stod(tok));
            return v;
        };
        if (a == "--naca") nacaCode = next();
        else if (a == "--dat") datFile = next();
        else if (a == "--re") Re = std::stod(next());
        else if (a == "--chord") chord = std::stod(next());
        else if (a == "--panels") panelsPerSide = std::stoi(next());
        else if (a == "--cl") targetCl = parseList(next());
        else if (a == "--alpha") targetAlpha = parseList(next());
        else if (a == "--help" || a == "-h") {
            std::cout << "See header comment in AirfoilPolar.cpp for usage.\n";
            return 0;
        } else { std::cerr << "Unknown argument: " << a << "\n"; return 1; }
    }
    if (nacaCode.empty() && datFile.empty()) {
        std::cerr << "Need --naca <4digits> or --dat <file>\n";
        return 1;
    }
    if (targetCl.empty() && targetAlpha.empty()) {
        std::cerr << "Need --cl <list> or --alpha <list>\n";
        return 1;
    }

    std::vector<Vec2> loop = nacaCode.empty() ? ReadDatFile(datFile) : Naca4(nacaCode, panelsPerSide, chord);
    std::vector<Panel> panels = BuildPanels(loop);

    double Vinf = 20.0; // reference speed; only Re (via nu) and chord matter physically for cd(cl)
    double nu = Vinf * chord / Re;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Airfoil: " << (nacaCode.empty() ? datFile : ("NACA " + nacaCode))
              << "   chord=" << chord << " m   Re=" << std::scientific << Re << std::fixed
              << "   panels=" << panels.size() << "\n\n";
    std::cout << std::left << std::setw(10) << "alpha" << std::setw(10) << "cl" << std::setw(10) << "cd"
              << std::setw(10) << "cl/cd" << std::setw(14) << "xtr_upper" << std::setw(14) << "xtr_lower"
              << "note" << "\n";

    // linear-theory: two solves give Cl(alpha) slope+intercept exactly (panel method is linear in alpha)
    double a1 = 2.0 * Pi / 180.0, a2 = -2.0 * Pi / 180.0;
    SolveResult r1 = Solve(panels, a1), r2 = Solve(panels, a2);
    double clSlope = (r1.Cl - r2.Cl) / (a1 - a2);
    double clAtZero = r1.Cl - clSlope * a1;

    std::vector<double> alphasToRun = targetAlpha;
    for (double clT : targetCl) {
        double alphaRad = (clT - clAtZero) / clSlope;
        alphasToRun.push_back(alphaRad * 180.0 / Pi);
    }
    std::sort(alphasToRun.begin(), alphasToRun.end());

    std::ofstream csv("polar.csv");
    csv << "alpha_deg,cl,cd,xtr_upper,xtr_lower,separated\n";

    for (double aDeg : alphasToRun) {
        double aRad = aDeg * Pi / 180.0;
        try {
            PolarPoint pp = EvaluateAtAlpha(panels, aRad, Vinf, nu, chord);
            std::string note = pp.Separated ? "SEPARATION WARNING" : "";
            std::cout << std::left << std::setw(10) << pp.alphaDeg << std::setw(10) << pp.cl
                      << std::setw(10) << pp.cd << std::setw(10) << (pp.cd > 1e-9 ? pp.cl / pp.cd : 0.0)
                      << std::setw(14) << pp.sTransUpper << std::setw(14) << pp.sTransLower
                      << note << "\n";
            csv << pp.alphaDeg << "," << pp.cl << "," << pp.cd << "," << pp.sTransUpper << ","
                << pp.sTransLower << "," << (pp.Separated ? 1 : 0) << "\n";
        } catch (const std::exception& e) {
            std::cerr << "alpha=" << aDeg << ": failed (" << e.what() << ")\n";
        }
    }
    std::cout << "\nPolar written to polar.csv\n";
    std::cout << "(xtr_upper/xtr_lower are transition arc-length from stagnation, in meters; -1 = stayed laminar to TE)\n";

    return 0;
}
