// GeometryVlm.cpp — aircraft-from-CAD (STL) VLM driver
//
// Build (from the repo root):
//   g++ -std=c++17 -O2 -Iinclude -o geometry_vlm src/GeometryVlm.cpp \
//       -llapack -lblas
// (or use the top-level CMakeLists.txt, which wires up the LAPACK/OpenBLAS
// backend for you)
//
// Workflow:
//   Export each lifting surface as its own STL file from your CAD tool
//   (recommended -- see the assumptions in MeshSlice.h for why
//   per-component files are much more reliable than auto-segmenting one
//   big STL), then point this program at them.
//
// Usage:
//   ./geometry_vlm --alpha 5 --vinf 25 --rho 1.225
//       --surface name=wing,file=wing.stl,span_axis=y,panels=24
//       --surface name=htail,file=htail.stl,span_axis=y,panels=12
//       --surface name=vtail,file=vtail.stl,span_axis=z,up=y,panels=10
//       --fuselage fuselage.stl
//
// Surface spec keys (comma-separated, no spaces):
//   name=<tag>              required, used for reporting and CSV output
//   file=<path.stl>         required
//   span_axis=x|y|z         default y   (axis the surface spans along)
//   chord_axis=x|y|z        default x   (axis representing chordwise position)
//   thickness_axis=x|y|z    default z   (remaining axis, splits upper/lower skin)
//   up=x|y|z or up=-x|-y|-z default z   (reference direction to disambiguate normal sign)
//   panels=<int>            default 20  (spanwise panels for this surface)
//   cosine=0|1               default 1
//
// The reference area used for CL/CD/CY is the area of the surface named
// "wing" (standard aircraft convention: coefficients are referenced to
// wing planform area even though tail surfaces also generate force and
// are included in the same linear system). If no surface is named "wing",
// the sum of all lifting-surface areas is used instead (a warning is
// printed).
//
// The fuselage, if given, is loaded and its bounding box reported, but it
// is NOT added to the vortex lattice -- classical VLM does not model a
// slender body as a lifting (camber) surface. If fuselage carryover /
// blockage effects matter for your case, that needs a different method
// (e.g. slender-body source/doublet line) layered on top of this.
//
// --- Viscous drag (see DragEstimate.h) ---
// Vlm::Solve() only gives INDUCED drag (CDi) -- a potential-flow method
// has no viscosity. Add --drag-component flags to get a component-buildup
// CD0 estimate, reported alongside a Total CD = CDi + CD0:
//
//   --drag-component surface=wing,tc=0.12,xcmax=0.3,sweep=8,laminar_frac=0,Q=1.0
//   --drag-component surface=htail,tc=0.10,sweep=0,Q=1.04
//   --drag-component surface=fuselage,body=1
//   --drag-component surface=other,label=gear,swet=0.3,length=0.2,body=1,fineness=3,Q=1.2
//
// Drag-component spec keys:
//   surface=<name>|fuselage|other   which loaded surface to attach to
//                                    ("other" = standalone, needs swet=,length=)
//   tc=<thickness ratio>            default 0.12 (airfoil-like components)
//   xcmax=<x/c of max thickness>    default 0.3
//   sweep=<deg>                     default 0 (sweep of max-thickness line)
//   body=0|1                        1 = use fineness-ratio form factor instead of airfoil
//   fineness=<L/D>                  only for surface=other with body=1
//   Q=<interference factor>         default 1.0
//   laminar_frac=<0..1>             default 0 (fully turbulent, conservative)
//   swet=<m^2>, length=<m>          only for surface=other
//   label=<name>                    only for surface=other (display name)
//
// Other drag-related flags: --misc-drag <fraction> (default 0.03, excrescence
// allowance), --mach <M> (compressibility correction to skin friction),
// --mu <Pa*s> (air dynamic viscosity, default sea-level ISA).

#include "Aeolion/Vlm.h"
#include "Aeolion/Mesh.h"
#include "Aeolion/MeshSlice.h"
#include "Aeolion/DragEstimate.h"
#include "Aeolion/VizExport.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Aeolion;
using namespace Aeolion::Vlm;
using MeshIO::PartMesh;

static int AxisFromChar(char c) {
    if (c == 'x' || c == 'X') return 0;
    if (c == 'y' || c == 'Y') return 1;
    if (c == 'z' || c == 'Z') return 2;
    throw std::runtime_error(std::string("bad axis character: ") + c);
}

static Vec3 ParseUp(const std::string& s) {
    bool neg = (!s.empty() && s[0] == '-');
    char c = neg ? s[1] : s[0];
    int ax = AxisFromChar(c);
    Vec3 v(0, 0, 0);
    double val = neg ? -1.0 : 1.0;
    if (ax == 0) v.x = val; else if (ax == 1) v.y = val; else v.z = val;
    return v;
}

struct SurfaceSpec {
    std::string Name;
    std::string File;
    int SpanAxis = 1, ChordAxis = 0, ThicknessAxis = 2;
    Vec3 Up = Vec3(0, 0, 1);
    int Panels = 20;
    bool Cosine = true;
};

static SurfaceSpec ParseSurfaceSpec(const std::string& arg) {
    SurfaceSpec s;
    std::istringstream iss(arg);
    std::string kv;
    while (std::getline(iss, kv, ',')) {
        auto eq = kv.find('=');
        if (eq == std::string::npos) continue;
        std::string k = kv.substr(0, eq), v = kv.substr(eq + 1);
        if (k == "name") s.Name = v;
        else if (k == "file") s.File = v;
        else if (k == "span_axis") s.SpanAxis = AxisFromChar(v[0]);
        else if (k == "chord_axis") s.ChordAxis = AxisFromChar(v[0]);
        else if (k == "thickness_axis") s.ThicknessAxis = AxisFromChar(v[0]);
        else if (k == "up") s.Up = ParseUp(v);
        else if (k == "panels") s.Panels = std::stoi(v);
        else if (k == "cosine") s.Cosine = (v != "0");
        else std::cerr << "warning: unknown surface key '" << k << "'\n";
    }
    if (s.Name.empty() || s.File.empty())
        throw std::runtime_error("surface spec needs at least name=... and file=...");
    return s;
}

struct DragComponentSpec {
    std::string SurfaceName; // matches a loaded --surface name, or "fuselage", or "other"
    double ThicknessRatio = 0.12;
    double xcMaxThickness = 0.3;
    double SweepDeg = 0.0;
    double Q = 1.0;
    double LaminarFraction = 0.0;
    bool IsBody = false;
    double Fineness = 6.0;         // only used if SurfaceName=="other" (standalone component)
    double SwetOverride = -1.0;    // only used if SurfaceName=="other"
    double CharLenOverride = -1.0; // only used if SurfaceName=="other"
    std::string Label;             // only used if SurfaceName=="other"
};

static DragComponentSpec ParseDragSpec(const std::string& arg) {
    DragComponentSpec d;
    std::istringstream iss(arg);
    std::string kv;
    while (std::getline(iss, kv, ',')) {
        auto eq = kv.find('=');
        if (eq == std::string::npos) continue;
        std::string k = kv.substr(0, eq), v = kv.substr(eq + 1);
        if (k == "surface") d.SurfaceName = v;
        else if (k == "tc") d.ThicknessRatio = std::stod(v);
        else if (k == "xcmax") d.xcMaxThickness = std::stod(v);
        else if (k == "sweep") d.SweepDeg = std::stod(v);
        else if (k == "Q") d.Q = std::stod(v);
        else if (k == "laminar_frac") d.LaminarFraction = std::stod(v);
        else if (k == "body") d.IsBody = (v != "0");
        else if (k == "fineness") d.Fineness = std::stod(v);
        else if (k == "swet") d.SwetOverride = std::stod(v);
        else if (k == "length") d.CharLenOverride = std::stod(v);
        else if (k == "label") d.Label = v;
        else std::cerr << "warning: unknown drag-component key '" << k << "'\n";
    }
    if (d.SurfaceName.empty()) throw std::runtime_error("drag-component spec needs surface=...");
    return d;
}

int main(int argc, char** argv) {
    double alphaDeg = 5.0, betaDeg = 0.0, Vinf = 20.0, rho = 1.225;
    double pRate = 0.0, qRate = 0.0, rRate = 0.0;
    Vec3 momentRefPoint(0, 0, 0);
    bool wantDerivatives = false;
    std::vector<SurfaceSpec> specs;
    std::string fuselageFile;
    std::string csvPath = "spanwise_loading.csv";
    std::vector<DragComponentSpec> dragSpecs;
    double miscDragFraction = 0.03;
    double mach = 0.0;
    double airMu = 1.789e-5;
    std::string jsonPath;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : throw std::runtime_error("missing value after " + a); };
        if (a == "--alpha") alphaDeg = std::stod(next());
        else if (a == "--beta") betaDeg = std::stod(next());
        else if (a == "--vinf") Vinf = std::stod(next());
        else if (a == "--rho") rho = std::stod(next());
        else if (a == "--p") pRate = std::stod(next());
        else if (a == "--q") qRate = std::stod(next());
        else if (a == "--r") rRate = std::stod(next());
        else if (a == "--cg") {
            std::string v = next();
            std::istringstream iss(v); std::string tok; std::vector<double> vals;
            while (std::getline(iss, tok, ',')) vals.push_back(std::stod(tok));
            if (vals.size() != 3) throw std::runtime_error("--cg expects x,y,z");
            momentRefPoint = Vec3(vals[0], vals[1], vals[2]);
        }
        else if (a == "--derivatives") wantDerivatives = true;
        else if (a == "--drag-component") dragSpecs.push_back(ParseDragSpec(next()));
        else if (a == "--misc-drag") miscDragFraction = std::stod(next());
        else if (a == "--mach") mach = std::stod(next());
        else if (a == "--mu") airMu = std::stod(next());
        else if (a == "--json") jsonPath = next();
        else if (a == "--surface") specs.push_back(ParseSurfaceSpec(next()));
        else if (a == "--fuselage") fuselageFile = next();
        else if (a == "--csv") csvPath = next();
        else if (a == "--help" || a == "-h") {
            std::cout << "See header comment in GeometryVlm.cpp for full usage.\n";
            return 0;
        } else {
            std::cerr << "Unknown argument: " << a << "\n";
            return 1;
        }
    }

    if (specs.empty()) {
        std::cerr << "No --surface given. Example:\n"
                     "  ./geometry_vlm --surface name=wing,file=wing.stl,span_axis=y,panels=24\n";
        return 1;
    }

    std::vector<Panel> allPanels;
    double wingArea = -1.0;
    double wingSpan = -1.0;
    double maxSpanExtent = 0.0;
    std::map<std::string, double> swetBySurface, charLenBySurface;

    for (const auto& spec : specs) {
        std::vector<PartMesh> parts;
        try {
            parts = MeshIO::ReadStl(spec.File, spec.Name);
        } catch (const std::exception& e) {
            std::cerr << "Failed to read '" << spec.File << "': " << e.what() << "\n";
            return 1;
        }
        if (parts.size() > 1) {
            std::cerr << "Note: '" << spec.File << "' contains " << parts.size()
                      << " named solids; using the first ('" << parts[0].Name << "'). "
                      << "Split multi-solid STL files if you need the others.\n";
        }
        const PartMesh& part = parts[0];

        Vec3 lo, hi; part.Bbox(lo, hi);
        double extent = MeshSlice::AxisComponent(hi, spec.SpanAxis) - MeshSlice::AxisComponent(lo, spec.SpanAxis);
        maxSpanExtent = std::max(maxSpanExtent, extent);

        MeshSlice::SurfaceMeshParams sp;
        sp.SurfaceName = spec.Name;
        sp.SpanAxis = spec.SpanAxis;
        sp.ChordAxis = spec.ChordAxis;
        sp.ThicknessAxis = spec.ThicknessAxis;
        sp.NPanelsSpan = spec.Panels;
        sp.CosineSpacing = spec.Cosine;
        sp.ReferenceUp = spec.Up;

        std::vector<Panel> panels;
        try {
            panels = MeshSlice::MeshToPanels(part, sp);
        } catch (const std::exception& e) {
            std::cerr << "Failed to build panels for surface '" << spec.Name << "': " << e.what() << "\n";
            return 1;
        }

        double area = 0.0;
        for (auto& p : panels) area += p.Area;
        std::cout << "Surface '" << spec.Name << "': " << panels.size() << " panels, area = " << area << " m^2\n";
        if (spec.Name == "wing") { wingArea = area; wingSpan = extent; }

        swetBySurface[spec.Name] = DragEstimate::WettedArea(part);
        charLenBySurface[spec.Name] = (extent > 1e-9) ? area / extent : area; // mean geometric chord

        allPanels.insert(allPanels.end(), panels.begin(), panels.end());
    }

    double fuselageSwet = 0.0, fuselageLength = 0.0, fuselageFineness = 0.0;
    if (!fuselageFile.empty()) {
        try {
            auto fParts = MeshIO::ReadStl(fuselageFile, "fuselage");
            Vec3 lo, hi; fParts[0].Bbox(lo, hi);
            std::cout << "Fuselage '" << fuselageFile << "' loaded (" << fParts[0].Tris.size()
                      << " tris) for reference only -- NOT added to the lattice. "
                      << "bbox (" << lo.x << "," << lo.y << "," << lo.z << ") to ("
                      << hi.x << "," << hi.y << "," << hi.z << ")\n";
            Vec3 ext(hi.x - lo.x, hi.y - lo.y, hi.z - lo.z);
            double extents[3] = {ext.x, ext.y, ext.z};
            int lenAxis = std::max_element(extents, extents + 3) - extents;
            fuselageLength = extents[lenAxis];
            double d1 = extents[(lenAxis + 1) % 3], d2 = extents[(lenAxis + 2) % 3];
            double equivDiameter = std::sqrt(std::max(d1 * d2, 1e-9));
            fuselageFineness = (equivDiameter > 1e-9) ? fuselageLength / equivDiameter : 6.0;
            fuselageSwet = DragEstimate::WettedArea(fParts[0]);
            std::cout << "  length=" << fuselageLength << "  equiv.diameter=" << equivDiameter
                      << "  fineness=" << fuselageFineness << "  Swet=" << fuselageSwet << " m^2\n";
        } catch (const std::exception& e) {
            std::cerr << "Warning: failed to read fuselage file: " << e.what() << "\n";
        }
    }

    double referenceArea = wingArea;
    if (referenceArea <= 0.0) {
        std::cerr << "Warning: no surface named 'wing' found; using total lifting-surface area as reference.\n";
        double S = 0.0; for (auto& p : allPanels) S += p.Area;
        referenceArea = S;
    }
    double trail = 50.0 * std::max(maxSpanExtent, 1e-6);

    ReferenceGeometry ref;
    ref.Area = referenceArea;
    ref.Span = wingSpan > 0.0 ? wingSpan : maxSpanExtent;
    ref.Chord = ref.Span > 0.0 ? referenceArea / ref.Span : referenceArea;

    FreestreamConditions fc;
    fc.Vinf = Vinf; fc.alphaDeg = alphaDeg; fc.betaDeg = betaDeg; fc.rho = rho;
    fc.p = pRate; fc.q = qRate; fc.r = rRate;
    fc.RefPoint = momentRefPoint;

    SolveResult res = Solve(allPanels, fc, ref, trail);

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "\n=== Aircraft VLM result (alpha=" << alphaDeg << " deg, Vinf=" << Vinf << " m/s) ===\n";
    std::cout << "Reference area (wing) : " << res.ReferenceArea << " m^2\n";
    std::cout << "Total panels           : " << res.gamma.size() << "\n";
    std::cout << "CL                     : " << res.CL << "\n";
    std::cout << "CDi                    : " << res.CDi << "\n";
    std::cout << "CY                     : " << res.CY << "\n";

    std::cout << "\nPer-surface breakdown (forces in N; CL_contrib uses the wing reference area/q):\n";
    double q = 0.5 * rho * Vinf * Vinf;
    for (auto& kv : res.LiftBySurface) {
        const std::string& name = kv.first;
        double L = kv.second;
        double D = res.DragBySurface.count(name) ? res.DragBySurface.at(name) : 0.0;
        double A = res.AreaBySurface.count(name) ? res.AreaBySurface.at(name) : 0.0;
        std::cout << "  " << std::left << std::setw(12) << name << std::right
                  << "  area=" << std::setw(9) << A
                  << "  L=" << std::setw(10) << L
                  << "  Di=" << std::setw(10) << D
                  << "  CL_contrib=" << std::setw(9) << (q * res.ReferenceArea > 1e-9 ? L / (q * res.ReferenceArea) : 0.0)
                  << "\n";
    }

    std::ofstream csv(csvPath);
    csv << "surface,y,gamma,lift_per_span_N_m,cl_local\n";
    csv << std::setprecision(8);
    for (const auto& s : res.Stations) {
        csv << s.Surface << "," << s.y << "," << s.gamma << "," << s.LiftPerSpan << "," << s.cl_local << "\n";
    }
    std::cout << "\nSpanwise loading (all surfaces) written to " << csvPath << "\n";

    if (!jsonPath.empty()) {
        std::ofstream jf(jsonPath);
        jf << VizExport::PanelsToJson(allPanels, res, Vinf, rho);
        std::cout << "Panel geometry + scalars written to " << jsonPath << " (for the WebGL viewer)\n";
    }

    if (!dragSpecs.empty()) {
        DragEstimate::AirProperties air; air.rho = rho; air.mu = airMu;
        std::vector<DragEstimate::ComponentSpec> comps;
        for (const auto& ds : dragSpecs) {
            DragEstimate::ComponentSpec c;
            if (ds.SurfaceName == "fuselage") {
                if (fuselageSwet <= 0.0) { std::cerr << "warning: --drag-component surface=fuselage given but no --fuselage loaded; skipping.\n"; continue; }
                c.Name = "fuselage";
                c.Swet = fuselageSwet;
                c.CharLength = fuselageLength;
                c.IsBody = ds.IsBody || true; // fuselage is body-like by default
                c.Fineness = fuselageFineness;
            } else if (ds.SurfaceName == "other") {
                c.Name = ds.Label.empty() ? "other" : ds.Label;
                c.Swet = ds.SwetOverride;
                c.CharLength = ds.CharLenOverride;
                c.IsBody = ds.IsBody;
                c.Fineness = ds.Fineness;
            } else {
                if (!swetBySurface.count(ds.SurfaceName)) { std::cerr << "warning: --drag-component surface='" << ds.SurfaceName << "' doesn't match any loaded --surface; skipping.\n"; continue; }
                c.Name = ds.SurfaceName;
                c.Swet = swetBySurface[ds.SurfaceName];
                c.CharLength = charLenBySurface[ds.SurfaceName];
                c.IsBody = ds.IsBody;
                c.Fineness = ds.Fineness;
            }
            c.ThicknessRatio = ds.ThicknessRatio;
            c.xcMaxThickness = ds.xcMaxThickness;
            c.SweepDeg = ds.SweepDeg;
            c.Q = ds.Q;
            c.LaminarFraction = ds.LaminarFraction;
            comps.push_back(c);
        }

        DragEstimate::BuildupResult drag = DragEstimate::EstimateCD0(comps, Vinf, res.ReferenceArea, air, miscDragFraction, mach);

        std::cout << "\n=== Viscous drag buildup (component method -- conceptual-design accuracy) ===\n";
        std::cout << std::left << std::setw(12) << "component" << std::right
                  << std::setw(14) << "Re" << std::setw(10) << "Cf" << std::setw(10) << "FF"
                  << std::setw(8) << "Q" << std::setw(10) << "Swet" << std::setw(14) << "CD0_contrib" << "\n";
        for (const auto& cr : drag.Components) {
            std::cout << std::left << std::setw(12) << cr.Name << std::right
                      << std::setw(14) << std::scientific << std::setprecision(3) << cr.Re
                      << std::fixed << std::setprecision(5)
                      << std::setw(10) << cr.Cf << std::setw(10) << cr.FF
                      << std::setw(8) << cr.Q << std::setw(10) << cr.Swet
                      << std::setw(14) << cr.CD0_contribution << "\n";
        }
        std::cout << "CD0 (clean, before misc allowance) : " << drag.CD0_clean << "\n";
        std::cout << "CD0 (with " << (miscDragFraction * 100) << "% misc/excrescence allowance) : " << drag.CD0 << "\n";
        std::cout << "\nTotal CD = CDi + CD0 = " << res.CDi << " + " << drag.CD0 << " = " << (res.CDi + drag.CD0) << "\n";
        std::cout << "L/D (total) = " << (res.CL / (res.CDi + drag.CD0)) << "\n";
    } else {
        std::cout << "\n(No --drag-component given -- CDi above is INDUCED drag only; "
                     "see the --drag-component flag to add a viscous CD0 estimate for total CD.)\n";
    }

    if (wantDerivatives) {
        std::cout << "\n=== Stability & control derivatives about CG=("
                  << momentRefPoint.x << "," << momentRefPoint.y << "," << momentRefPoint.z
                  << "), baseline alpha=" << alphaDeg << " deg, beta=" << betaDeg << " deg ===\n";
        StabilityDerivatives d = ComputeDerivatives(allPanels, fc, ref, trail);
        std::cout << "Baseline: CL0=" << d.CL0 << "  CDi0=" << d.CDi0 << "  CY0=" << d.CY0
                  << "  Cm0=" << d.Cm0 << "  Croll0=" << d.Croll0 << "  Cn0=" << d.Cn0 << "\n";
        std::cout << "Longitudinal (per rad):  CL_alpha=" << d.CL_alpha
                  << "  CDi_alpha=" << d.CDi_alpha << "  Cm_alpha=" << d.Cm_alpha << "\n";
        std::cout << "  -> static margin sign: " << (d.Cm_alpha < 0 ? "stable (Cm_alpha < 0)" : "UNSTABLE (Cm_alpha >= 0)") << "\n";
        std::cout << "Lateral-directional (per rad):  CY_beta=" << d.CY_beta
                  << "  Croll_beta=" << d.Croll_beta << "  Cn_beta=" << d.Cn_beta << "\n";
        std::cout << "  -> Cn_beta " << (d.Cn_beta > 0 ? "> 0 (weathercock-stable)" : "<= 0 (weathercock-UNSTABLE)")
                  << ",  Croll_beta " << (d.Croll_beta < 0 ? "< 0 (dihedral-stable)" : ">= 0 (dihedral-destabilizing)") << "\n";
        std::cout << "Pitch rate:     Cm_q=" << d.Cm_q << " /(rad/s)   Cm_q_nd=" << d.Cm_q_nd << " (x cbar/2V)\n";
        std::cout << "Roll rate:      Croll_p=" << d.Croll_p << " /(rad/s)   Croll_p_nd=" << d.Croll_p_nd << " (x b/2V)  [roll damping, expect < 0]\n";
        std::cout << "Yaw rate:       Cn_r=" << d.Cn_r << " /(rad/s)   Cn_r_nd=" << d.Cn_r_nd << " (x b/2V)  [yaw damping, expect < 0]\n";
        std::cout << "Cross terms:    Cn_p=" << d.Cn_p << "   Croll_r=" << d.Croll_r << "   CY_p=" << d.CY_p << "   CY_r=" << d.CY_r << "\n";
    }

    return 0;
}
