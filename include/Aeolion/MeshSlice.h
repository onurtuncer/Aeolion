// MeshSlice.h
//
// Converts a triangulated CAD/mesh lifting surface (wing, horizontal tail,
// vertical tail -- read from STL, or from STEP via the step2mesh
// converter) into a VLM::Panel lattice, by:
//
//   1. Slicing the mesh with a series of planes perpendicular to the
//      surface's dominant span axis.
//   2. At each station, collecting the mesh/plane intersection points and
//      splitting them into "upper skin" / "lower skin" by height relative
//      to the station's own mean height (robust for the thin, single
//      closed-loop cross-sections typical of a wing/tail skin -- this
//      deliberately avoids reconstructing full loop topology, which is a
//      much harder problem on imperfect meshes).
//   3. Interpolating upper/lower skin height at the local quarter-chord
//      and three-quarter-chord chordwise stations and averaging them to
//      get the mean camber line -- this is where sweep, taper, and twist
//      all fall out automatically, because they're just read from the
//      actual geometry rather than re-derived analytically.
//   4. Building one horseshoe-vortex panel per pair of adjacent stations,
//      exactly like the parametric BuildWing() path, so the rest of
//      VLM.h (influence matrix, solve, near-field forces) is untouched.
//
// IMPORTANT ASSUMPTIONS (read before trusting results on a new geometry):
//   - The component is a single, simply-connected, watertight-ish thin
//     shell (e.g. a wing's OML skin) with one closed cross-section loop
//     per span station. Multi-element surfaces (e.g. a wing with a
//     separate flap/slat body) will NOT slice correctly as one part --
//     mesh them and pass them as separate PartMesh/surfaces instead.
//   - The fuselage is NOT a lifting surface in classical VLM and should
//     not be passed through this function; load it (for bookkeeping /
//     visualization) but exclude it from the panel lattice.
//   - Camber is captured (via upper/lower averaging) but this is still a
//     thin-surface, linear, attached-flow method -- no thickness effects,
//     no separation, no compressibility.

#pragma once
#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <numbers>
#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Types/Panel.h"
#include "Aeolion/Mesh.h"

namespace Aeolion::MeshSlice {

using VLM::Vec3;
using MeshIO::PartMesh;

inline double AxisComponent(const Vec3& v, int axis) {
    return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
}
inline void SetAxisComponent(Vec3& v, int axis, double val) {
    if (axis == 0) v.x = val; else if (axis == 1) v.y = val; else v.z = val;
}

struct SurfaceMeshParams {
    std::string SurfaceName = "wing";
    int SpanAxis = 1;         // 0=x, 1=y, 2=z -- axis the surface spans along
    int ChordAxis = 0;        // axis representing chordwise position
    int ThicknessAxis = 2;    // remaining axis, used to split upper/lower skin
    int NPanelsSpan = 20;     // number of spanwise panels (stations = NPanelsSpan+1)
    bool CosineSpacing = true;
    double TipClipFraction = 0.003; // stay this far inside the raw bbox extent to avoid tip-cap facets
    Vec3 ReferenceUp = Vec3(0, 0, 1); // disambiguates panel-normal sign (e.g. use (0,1,0) for a vertical tail)
};

struct StationCut {
    double Coord = 0.0;         // SpanAxis coordinate of this station
    Vec3 QuarterChordPt;
    Vec3 ThreeQuarterChordPt;
    double Chord = 0.0;
    bool Ok = false;
};

// Slice the mesh at SpanAxis == sVal and extract the mean-camber quarter-
// and three-quarter-chord points.
inline StationCut SliceStation(const PartMesh& part, const SurfaceMeshParams& sp, double sVal) {
    StationCut result;
    result.Coord = sVal;

    struct Pt { double chord, thick; };
    std::vector<Pt> pts;

    auto tryEdge = [&](const Vec3& p0, const Vec3& p1) {
        double a = AxisComponent(p0, sp.SpanAxis) - sVal;
        double b = AxisComponent(p1, sp.SpanAxis) - sVal;
        if ((a > 0 && b > 0) || (a < 0 && b < 0)) return; // no crossing
        if (a == b) return; // degenerate (edge lies in-plane); ignore, other edges of same tri will catch it
        double t = a / (a - b);
        Vec3 p = p0 + (p1 - p0) * t;
        pts.push_back({AxisComponent(p, sp.ChordAxis), AxisComponent(p, sp.ThicknessAxis)});
    };

    for (const auto& tri : part.Tris) {
        tryEdge(tri.v0, tri.v1);
        tryEdge(tri.v1, tri.v2);
        tryEdge(tri.v2, tri.v0);
    }

    if (pts.size() < 4) { result.Ok = false; return result; } // need at least ~2 upper + ~2 lower

    double chordMin = pts[0].chord, chordMax = pts[0].chord;
    for (auto& p : pts) { chordMin = std::min(chordMin, p.chord); chordMax = std::max(chordMax, p.chord); }
    double chord = chordMax - chordMin;
    if (chord < 1e-9) { result.Ok = false; return result; }

    // Least-squares line thick ~= m*chord + b through ALL points first, to
    // remove the mean incidence/twist trend before splitting into skins.
    // A flat mean (thick >= meanThick) misclassifies points once the
    // section's rigid tilt becomes comparable to its thickness (exactly
    // the regime near the leading/trailing edge of a twisted section,
    // where thickness itself goes to ~0).
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (auto& p : pts) { sx += p.chord; sy += p.thick; sxx += p.chord * p.chord; sxy += p.chord * p.thick; }
    double n = (double)pts.size();
    double denom = n * sxx - sx * sx;
    double m = (std::fabs(denom) > 1e-12) ? (n * sxy - sx * sy) / denom : 0.0;
    double b = (sy - m * sx) / n;

    std::vector<Pt> upper, lower;
    for (auto& p : pts) {
        double trend = m * p.chord + b;
        (p.thick >= trend ? upper : lower).push_back(p);
    }
    if (upper.size() < 2 || lower.size() < 2) { result.Ok = false; return result; }

    std::sort(upper.begin(), upper.end(), [](const Pt& a, const Pt& b) { return a.chord < b.chord; });
    std::sort(lower.begin(), lower.end(), [](const Pt& a, const Pt& b) { return a.chord < b.chord; });

    auto interp = [](const std::vector<Pt>& skin, double x) -> double {
        if (x <= skin.front().chord) return skin.front().thick;
        if (x >= skin.back().chord) return skin.back().thick;
        for (size_t i = 0; i + 1 < skin.size(); ++i) {
            if (x >= skin[i].chord && x <= skin[i + 1].chord) {
                double dx = skin[i + 1].chord - skin[i].chord;
                if (dx < 1e-12) return skin[i].thick;
                double t = (x - skin[i].chord) / dx;
                return skin[i].thick + t * (skin[i + 1].thick - skin[i].thick);
            }
        }
        return skin.back().thick;
    };

    double xq = chordMin + 0.25 * chord;
    double x3 = chordMin + 0.75 * chord;

    double zq = 0.5 * (interp(upper, xq) + interp(lower, xq));
    double z3 = 0.5 * (interp(upper, x3) + interp(lower, x3));

    Vec3 pq, p3;
    SetAxisComponent(pq, sp.SpanAxis, sVal); SetAxisComponent(pq, sp.ChordAxis, xq); SetAxisComponent(pq, sp.ThicknessAxis, zq);
    SetAxisComponent(p3, sp.SpanAxis, sVal); SetAxisComponent(p3, sp.ChordAxis, x3); SetAxisComponent(p3, sp.ThicknessAxis, z3);

    result.QuarterChordPt = pq;
    result.ThreeQuarterChordPt = p3;
    result.Chord = chord;
    result.Ok = true;
    return result;
}

// Build the full panel lattice for one lifting-surface mesh.
inline std::vector<VLM::Panel> MeshToPanels(const PartMesh& part, const SurfaceMeshParams& sp) {
    Vec3 lo, hi;
    part.Bbox(lo, hi);
    double sMin = AxisComponent(lo, sp.SpanAxis);
    double sMax = AxisComponent(hi, sp.SpanAxis);
    double span = sMax - sMin;
    if (span < 1e-9) throw std::runtime_error("MeshToPanels: '" + part.Name + "' has ~zero extent along span axis " + std::to_string(sp.SpanAxis));

    double clip = sp.TipClipFraction * span;
    double s0 = sMin + clip, s1 = sMax - clip;

    int N = sp.NPanelsSpan;
    std::vector<double> stationCoord(N + 1);
    for (int i = 0; i <= N; ++i) {
        double u = (double)i / N; // 0..1
        if (sp.CosineSpacing) {
            double theta = u * std::numbers::pi;
            double uc = 0.5 * (1.0 - std::cos(theta)); // clusters near both ends
            stationCoord[i] = s0 + uc * (s1 - s0);
        } else {
            stationCoord[i] = s0 + u * (s1 - s0);
        }
    }

    std::vector<StationCut> cuts(N + 1);
    for (int i = 0; i <= N; ++i) {
        cuts[i] = SliceStation(part, sp, stationCoord[i]);
        if (!cuts[i].Ok) {
            std::ostringstream oss;
            oss << "MeshToPanels: slice failed for surface '" << sp.SurfaceName
                << "' at station " << i << " (SpanAxis coord = " << stationCoord[i] << "). "
                << "Check SpanAxis/ChordAxis/ThicknessAxis assignment, or that this component "
                << "is a single simply-connected thin-shell surface.";
            throw std::runtime_error(oss.str());
        }
    }

    std::vector<VLM::Panel> panels;
    panels.reserve(N);
    for (int i = 0; i < N; ++i) {
        const StationCut& c0 = cuts[i];
        const StationCut& c1 = cuts[i + 1];

        VLM::Panel p;
        p.A = c0.QuarterChordPt;
        p.B = c1.QuarterChordPt;
        p.ControlPoint = (c0.ThreeQuarterChordPt + c1.ThreeQuarterChordPt) * 0.5;

        Vec3 chordTan0 = (c0.ThreeQuarterChordPt - c0.QuarterChordPt).Normalized();
        Vec3 chordTan1 = (c1.ThreeQuarterChordPt - c1.QuarterChordPt).Normalized();
        Vec3 chordTan = (chordTan0 + chordTan1).Normalized();
        Vec3 spanTan = (p.B - p.A).Normalized();

        Vec3 n = VLM::Cross(spanTan, chordTan).Normalized();
        if (VLM::Dot(n, sp.ReferenceUp) < 0.0) n = -n;
        p.Normal = n;

        p.TrailDirA = Vec3(1, 0, 0);
        p.TrailDirB = Vec3(1, 0, 0);
        p.SpanwiseWidth = std::fabs(stationCoord[i + 1] - stationCoord[i]);
        p.Area = 0.5 * (c0.Chord + c1.Chord) * p.SpanwiseWidth;
        p.Surface = sp.SurfaceName;

        panels.push_back(p);
    }
    return panels;
}

} // namespace Aeolion::MeshSlice
