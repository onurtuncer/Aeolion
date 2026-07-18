// VizExport.h
//
// Exports a solved VLM panel set as JSON: one quad per panel (reconstructed
// from the bound-vortex line + local chordwise direction, since Panel only
// stores the quarter-chord line and control point, not full LE/TE corners)
// plus per-panel scalar fields for coloring. No VTK/ParaView dependency --
// meant to feed a raw WebGL viewer directly.
//
// Quad reconstruction: each panel's bound vortex (A->B) sits at the local
// quarter-chord; the control point sits at the local three-quarter-chord.
// chordDir = normalize(ControlPoint - midpoint(A,B)) approximates the local
// chordwise (aft) direction, and chord = panel.Area / panel.SpanwiseWidth
// is the local mean chord. LE and TE corners are then A/B offset by
// -0.25*chord and +0.75*chord along chordDir. This is an approximation
// (treats chord direction/length as uniform across the panel's span) but
// is accurate enough for visualization.

#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "Aeolion/Vlm.h"

namespace Aeolion { namespace VizExport {

using Vlm::Vec3;
using Vlm::Panel;
using Vlm::SolveResult;

inline std::string PanelsToJson(const std::vector<Panel>& panels, const SolveResult& res,
                                 double Vinf, double rho) {
    std::ostringstream j;
    j << std::setprecision(7);
    j << "{\n \"panels\": [\n";

    double q = 0.5 * rho * Vinf * Vinf;

    for (size_t i = 0; i < panels.size(); ++i) {
        const Panel& p = panels[i];
        Vec3 mid = (p.A + p.B) * 0.5;
        Vec3 chordDir = (p.ControlPoint - mid).Normalized();
        double chord = (p.SpanwiseWidth > 1e-9) ? p.Area / p.SpanwiseWidth : 0.0;

        Vec3 leA = p.A - chordDir * (0.25 * chord);
        Vec3 teA = p.A + chordDir * (0.75 * chord);
        Vec3 leB = p.B - chordDir * (0.25 * chord);
        Vec3 teB = p.B + chordDir * (0.75 * chord);

        double gamma = (i < res.gamma.size()) ? res.gamma[i] : 0.0;
        // approximate local Cp-like loading indicator: normalize circulation
        // by a reference (Vinf * mean chord) to get a dimensionless number
        // comparable across panels of different size.
        double refChord = (chord > 1e-9) ? chord : 1.0;
        double loadingCoeff = (Vinf > 1e-9) ? gamma / (Vinf * refChord) : 0.0;
        (void)q;

        j << "  {\"i\":" << i << ",\"surface\":\"" << p.Surface << "\","
          << "\"le_a\":[" << leA.x << "," << leA.y << "," << leA.z << "],"
          << "\"te_a\":[" << teA.x << "," << teA.y << "," << teA.z << "],"
          << "\"te_b\":[" << teB.x << "," << teB.y << "," << teB.z << "],"
          << "\"le_b\":[" << leB.x << "," << leB.y << "," << leB.z << "],"
          << "\"normal\":[" << p.Normal.x << "," << p.Normal.y << "," << p.Normal.z << "],"
          << "\"gamma\":" << gamma << ",\"loading\":" << loadingCoeff
          << "}" << (i + 1 < panels.size() ? "," : "") << "\n";
    }
    j << " ],\n";

    j << " \"stations\": [\n";
    for (size_t k = 0; k < res.Stations.size(); ++k) {
        const auto& s = res.Stations[k];
        j << "  {\"panel\":" << s.PanelIndex << ",\"y\":" << s.y
          << ",\"cl_local\":" << s.cl_local << ",\"lift_per_span\":" << s.LiftPerSpan
          << ",\"surface\":\"" << s.Surface << "\"}" << (k + 1 < res.Stations.size() ? "," : "") << "\n";
    }
    j << " ],\n";

    j << " \"summary\": {\"CL\":" << res.CL << ",\"CDi\":" << res.CDi << ",\"CY\":" << res.CY
      << ",\"referenceArea\":" << res.ReferenceArea << "}\n";
    j << "}\n";
    return j.str();
}

}} // namespace Aeolion::VizExport
