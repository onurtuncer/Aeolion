// LatticeFigureExport.cpp -- dump the paper's reference configuration
// (two-blade rotor, snug shroud, four-vane cruciform; the same case
// TestVaneCascade's coupled solve pins) as JSON for the publication
// figures: panel quads plus the solved circulation, rendered offline
// instead of screenshotted from the viewer.
//
// The panel quads are reconstructed exactly the way the viewer draws
// them (LatticeRenderer.cpp): a panel stores its bound quarter-chord
// segment A-B, planform area, and spanwise width; the chordwise
// footprint runs a quarter chord ahead of the bound vortex and three
// quarters aft along the control-point direction.

#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/RotorVaneCoupling.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <fstream>
#include <iostream>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

using namespace Aeolion;

namespace {

constexpr double Rpm = 6000.0;
constexpr double Rho = 1.225;
constexpr double HoverFloorSpeed = 0.01;
constexpr double Omega = Rpm * 2.0 * std::numbers::pi / 60.0;

Geometry::Propeller MakeReferenceProp() {
    Geometry::Propeller prop;
    prop.BladeCount = 2;
    prop.Radius = 0.15;
    prop.HubRadius = 0.02;
    for (int i = 0; i < 12; ++i) {
        const double r = prop.HubRadius + (prop.Radius - prop.HubRadius) * (i + 0.5) / 13.0;
        const double eta = r / prop.Radius;
        prop.Stations.push_back({r, 0.025 * (1.0 - 0.3 * eta), 30.0 - 20.0 * eta});
    }
    return prop;
}

std::vector<Geometry::ControlSurface> MakeVaneSet() {
    std::vector<Geometry::ControlSurface> surfaces;
    const Math::Vec3 axes[] = {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, -1.0}};
    for (const Math::Vec3& axis : axes) {
        Geometry::ControlSurface vane;
        vane.Name = "vane";
        vane.Binding = Geometry::ControlSurfaceBinding::DuctJet;
        vane.ChordFraction = 1.0;
        vane.EtaStart = 0.4;
        vane.EtaEnd = 1.0;
        vane.HingeAxis = axis;
        surfaces.push_back(vane);
    }
    return surfaces;
}

void WriteVec3(std::ofstream& out, const Math::Vec3& v) {
    out << '[' << v.x << ',' << v.y << ',' << v.z << ']';
}

// The viewer's quad reconstruction (see LatticeRenderer.cpp): the drawn
// quad is the panel's planform footprint, so it is built from the
// projected area -- Area is the curved surface and would over-length the
// chord.
void WriteLiftingQuad(std::ofstream& out, const Solver::Panel& panel, double gamma, bool first) {
    const double chord =
        (panel.SpanwiseWidth > 0.0) ? panel.PlanformArea / panel.SpanwiseWidth : 0.0;
    Math::Vec3 chordDir = panel.ControlPoint - (panel.A + panel.B) * 0.5;
    const double len = chordDir.Norm();
    chordDir = (len > 1e-12) ? chordDir * (1.0 / len) : Math::Vec3(1.0, 0.0, 0.0);

    const Math::Vec3 leadingA = panel.A - chordDir * (0.25 * chord);
    const Math::Vec3 leadingB = panel.B - chordDir * (0.25 * chord);
    const Math::Vec3 trailingA = panel.A + chordDir * (0.75 * chord);
    const Math::Vec3 trailingB = panel.B + chordDir * (0.75 * chord);

    if (!first) out << ",\n";
    out << R"(  {"surface":")" << panel.Surface << R"(","gamma":)" << gamma << R"(,"normal":)";
    WriteVec3(out, panel.Normal);
    out << R"(,"corners":[)";
    WriteVec3(out, leadingA);
    out << ',';
    WriteVec3(out, leadingB);
    out << ',';
    WriteVec3(out, trailingB);
    out << ',';
    WriteVec3(out, trailingA);
    out << "]}";
}

void WriteSourceQuad(std::ofstream& out, const Lattice::SourcePanel& panel, bool first) {
    if (!first) out << ",\n";
    out << R"(  {"surface":")" << panel.Surface << R"(","normal":)";
    WriteVec3(out, panel.Normal);
    out << R"(,"corners":[)";
    for (int i = 0; i < 4; ++i) {
        if (i) out << ',';
        WriteVec3(out, panel.Corners[i]);
    }
    out << "]}";
}

} // namespace

int main(int argc, char** argv) {
    const std::string outPath = (argc > 1) ? argv[1] : "lattice-solution.json";

    const auto prop = MakeReferenceProp();
    const double shroudInner = 1.04 * prop.Radius;
    const double shroudOuter = 1.2 * prop.Radius;
    const double shroudChord = 0.5 * prop.Radius;
    const auto panels = PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega);
    const auto strips = PanelBuilder::BuildPropellerStrips(prop);
    const auto duct = PanelBuilder::BuildPropellerDuct(shroudInner, shroudOuter, shroudChord);
    const auto surfaces = MakeVaneSet();

    Solver::FreestreamConditions fc;
    fc.Vinf = HoverFloorSpeed;
    fc.alphaDeg = 0.0;
    fc.betaDeg = 0.0;
    fc.rho = Rho;
    fc.p = Omega;
    fc.RefPoint = Solver::Vec3(0.0, 0.0, 0.0);

    Solver::ReferenceGeometry ref;
    for (const auto& panel : panels) ref.Area += panel.PlanformArea;
    ref.Span = 2.0 * prop.Radius;
    ref.Chord = ref.Area / ref.Span;
    const double trail = Solver::DefaultTrailSpanFactor * 2.0 * prop.Radius;

    Solver::FreestreamConditions vaneFc = fc;
    vaneFc.p = 0.0;
    Solver::ReferenceGeometry vaneRef;
    vaneRef.Area = 1.0;
    vaneRef.Span = 2.0 * shroudInner;
    vaneRef.Chord = shroudChord;

    const std::vector<double> deflections; // neutral cruciform
    const Solver::VaneBuilder vaneBuilder =
        [&](const std::function<Solver::Vec3(const Solver::Vec3&)>& flow)
        -> std::pair<std::vector<Solver::Panel>, std::vector<Solver::StripSection>> {
        return {PanelBuilder::BuildDuctVanes(surfaces, shroudInner, 0.5 * shroudChord, shroudChord,
                                             deflections, flow),
                PanelBuilder::BuildDuctVaneStrips(surfaces, shroudInner, shroudChord, deflections)};
    };
    const Solver::RotorBuilder rotorBuilder = [&](const std::function<double(double)>& inflow) {
        return PanelBuilder::BuildPropellerLattice(prop, 0.0, Omega, inflow);
    };

    const auto result = Solver::SolveRotorVaneCoupled(
        panels, strips, fc, ref, trail, Solver::AnalyticSectionModel{}, duct, vaneBuilder, vaneFc,
        vaneRef, Solver::DefaultTrailSpanFactor * 2.0 * shroudInner, Solver::AnalyticSectionModel{},
        0.0, Rho, {}, rotorBuilder);

    const double thrust = -(result.Rotor.Base.Di + result.Vanes.Base.Di);
    std::cout << "coupled cascade hover: converged=" << result.Converged
              << " outer=" << result.OuterIterations << " residual=" << result.Residual
              << " T=" << thrust << " N\n";
    if (!result.Converged) {
        std::cerr << "refusing to export an unconverged solution\n";
        return 1;
    }

    const auto& bladePanels = result.BladePanels.empty() ? panels : result.BladePanels;
    const auto& bladeGamma = result.Rotor.Base.gamma;
    const auto& vaneGamma = result.Vanes.Base.gamma;

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot open " << outPath << " for writing\n";
        return 1;
    }
    out.precision(9);

    out << "{\n\"meta\":{\"rpm\":" << Rpm << ",\"rho\":" << Rho
        << ",\"bladeCount\":" << prop.BladeCount << ",\"radius\":" << prop.Radius
        << ",\"hubRadius\":" << prop.HubRadius << ",\"shroudInner\":" << shroudInner
        << ",\"shroudOuter\":" << shroudOuter << ",\"shroudChord\":" << shroudChord
        << ",\"thrustN\":" << thrust << ",\"outerIterations\":" << result.OuterIterations
        << "},\n";

    out << "\"bladePanels\":[\n";
    for (std::size_t i = 0; i < bladePanels.size(); ++i)
        WriteLiftingQuad(out, bladePanels[i], (i < bladeGamma.size()) ? bladeGamma[i] : 0.0, i == 0);
    out << "\n],\n\"vanePanels\":[\n";
    for (std::size_t i = 0; i < result.VanePanels.size(); ++i)
        WriteLiftingQuad(out, result.VanePanels[i], (i < vaneGamma.size()) ? vaneGamma[i] : 0.0,
                         i == 0);
    out << "\n],\n\"ductPanels\":[\n";
    for (std::size_t i = 0; i < duct.size(); ++i) WriteSourceQuad(out, duct[i], i == 0);
    out << "\n]\n}\n";

    std::cout << "wrote " << outPath << " (" << bladePanels.size() << " blade, "
              << result.VanePanels.size() << " vane, " << duct.size() << " duct panels)\n";
    return 0;
}
