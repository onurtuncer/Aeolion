// FanInductionExport.cpp -- the aft-fan induction paper's figure data.
//
// Exports three things the paper needs, from the same vortex-cylinder model
// the solve uses (Solver/DiskInduction.h), on the real geometry handoff:
//
//   1. The AXIS profile, discretized against the closed form. This is the
//      verification figure: the model's only exactly-known answer.
//   2. The chordwise induction gradient across the span, at several points
//      along a transition. This is the mechanism figure -- the gradient is
//      what delays separation, and it is what varies across the span.
//   3. The induction field on a meridional slice, for the portrait.
//
// Geometry comes from the contract, converted to solver axes (x aft), with
// the disk annular between the body radius at the duct station and the
// blade tip. Rerun after any change to DiskInduction.h.

#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/DiskInduction.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace Aeolion;
namespace S = Aeolion::Solver;

namespace {

constexpr double Rho = 1.225;

struct Operating {
    const char* Name;
    double Thrust; // [N]
    double Vinf;   // [m/s]
};

// Four points along a transition, from near-hover to cruise. Thrust falls
// as the wing takes over the lift, which is what a real transition does.
const Operating Points[] = {
    {"hover", 25.0, 5.0},
    {"early", 25.0, 12.0},
    {"late", 20.0, 18.0},
    {"cruise", 12.0, 25.0},
};

S::ActuatorDisk MakeDisk(const Geometry::HandoffContract& contract, double thrust, double vinf) {
    S::ActuatorDisk disk;
    const double contractX = contract.Duct.Center.x;
    disk.Center = S::Vec3(-contractX, contract.Duct.Center.y, -contract.Duct.Center.z);
    disk.Axis = S::Vec3(1.0, 0.0, 0.0);
    disk.Radius = contract.Propulsion.DiskRadius;
    disk.HubRadius = Geometry::RadiusAt(contract.Body, contractX);
    if (!(disk.HubRadius < disk.Radius)) disk.HubRadius = 0.0;
    disk.InducedVelocity = S::InducedVelocityFromThrust(thrust, Rho, disk.Area(), vinf);
    return disk;
}

} // namespace

int main(int argc, char** argv) {
    const std::string handoff =
        (argc > 1) ? argv[1] : "tests/Data/AeolionGeometryHandoff-1.8.0.json";
    const std::string outPath = (argc > 2) ? argv[2] : "fan-induction.json";

    Geometry::HandoffContract contract;
    try {
        contract = Geometry::LoadHandoff(handoff);
    } catch (const std::exception& error) {
        std::cerr << "cannot load " << handoff << ": " << error.what() << "\n";
        return 1;
    }

    const double wingLE = -contract.Placement.RootLeadingEdge.x;
    const double chord = contract.Stations.empty() ? 0.1771 : contract.Stations.front().Chord;
    const double semiSpan = 0.5 * contract.Span;

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot open " << outPath << "\n";
        return 1;
    }
    out.precision(9);

    const S::ActuatorDisk reference = MakeDisk(contract, Points[1].Thrust, Points[1].Vinf);
    out << "{\n\"meta\":{\"radius\":" << reference.Radius << ",\"hubRadius\":" << reference.HubRadius
        << ",\"area\":" << reference.Area() << ",\"diskX\":" << reference.Center.x
        << ",\"wingLE\":" << wingLE << ",\"chord\":" << chord << ",\"semiSpan\":" << semiSpan
        << ",\"wingTE\":" << (wingLE + chord) << "},\n";

    // --- 1. axis profile, discretized against the closed form --------------
    {
        S::ActuatorDisk solid = reference;
        solid.HubRadius = 0.0;
        solid.InducedVelocity = 1.0; // report u/vi
        const auto meshSolid = S::BuildVortexCylinder(solid);
        S::ActuatorDisk annulus = reference;
        annulus.InducedVelocity = 1.0;
        const auto meshAnn = S::BuildVortexCylinder(annulus);

        out << "\"axis\":[";
        bool first = true;
        for (double sOverR = -6.0; sOverR <= 10.001; sOverR += 0.05) {
            const double s = sOverR * solid.Radius;
            if (!first) out << ',';
            first = false;
            out << '[' << sOverR << ',' << S::DiskAxisInducedVelocity(solid, s) << ','
                << S::VortexCylinderVelocity(meshSolid, S::Vec3(solid.Center.x + s, 0, 0)).x << ','
                << S::DiskAxisInducedVelocity(annulus, s) << ','
                << S::VortexCylinderVelocity(meshAnn, S::Vec3(annulus.Center.x + s, 0, 0)).x << ']';
        }
        out << "],\n";
    }

    // --- 2. chordwise gradient across the span, per operating point --------
    out << "\"operating\":[\n";
    for (std::size_t k = 0; k < std::size(Points); ++k) {
        const Operating& op = Points[k];
        const S::ActuatorDisk disk = MakeDisk(contract, op.Thrust, op.Vinf);
        const auto mesh = S::BuildVortexCylinder(disk);
        const double vh = std::sqrt(op.Thrust / (2.0 * Rho * disk.Area()));

        if (k) out << ",\n";
        out << " {\"name\":\"" << op.Name << "\",\"thrust\":" << op.Thrust
            << ",\"Vinf\":" << op.Vinf << ",\"vi\":" << disk.InducedVelocity << ",\"vh\":" << vh
            << ",\"mu\":" << op.Vinf / vh << ",\"span\":[";
        bool first = true;
        for (double eta = 0.06; eta <= 1.0001; eta += 0.02) {
            const double y = eta * semiSpan;
            const double uLE = S::VortexCylinderVelocity(mesh, S::Vec3(wingLE, y, 0.0)).x;
            const double uTE = S::VortexCylinderVelocity(mesh, S::Vec3(wingLE + chord, y, 0.0)).x;
            if (!first) out << ',';
            first = false;
            out << '[' << eta << ',' << uLE << ',' << uTE << ',' << (uTE - uLE) / op.Vinf << ']';
        }
        out << "]}";
    }
    out << "\n],\n";

    // --- 3. the coupled configuration itself -------------------------------
    // Exported so the figure can SHOW the two geometric claims the paper
    // rests on -- the fan sits aft of the wing, and its radius covers the
    // span station that separates first -- rather than only asserting them.
    {
        namespace PB = Aeolion::PanelBuilder;
        PB::LatticeOptions options;
        options.BodyCircumferentialPanels = 16;
        Geometry::HandoffContract single = contract;
        single.Mesh.ChordwisePanels = 1;
        PB::LatticeBuilder builder(single, options);
        const auto wing = builder.Build();

        out << "\"configuration\":{\"trimEta\":" << builder.TrimEta() << ",\"wing\":[";
        bool first = true;
        for (const auto& panel : wing) {
            const double c = (panel.SpanwiseWidth > 0.0)
                                 ? panel.PlanformArea / panel.SpanwiseWidth : 0.0;
            S::Vec3 dir = panel.ControlPoint - (panel.A + panel.B) * 0.5;
            const double n = dir.Norm();
            dir = (n > 1e-12) ? dir * (1.0 / n) : S::Vec3(1, 0, 0);
            const S::Vec3 le = (panel.A + panel.B) * 0.5 - dir * (0.25 * c);
            if (!first) out << ',';
            first = false;
            out << '[' << 0.5 * (panel.A.y + panel.B.y) << ',' << le.x << ',' << (le.x + c) << ']';
        }
        out << "],\"body\":[";
        first = true;
        for (const auto& st : contract.Body.Stations) {
            if (!first) out << ',';
            first = false;
            out << '[' << -st.x << ',' << st.Radius << ']';
        }
        out << "],\"duct\":{\"x0\":" << (-contract.Duct.Center.x - 0.5 * contract.Duct.Chord)
            << ",\"x1\":" << (-contract.Duct.Center.x + 0.5 * contract.Duct.Chord)
            << ",\"rInner\":" << 0.5 * contract.Duct.InnerDiameter
            << ",\"rOuter\":" << 0.5 * contract.Duct.OuterDiameter << "}},\n";

        // The same geometry as QUADS, for the three-dimensional view: wing
        // panels reconstructed the way the viewer draws them (a quarter
        // chord ahead of the bound segment, three quarters aft), and the
        // body and duct as the source panels the solve actually carries.
        const auto body3d = builder.BuildBody();
        const auto duct3d = builder.BuildDuct();
        const auto quad = [&out](const S::Vec3 c[4], bool first) {
            if (!first) out << ',';
            out << '[';
            for (int i = 0; i < 4; ++i)
                out << (i ? "," : "") << '[' << c[i].x << ',' << c[i].y << ',' << c[i].z << ']';
            out << ']';
        };

        out << "\"geometry3d\":{\"wing\":[";
        bool firstQuad = true;
        for (const auto& panel : wing) {
            const double c = (panel.SpanwiseWidth > 0.0)
                                 ? panel.PlanformArea / panel.SpanwiseWidth : 0.0;
            S::Vec3 dir = panel.ControlPoint - (panel.A + panel.B) * 0.5;
            const double n = dir.Norm();
            dir = (n > 1e-12) ? dir * (1.0 / n) : S::Vec3(1, 0, 0);
            const S::Vec3 corners[4] = {panel.A - dir * (0.25 * c), panel.B - dir * (0.25 * c),
                                        panel.B + dir * (0.75 * c), panel.A + dir * (0.75 * c)};
            quad(corners, firstQuad);
            firstQuad = false;
        }
        out << "],\"body\":[";
        firstQuad = true;
        for (const auto& panel : body3d) {
            const S::Vec3 corners[4] = {panel.Corners[0], panel.Corners[1], panel.Corners[2],
                                        panel.Corners[3]};
            quad(corners, firstQuad);
            firstQuad = false;
        }
        out << "],\"duct\":[";
        firstQuad = true;
        for (const auto& panel : duct3d) {
            const S::Vec3 corners[4] = {panel.Corners[0], panel.Corners[1], panel.Corners[2],
                                        panel.Corners[3]};
            quad(corners, firstQuad);
            firstQuad = false;
        }
        out << "]},\n";
    }

    // --- 4. meridional slice of the induction field ------------------------
    {
        const S::ActuatorDisk disk = MakeDisk(contract, Points[1].Thrust, Points[1].Vinf);
        const auto mesh = S::BuildVortexCylinder(disk);
        out << "\"field\":{\"vi\":" << disk.InducedVelocity << ",\"points\":[";
        bool first = true;
        for (double x = 0.05; x <= 0.75001; x += 0.01) {
            for (double y = -0.55; y <= 0.55001; y += 0.01) {
                const S::Vec3 v = S::VortexCylinderVelocity(mesh, S::Vec3(x, y, 0.0));
                if (!first) out << ',';
                first = false;
                out << '[' << x << ',' << y << ',' << v.x << ']';
            }
        }
        out << "]}\n}\n";
    }

    std::cout << "wrote " << outPath << " (disk r=" << reference.Radius
              << " hub=" << reference.HubRadius << " at x=" << reference.Center.x << ")\n";
    return 0;
}
