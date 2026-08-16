// ParasiteDragExport.cpp -- the aeroCD0 table of the DAVE-ML flight model
// (models/README.md): the parasite drag a potential-flow method cannot
// produce, over the model's angle-of-attack grid.
//
// WHAT THIS DOES AND DOES NOT COVER -- the whole point of the driver, and
// the thing to get wrong:
//
// The Level-2 coupled solve already carries TWO of the three drag
// contributions. Induced drag comes from the circulation, and the wing's
// PROFILE drag comes from the section model's cd, integrated per strip
// and folded into the same force vector (see the comment at
// PostStallSweepExport.cpp's CD, "induced + profile"). What no
// potential-flow solve can produce is the parasite drag of the closed
// bodies: a source-panelled body carries zero net force in uniform flow
// to machine precision, which is d'Alembert behaving correctly.
//
// So this driver computes the BODY and DUCT contribution ONLY. It must
// never be given the wing: adding a conventional whole-airframe buildup
// on top of the coupled solve's tables double-counts the wing's skin
// friction. That is the single most likely misuse and it would be
// invisible in the totals.
//
// TWO REGIMES, because the model's alpha grid runs to 90 degrees:
//
//   1. Attached branch -- the classic component buildup (Raymer/Hoerner,
//      DragEstimate/DragEstimate.h): skin friction at the component
//      Reynolds number, times a form factor, times an interference
//      factor, per wetted area. Weakly attitude-dependent, and treated
//      here as constant in alpha.
//
//   2. Crossflow branch -- at incidence a body of revolution sheds a
//      crossflow wake whose drag has nothing to do with skin friction.
//      The Allen-Perkins / Jorgensen slender-body form is used:
//
//          CD_cross(alpha) = eta * Cdc * (Splan / Sref) * sin^3(alpha)
//
//      with Splan the body's SIDE-projected area, Cdc the circular
//      cylinder's crossflow drag coefficient, and eta the finite-length
//      proportionality factor. A constant CD0 -- the conventional choice
//      for a cruise-envelope model -- would omit this entirely, and it
//      is the dominant parasite term at high incidence: on this airframe
//      the crossflow term at 90 degrees is several times the friction
//      term.
//
// NOT covered, and stated rather than hidden: the duct's own separated
// drag at incidence (an annular ring at 90 degrees is a bluff body, but
// the slender-body crossflow form does not apply to it), and any
// interference between the body wake and the wing. Both make this an
// UNDER-estimate at high alpha.
//
// Usage:
//   aeolion_parasite_drag <handoff.json> <out.json> [Vinf] [laminarFraction]

#include "Aeolion/DragEstimate/DragEstimate.h"
#include "Aeolion/Geometry/HandoffContract.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using namespace Aeolion;
namespace DE = Aeolion::DragEstimate;

namespace {

constexpr double FlightSpeed = 25.0; // the reference condition of every sweep
constexpr double Rho = 1.225;

// --- crossflow constants ----------------------------------------------------
// Circular cylinder crossflow drag at subcritical Reynolds number. The
// crossflow Reynolds number here is Vinf*sin(alpha)*d/nu ~ 1.7e5 at 90
// degrees, comfortably below the drag crisis, so the subcritical value
// stands across the whole grid.
constexpr double CylinderCrossflowCd = 1.2;
// Jorgensen's finite-length proportionality factor: the fraction of the
// infinite-cylinder crossflow drag a body of this fineness actually
// develops. ~0.65 at fineness 5; it tends to 1 only for very slender
// bodies. Stated rather than tuned.
constexpr double CrossflowEta = 0.65;

// The model's own alpha breakpoints (models/README.md alphaBp), so the
// exported table lands on the grid the assembler expects without
// resampling.
std::vector<double> BuildAlphaGrid() {
    std::vector<double> alphas;
    for (double a = -4.0; a <= 26.0 + 1e-9; a += 2.0) alphas.push_back(a);
    for (const double a : {30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0})
        alphas.push_back(a);
    return alphas;
}

/** Wetted and side-projected area of the body of revolution, by frusta. */
struct BodyAreas {
    double Swet = 0.0;    ///< [m^2] surface of revolution
    double Splan = 0.0;   ///< [m^2] side-projected (crossflow reference)
    double MaxRadius = 0.0;
};

BodyAreas RevolveBody(const Geometry::BodyGeometry& body) {
    BodyAreas areas;
    const auto& st = body.Stations;
    for (std::size_t i = 0; i + 1 < st.size(); ++i) {
        const double dx = std::fabs(st[i].x - st[i + 1].x);
        const double r1 = st[i].Radius, r2 = st[i + 1].Radius;
        // Frustum lateral area pi*(r1+r2)*slant, and the side-view area
        // of the same frustum, (r1+r2)*dx -- the integral of the diameter.
        areas.Swet += std::numbers::pi * (r1 + r2) * std::hypot(dx, r2 - r1);
        areas.Splan += (r1 + r2) * dx;
        areas.MaxRadius = std::max({areas.MaxRadius, r1, r2});
    }
    return areas;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: aeolion_parasite_drag <handoff.json> <out.json>"
                     " [Vinf] [laminarFraction]\n";
        return 1;
    }
    const std::string handoffPath = argv[1];
    const std::string outPath = argv[2];
    const double Vinf = (argc > 3) ? std::atof(argv[3]) : FlightSpeed;
    const double laminarFraction = (argc > 4) ? std::atof(argv[4]) : 0.0;

    Geometry::HandoffContract contract;
    try {
        contract = Geometry::LoadHandoff(handoffPath);
    } catch (const std::exception& error) {
        std::cerr << "cannot load " << handoffPath << ": " << error.what() << "\n";
        return 1;
    }
    if (!contract.Body.IsPresent()) {
        std::cerr << "contract states no body; nothing to estimate\n";
        return 1;
    }

    // Reference area: the wing's gross planform, the same Sref every
    // airframe coefficient in the model is normalised by. Rectangular
    // unswept planform, so span times root chord.
    const double rootChord =
        contract.Stations.empty() ? 0.0 : contract.Stations.front().Chord;
    const double Sref = contract.Span * rootChord;
    if (!(Sref > 0.0)) {
        std::cerr << "degenerate reference area\n";
        return 1;
    }

    const BodyAreas body = RevolveBody(contract.Body);
    const double bodyLength = contract.Body.Length;
    const double fineness = (body.MaxRadius > 0.0) ? bodyLength / (2.0 * body.MaxRadius) : 0.0;

    std::vector<DE::ComponentSpec> specs;

    DE::ComponentSpec bodySpec;
    bodySpec.Name = "body";
    bodySpec.Swet = body.Swet;
    bodySpec.CharLength = bodyLength;
    bodySpec.IsBody = true;
    bodySpec.Fineness = fineness;
    bodySpec.Q = 1.0; // the wing-body junction's interference is the wing's to carry
    bodySpec.LaminarFraction = laminarFraction;
    specs.push_back(bodySpec);

    // The duct ring: two cylindrical walls of the duct's own chord. Its
    // cross-section is a thin annular airfoil, so the airfoil form factor
    // applies with t/c from the wall thickness.
    double ductSwet = 0.0, ductTc = 0.0;
    if (contract.Duct.IsStated) {
        const double Do = contract.Duct.OuterDiameter, Di = contract.Duct.InnerDiameter;
        const double chord = contract.Duct.Chord;
        ductSwet = std::numbers::pi * (Do + Di) * chord;
        ductTc = (chord > 0.0) ? (0.5 * (Do - Di)) / chord : 0.0;

        DE::ComponentSpec duct;
        duct.Name = "duct";
        duct.Swet = ductSwet;
        duct.CharLength = chord;
        duct.IsBody = false;
        duct.ThicknessRatio = ductTc;
        duct.xcMaxThickness = DE::DefaultMaxThicknessLoc;
        duct.SweepDeg = 0.0;
        duct.Q = 1.0;
        duct.LaminarFraction = laminarFraction;
        specs.push_back(duct);
    }

    const DE::AirProperties air{Rho, DE::SeaLevelViscosity};
    const DE::BuildupResult buildup = DE::EstimateCD0(specs, Vinf, Sref, air);

    // The crossflow coefficient, referred to the same Sref.
    const double crossflowCoeff = CrossflowEta * CylinderCrossflowCd * (body.Splan / Sref);

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot open " << outPath << " for writing\n";
        return 1;
    }
    out.precision(9);

    out << "{\n\"meta\":{\"designId\":\"" << contract.DesignId << "\",\"schema\":\""
        << contract.SchemaVersion << "\",\"Vinf\":" << Vinf << ",\"rho\":" << Rho
        << ",\"Sref\":" << Sref << ",\"laminarFraction\":" << laminarFraction
        << ",\"bodyLength\":" << bodyLength << ",\"bodySwet\":" << body.Swet
        << ",\"bodySplan\":" << body.Splan << ",\"bodyFineness\":" << fineness
        << ",\"ductSwet\":" << ductSwet << ",\"ductThicknessRatio\":" << ductTc
        << ",\"crossflowEta\":" << CrossflowEta << ",\"crossflowCdc\":" << CylinderCrossflowCd
        << ",\"crossflowCoeff\":" << crossflowCoeff << ",\"CD0friction\":" << buildup.CD0
        << ",\"miscFraction\":" << buildup.MiscFraction
        << ",\"excludes\":\"wing (its profile drag is already inside the coupled"
           " solve's force tables); duct separated drag at incidence; body-wake"
           " interference\"},\n";

    out << "\"components\":[\n";
    for (std::size_t i = 0; i < buildup.Components.size(); ++i) {
        const DE::ComponentResult& c = buildup.Components[i];
        if (i) out << ",\n";
        out << R"( {"name":")" << c.Name << R"(","Re":)" << c.Re << ",\"Cf\":" << c.Cf
            << ",\"FF\":" << c.FF << ",\"Q\":" << c.Q << ",\"Swet\":" << c.Swet
            << ",\"CD0\":" << c.CD0_contribution << '}';
    }
    out << "\n],\n\"table\":[\n";

    const std::vector<double> alphas = BuildAlphaGrid();
    for (std::size_t i = 0; i < alphas.size(); ++i) {
        const double alphaRad = alphas[i] * std::numbers::pi / 180.0;
        const double s = std::sin(alphaRad);
        // sin^3, and |sin| so the branch is even about zero incidence:
        // crossflow drag does not know the sign of alpha.
        const double crossflow = crossflowCoeff * std::fabs(s * s * s);
        if (i) out << ",\n";
        out << R"( {"alphaDeg":)" << alphas[i] << R"(,"CD0friction":)" << buildup.CD0
            << R"(,"CD0crossflow":)" << crossflow << R"(,"CD0":)" << (buildup.CD0 + crossflow)
            << '}';
    }
    out << "\n]}\n";

    std::cout << "Sref=" << Sref << " m^2  body: Swet=" << body.Swet << " Splan=" << body.Splan
              << " fineness=" << fineness << "\n";
    for (const DE::ComponentResult& c : buildup.Components)
        std::cout << "  " << c.Name << ": Re=" << c.Re << " Cf=" << c.Cf << " FF=" << c.FF
                  << " CD0=" << c.CD0_contribution << "\n";
    std::cout << "CD0 friction (incl. " << 100.0 * buildup.MiscFraction << "% misc) = "
              << buildup.CD0 << "\n"
              << "crossflow coefficient = " << crossflowCoeff << " -> CD0(90 deg) = "
              << (buildup.CD0 + crossflowCoeff) << " ("
              << crossflowCoeff / std::max(buildup.CD0, 1e-12) << "x the friction term)\n"
              << "wrote " << outPath << "\n";
    return 0;
}
