// AeroMapExport.cpp -- the aero* half of the DAVE-ML flight model
// (models/README.md): the airframe's body-axis force and moment
// coefficients over the alpha x beta grid, plus the reduced-rate
// stability derivatives, exported for the assembler.
//
// TWO SOLVE PATHS, deliberately, because the model's two table families
// want different things:
//
//   BASELINE tables come from the Level-2 coupled solve with the ANCHORED
//   post-stall section model -- Kirchhoff attenuation driven by the
//   computed separation tables, Viterna-Corrigan deep stall, Rayleigh
//   section cm. This is the model's single source of truth across the
//   whole alpha range (models/README.md decision 7): it reduces to the
//   attached solution below separation by construction, so there is no
//   seam to justify in the middle of the flight-critical region.
//
//   RATE DERIVATIVES come from central differences on the INVISCID
//   prepared system, the path Part I validated and TestSolverCore pins
//   (Cl_p = -0.452 against the textbook -0.45 for this AR = 6 wing). They
//   are meaningful only while the flow is attached; the assembler applies
//   the declared taper past 20 degrees. Differencing the coupled solve
//   instead would difference limit-cycle means, which is not a derivative.
//
// BODY AXES. Everything exported is contract-frame FRD (x forward, y
// right, z down), the frame of the DAVE-ML file, converted from the
// solver's x-aft/z-up frame by the single 180-degree-about-y rotation at
// one site. Moments are about the contract's moment_reference_point,
// converted at ingest. Those are the only two conversion sites, which is
// the discipline the 2026-08-09 moment-arm bug earned.
//
// WHAT IS NOT HERE. Aileron increment tables (aeroDC*) are BLOCKED and
// deliberately not emitted: PanelBuilder's MinRowsToResolveHinge = 2
// collides with SolveViscousCoupled's one-row-per-strip contract, so a
// deflection through this path is exactly zero at every attitude --
// measured in app/AileronEffectivenessExport.cpp. Emitting them would
// ship a flight model with no roll control. Parasite drag is its own
// driver (aeolion_parasite_drag) since no solve produces it.
//
// Beta is swept one-sided; the assembler mirrors it using the unpowered
// airframe's symmetry (CY, Cl, Cn odd in beta; CX, CZ, Cm even), which
// the source sweeps confirmed numerically.
//
// Usage:
//   aeolion_aero_map <handoff.json> <out.json> [Vinf] [relaxation]
//                    [maxIterations] [beta|all]

#include "Aeolion/Geometry/CstSurface.h"
#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/BodyAxes.h"
#include "Aeolion/Solver/PostStallSection.h"
#include "Aeolion/Solver/SeparationTables.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace Aeolion;
namespace PB = Aeolion::PanelBuilder;
namespace S = Aeolion::Solver;

namespace {

constexpr double FlightSpeed = 25.0; // the reference condition of every sweep
constexpr double Rho = 1.225;
constexpr double TrailSpans = 50.0;
constexpr int BodySectors = 16; // the measured value; see AttachmentSweepExport.cpp

// The model's alphaBp: 2 degrees through the stall break, coarsening in
// the plate regime where the loading varies slowly with attitude.
std::vector<double> BuildAlphaGrid() {
    std::vector<double> alphas;
    for (double a = -4.0; a <= 26.0 + 1e-9; a += 2.0) alphas.push_back(a);
    for (const double a : {30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0})
        alphas.push_back(a);
    return alphas;
}

// One-sided; the assembler mirrors.
const double Betas[] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0};

// Rate derivatives are computed only where they mean something. Past the
// taper's start the assembler zeroes them anyway.
constexpr double RateDerivativeMaxAlphaDeg = 20.0;

// The TRUE CHORD frame, not the cambered panel axes -- the camber
// double-count trap documented in ViscousCoupling.h and measured twice in
// this repo (Phase 0 finding 2, and the same bug in AttachmentSweepExport).
// The handoff wing is rectangular, unswept and untwisted, so its chord
// frame is the solver frame itself; a swept or twisted wing needs a
// PanelBuilder-side strip builder that carries the section plane.
std::vector<S::StripSection> StripsFromPanels(const std::vector<S::Panel>& panels, double halfSpan,
                                              const std::vector<Geometry::AirfoilSection>& sections) {
    std::vector<S::StripSection> strips;
    strips.reserve(panels.size());
    for (const S::Panel& panel : panels) {
        const S::Vec3 mid = (panel.A + panel.B) * 0.5;
        S::StripSection strip;
        strip.ChordDir = S::Vec3(1.0, 0.0, 0.0);
        strip.LiftDir = S::Vec3(0.0, 0.0, 1.0);
        strip.Chord = (panel.SpanwiseWidth > 0.0) ? panel.PlanformArea / panel.SpanwiseWidth : 0.0;
        strip.Width = panel.SpanwiseWidth;
        strip.Eta = std::fabs(mid.y) / halfSpan;
        strip.Alpha0Deg = Geometry::SectionZeroLiftAngleDeg(sections, strip.Eta);
        strips.push_back(strip);
    }
    return strips;
}

void WriteVec3(std::ofstream& out, const S::Vec3& v) {
    out << '[' << v.x << ',' << v.y << ',' << v.z << ']';
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: aeolion_aero_map <handoff.json> <out.json> [Vinf]"
                     " [relaxation] [maxIterations] [beta|all]\n";
        return 1;
    }
    const std::string handoffPath = argv[1];
    const std::string outPath = argv[2];
    const double flightSpeed = (argc > 3) ? std::atof(argv[3]) : FlightSpeed;
    // Plain damped iteration: the Anderson-accelerated fixed point lands
    // on spanwise-checkerboard equilibria on this many tightly packed
    // strips (the measured Phase-0 finding, TODO.md 3c).
    const double relaxation = (argc > 4) ? std::atof(argv[4]) : 0.05;
    const int maxIterations = (argc > 5) ? std::atoi(argv[5]) : 1000;
    const std::string betaArg = (argc > 6) ? argv[6] : "all";
    const bool betaFiltered = betaArg != "all";
    const double betaOnly = betaFiltered ? std::atof(betaArg.c_str()) : 0.0;

    Geometry::HandoffContract contract;
    try {
        contract = Geometry::LoadHandoff(handoffPath);
    } catch (const std::exception& error) {
        std::cerr << "cannot load " << handoffPath << ": " << error.what() << "\n";
        return 1;
    }
    contract.Mesh.ChordwisePanels = 1; // the coupling's one-row-per-strip contract

    PB::LatticeOptions options;
    options.BodyCircumferentialPanels = BodySectors;
    options.CarryThroughLift = false; // every strip midpoint outside the body
    PB::LatticeBuilder builder(contract, options);

    const auto wing = builder.Build();
    const auto body = builder.BuildBody();
    const auto duct = builder.BuildDuct();
    std::vector<Lattice::SourcePanel> sources = body;
    sources.insert(sources.end(), duct.begin(), duct.end());

    const double halfSpan = 0.5 * contract.Span;
    const auto strips = StripsFromPanels(wing, halfSpan, contract.AirfoilSections);

    const double trail = TrailSpans * contract.Span;
    const auto prepared = S::Prepare(S::PanelSystem{wing, sources}, trail);

    S::ReferenceGeometry ref;
    ref.Area = builder.GrossPlanformArea();
    ref.Span = contract.Span;
    ref.Chord = ref.Area / contract.Span;
    const double aspectRatio = contract.Span * contract.Span / ref.Area;

    S::FreestreamConditions fc;
    fc.Vinf = flightSpeed;
    fc.rho = Rho;
    // Ingest site: the contract states its frame once and the consumer
    // converts explicitly (ADR-0016).
    if (contract.MomentReferencePointStated)
        fc.RefPoint = S::Vec3(-contract.MomentReferencePoint.x, contract.MomentReferencePoint.y,
                              -contract.MomentReferencePoint.z);

    // The anchored section model: Kirchhoff attenuation driven by the
    // computed separation tables.
    S::PostStallSectionModel anchored;
    anchored.AspectRatio = aspectRatio;
    {
        auto tables =
            S::BuildSeparationTables(prepared, wing, strips, contract.AirfoilSections, fc, ref);
        std::size_t resolved = 0;
        for (const S::StripSeparationTable& t : tables)
            if (!t.AlphaDeg.empty()) ++resolved;
        std::vector<double> etas;
        etas.reserve(strips.size());
        for (const S::StripSection& s : strips) etas.push_back(s.Eta);
        anchored.SeparationPoint = S::MakeSeparationFunction(std::move(tables), std::move(etas));
        std::cout << "separation tables on " << resolved << "/" << strips.size()
                  << " strips; emergent stall (deg from zero lift) at eta 0.15/0.5/0.9 = "
                  << anchored.StallAngleDeg(0.15) << "/" << anchored.StallAngleDeg(0.5) << "/"
                  << anchored.StallAngleDeg(0.9) << "\n";
    }
    const S::SectionModel model = anchored;

    S::ViscousCouplingOptions coupling;
    coupling.Relaxation = relaxation;
    coupling.AndersonDepth = 0;
    coupling.MaxIterations = maxIterations;

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot open " << outPath << " for writing\n";
        return 1;
    }
    out.precision(9);

    out << "{\n\"meta\":{\"designId\":\"" << contract.DesignId << "\",\"schema\":\""
        << contract.SchemaVersion << "\",\"Vinf\":" << flightSpeed << ",\"rho\":" << Rho
        << ",\"area\":" << ref.Area << ",\"span\":" << ref.Span << ",\"chord\":" << ref.Chord
        << ",\"aspectRatio\":" << aspectRatio << ",\"wingStrips\":" << wing.size()
        << ",\"bodyPanels\":" << body.size() << ",\"ductPanels\":" << duct.size()
        << ",\"sectionModel\":\"anchored\",\"relaxation\":" << relaxation
        << ",\"frames\":\"body FRD; moments about moment_reference_point\",\"refPoint\":";
    WriteVec3(out, fc.RefPoint);
    out << ",\"rateDerivativeMaxAlphaDeg\":" << RateDerivativeMaxAlphaDeg
        << ",\"excludes\":\"aileron increments (BLOCKED: hinge unrepresentable at one"
           " chordwise row -- see AileronEffectivenessExport.cpp); parasite drag"
           " (aeolion_parasite_drag)\"},\n";

    // --- rate derivatives, inviscid, beta = 0 --------------------------------
    out << "\"rates\":[\n";
    bool firstRate = true;
    for (const double alphaDeg : BuildAlphaGrid()) {
        if (alphaDeg > RateDerivativeMaxAlphaDeg + 1e-9) break;
        S::FreestreamConditions base = fc;
        base.alphaDeg = alphaDeg;
        base.betaDeg = 0.0;
        const S::BodyAxisRateDerivatives d =
            S::ComputeBodyAxisRateDerivatives(prepared, base, ref);
        if (!firstRate) out << ",\n";
        firstRate = false;
        out << R"( {"alphaDeg":)" << alphaDeg << R"(,"CZq":)" << d.CZq << R"(,"Cmq":)" << d.Cmq
            << R"(,"Clp":)" << d.Clp << R"(,"Cnp":)" << d.Cnp << R"(,"CYp":)" << d.CYp
            << R"(,"Clr":)" << d.Clr << R"(,"Cnr":)" << d.Cnr << R"(,"CYr":)" << d.CYr << '}';
        out.flush();
    }
    out << "\n],\n\"baseline\":[\n";

    std::cout << "\nalpha  beta      CX        CZ        Cm        Cl        Cn   iters\n";

    bool firstRow = true;
    int unconverged = 0;
    for (const double betaDeg : Betas) {
        if (betaFiltered && std::fabs(betaDeg - betaOnly) > 1e-9) continue;
        // Warm-start continuation UP each alpha column: the map is the
        // ASCENDING branch, deliberately (hysteresis is a separate study,
        // not an accident of cold starts landing either side of the fold).
        std::vector<double> warmStart;
        for (const double alphaDeg : BuildAlphaGrid()) {
            fc.alphaDeg = alphaDeg;
            fc.betaDeg = betaDeg;
            S::ViscousCouplingOptions opts = coupling;
            opts.InitialGamma = warmStart;

            const auto start = std::chrono::steady_clock::now();
            const auto res =
                S::SolveViscousCoupled(wing, strips, fc, ref, trail, model, opts, sources);
            const double seconds =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            warmStart = res.Base.gamma;

            const double q = 0.5 * Rho * flightSpeed * flightSpeed;
            const S::BodyAxisCoefficients w = S::BodyAxisFromCoupled(res, q, ref);
            if (!res.Converged) ++unconverged;

            if (!firstRow) out << ",\n";
            firstRow = false;
            out << R"( {"alphaDeg":)" << alphaDeg << R"(,"betaDeg":)" << betaDeg << R"(,"CX":)"
                << w.CX << R"(,"CY":)" << w.CY << R"(,"CZ":)" << w.CZ << R"(,"Cl":)" << w.Cl
                << R"(,"Cm":)" << w.Cm << R"(,"Cn":)" << w.Cn << R"(,"converged":)"
                << (res.Converged ? "true" : "false") << R"(,"iterations":)" << res.Iterations
                << R"(,"residual":)" << res.MaxResidual << R"(,"seconds":)" << seconds << '}';
            out.flush();

            std::cout << alphaDeg << "\t" << betaDeg << "\t" << w.CX << "\t" << w.CZ << "\t"
                      << w.Cm << "\t" << w.Cl << "\t" << w.Cn << "\t" << res.Iterations
                      << (res.Converged ? "" : "  (cycle-mean)") << std::endl;
        }
    }

    out << "\n]}\n";
    std::cout << "\nwrote " << outPath << " (" << unconverged
              << " conditions exported as cycle means)\n";
    return 0;
}
