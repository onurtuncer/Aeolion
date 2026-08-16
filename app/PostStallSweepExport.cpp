// PostStallSweepExport.cpp -- Phase 0 of the post-separation study: the
// baseline post-stall map of the coupled wing--body configuration, computed
// with the CURRENT analytic deep-stall blend (AnalyticSectionModel), over an
// attitude matrix that runs far past separation:
//
//     alpha in [-4, 90] deg  x  beta in [0, 30] deg.
//
// What this run is FOR, and what it is not. It establishes the scaffolding
// and the convergence METRICS of the flat-plate-limit study -- where in
// (alpha, beta) the configuration's loads collapse onto bluff-plate
// scaling -- and it shakes out the Level-2 coupling at attitudes far beyond
// anything the existing consumers run. The deep-stall numbers themselves
// are set by the four hand-tuned constants of the analytic blend
// (DeepStallStartDeg/EndDeg, PlateNormal, ClMax), which is exactly what the
// study's later phases replace with anchored models (Viterna's AR-aware
// CDmax, Kirchhoff attenuation from the computed separation point). This
// export is the baseline those phases are compared against, not a result
// to quote on its own -- the same reason TODO.md 3b declined to put
// post-stall coefficients in the second paper.
//
// The solve runs on the CLEAN lattice (no carry-through): a carry strip's
// bound midpoint sits at the centreline INSIDE the fuselage, where the
// source-panel field is the interior continuation -- meaningless as the
// onset flow of a section model. The cost, the root circulation the
// carry-through restores, is quantified per condition by the inviscid CARRY
// solve exported alongside (CLCarryInviscid vs CLCleanInviscid).
//
// Continuation: within each beta column alpha ascends, each solve warm
// started from the previous circulation (ViscousCouplingOptions::
// InitialGamma). The map is therefore the ASCENDING branch; the descending
// branch -- static hysteresis -- is a deliberate later study, not an
// accident of cold starts landing on either side of the fold.
//
// The convergence metrics, per condition:
//
//   #1 forceAngleDeg -- angle between the total force and the chord-plane
//      normal (+z). Attached flow with full leading-edge suction carries
//      the force perpendicular to the FREESTREAM (angle ~ alpha, exported
//      as attachedAngleDeg for the comparison); a stalled plate carries it
//      perpendicular to the CHORD PLANE (angle -> 0). The ratio of the two
//      is the leading-edge suction fraction still alive.
//   #2 sigmaDeg -- total inclination of the freestream to the chord plane,
//      sin(sigma) = sin(alpha) cos(beta): the collapse variable. Bluff-
//      plate loading depends on attitude only through sigma, so post-stall
//      CN(alpha, beta) should merge into one curve CN(sigma).
//   #3 xcp -- centre of pressure from the moment balance: walks from the
//      attached quarter-chord toward the plate's mid-chord.
//   #4 CN at alpha = 90: the finite-plate anchor, to be compared against
//      Viterna's CDmax = 1.11 + 0.018 AR (exported in meta).
//
// Per-strip state ([eta, alpha_eff, cl, cd, residual]) and the coupling
// diagnostics ride along: deep post-stall strips exit as limit-cycle means
// (Converged stays false by design), and that has to be visible in the
// data, not laundered.

#include "Aeolion/Geometry/CstSurface.h"
#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/AttachmentBoundaryLayer.h"
#include "Aeolion/Solver/AttachmentLine.h"
#include "Aeolion/Solver/PostStallSection.h"
#include "Aeolion/Solver/SeparationTables.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace Aeolion;
namespace PB = Aeolion::PanelBuilder;
namespace S = Aeolion::Solver;

namespace {

constexpr double FlightSpeed = 25.0; // [m/s] -- same condition as the paper's sweep
constexpr double Rho = 1.225;
constexpr double TrailSpans = 50.0;
constexpr int BodySectors = 16; // see AttachmentSweepExport.cpp's measured warning

// Attached-flow attitudes step 2 degrees for continuity with the paper's
// table; past 20 the analytic blend is fully in charge and 5-degree steps
// resolve everything it can express.
std::vector<double> BuildAlphaGrid() {
    std::vector<double> alphas;
    for (double a = -4.0; a <= 20.0 + 1e-9; a += 2.0) alphas.push_back(a);
    for (double a = 25.0; a <= 90.0 + 1e-9; a += 5.0) alphas.push_back(a);
    return alphas;
}

// One-sided in beta: the unpowered configuration is symmetric and the
// paper's sweep confirmed the +-beta mirror numerically.
const double Betas[] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0};

// One StripSection per single-row panel, with the CST camber's zero-lift
// angle -- without it the section model would lose the camber the lattice
// geometry carries (see the StripSection comment in ViscousCoupling.h).
//
// The frame is the TRUE CHORD frame, not the panel's own: the Weissinger
// panel's Normal carries the camber slope at the control point, so an
// alpha_eff measured against the panel frame already absorbs the lattice's
// zero-lift shift -- and Alpha0Deg would then subtract the same camber a
// second time. Measured (first run of this driver, camber-tilted frames):
// the coupled zero-lift landed at alpha = -7.8 deg, within 0.3 deg of the
// sum of the lattice's -3.9 and the thin-airfoil -4.2 -- the double count,
// exactly the failure mode the ViscousCoupling.h header warns about. The
// handoff wing is rectangular, unswept and untwisted (established in the
// journal-of-aircraft paper), so its true chord frame is the solver frame
// itself; a swept or twisted wing needs a PanelBuilder-side strip builder
// that carries the section plane, which is where this belongs long-term.
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

void WriteVec3(std::ofstream& out, const Math::Vec3& v) {
    out << '[' << v.x << ',' << v.y << ',' << v.z << ']';
}

} // namespace

int main(int argc, char** argv) {
    const std::string handoffPath =
        (argc > 1) ? argv[1] : "tests/Data/AeolionGeometryHandoff-1.8.0.json";
    const std::string outPath = (argc > 2) ? argv[2] : "poststall-sweep.json";
    const double flightSpeed = (argc > 3) ? std::atof(argv[3]) : FlightSpeed;
    // Coupling knobs, exposed for the stabilization experiments this run
    // exists to perform: the wing lattice's 44 tightly packed strips couple
    // through their shared trailing legs far more stiffly than any propeller
    // consumer's 12, and the first full run showed the accelerated fixed
    // point both limit-cycling and landing on spanwise-checkerboard branches
    // (a known spurious-equilibrium family of collocation lifting-line).
    // Plain heavily-damped iteration suppresses the checkerboard mode at the
    // cost of many more iterations -- measured here, decided on the data.
    const double relaxation = (argc > 4) ? std::atof(argv[4]) : S::DefaultCouplingRelaxation;
    const int andersonDepth = (argc > 5) ? std::atoi(argv[5]) : 4;
    const int maxIterations = (argc > 6) ? std::atoi(argv[6]) : S::DefaultCouplingMaxIterations;
    // argv[7]: a single beta to run, or "all" (the default). argv[8]: the
    // section model -- "analytic" (Phase-0 baseline, hand-tuned deep-stall
    // blend) or "anchored" (Phase 1: Kirchhoff attenuation from the
    // computed separation tables + Viterna deep stall + Rayleigh cm).
    const std::string betaArg = (argc > 7) ? argv[7] : "all";
    const bool betaFiltered = betaArg != "all";
    const double betaOnly = betaFiltered ? std::atof(betaArg.c_str()) : 0.0;
    const std::string modelName = (argc > 8) ? argv[8] : "analytic";
    if (modelName != "analytic" && modelName != "anchored") {
        std::cerr << "unknown model '" << modelName << "' (analytic | anchored)\n";
        return 1;
    }

    Geometry::HandoffContract contract;
    try {
        contract = Geometry::LoadHandoff(handoffPath);
    } catch (const std::exception& error) {
        std::cerr << "cannot load " << handoffPath << ": " << error.what() << "\n";
        return 1;
    }
    contract.Mesh.ChordwisePanels = 1; // the coupling's one-Weissinger-row-per-strip contract

    PB::LatticeOptions cleanOptions;
    cleanOptions.BodyCircumferentialPanels = BodySectors;
    cleanOptions.CarryThroughLift = false; // every strip midpoint outside the body
    PB::LatticeBuilder cleanBuilder(contract, cleanOptions);
    PB::LatticeOptions carryOptions = cleanOptions;
    carryOptions.CarryThroughLift = true; // the physical inviscid reference
    PB::LatticeBuilder carryBuilder(contract, carryOptions);

    const auto wing = cleanBuilder.Build();
    const auto wingCarry = carryBuilder.Build();
    const auto body = cleanBuilder.BuildBody();
    const auto duct = cleanBuilder.BuildDuct();
    std::vector<Lattice::SourcePanel> sources = body;
    sources.insert(sources.end(), duct.begin(), duct.end());

    const double halfSpan = 0.5 * contract.Span;
    const auto strips = StripsFromPanels(wing, halfSpan, contract.AirfoilSections);

    const double trail = TrailSpans * contract.Span;
    const auto preparedClean = S::Prepare(S::PanelSystem{wing, sources}, trail);
    const auto preparedCarry = S::Prepare(S::PanelSystem{wingCarry, sources}, trail);

    S::ReferenceGeometry ref;
    ref.Area = carryBuilder.GrossPlanformArea();
    ref.Span = contract.Span;
    ref.Chord = ref.Area / contract.Span;
    const double aspectRatio = contract.Span * contract.Span / ref.Area;
    const double viternaCDmax = 1.11 + 0.018 * aspectRatio;

    S::FreestreamConditions fc;
    fc.Vinf = flightSpeed;
    fc.rho = Rho;
    // Contract frame is x-forward / z-down, solver x-aft / z-up (see the
    // moment-arm bug note in AttachmentSweepExport.cpp).
    if (contract.MomentReferencePointStated)
        fc.RefPoint = S::Vec3(-contract.MomentReferencePoint.x, contract.MomentReferencePoint.y,
                              -contract.MomentReferencePoint.z);

    // The wing plane's height, for the xcp moment balance: My picks up a
    // (z_wing - z_ref) * Fx arm once the chordwise force grows post-stall.
    const double wingPlaneZ = -contract.Placement.RootLeadingEdge.z;
    const double wingLeadingEdgeX = -contract.Placement.RootLeadingEdge.x;
    const double rootChord = contract.Stations.empty() ? ref.Chord : contract.Stations.front().Chord;

    // The Phase-0 baseline blend stays constructed either way: its
    // DeepStallStartDeg doubles as the diagnostic threshold for counting
    // deep-stall strips, independent of which model actually solves.
    const S::AnalyticSectionModel analytic{};
    S::SectionModel model = analytic;
    S::PostStallSectionModel anchored;
    if (modelName == "anchored") {
        anchored.AspectRatio = aspectRatio;
        auto tables =
            S::BuildSeparationTables(preparedClean, wing, strips, contract.AirfoilSections, fc, ref);
        std::size_t resolved = 0;
        for (const S::StripSeparationTable& table : tables)
            if (!table.AlphaDeg.empty()) ++resolved;
        std::vector<double> etas;
        etas.reserve(strips.size());
        for (const S::StripSection& strip : strips) etas.push_back(strip.Eta);
        anchored.SeparationPoint = S::MakeSeparationFunction(std::move(tables), std::move(etas));
        model = anchored;
        std::cout << "anchored model: separation tables on " << resolved << "/" << strips.size()
                  << " strips; emergent stall (deg from zero lift) at eta 0.15/0.5/0.9 = "
                  << anchored.StallAngleDeg(0.15) << "/" << anchored.StallAngleDeg(0.5) << "/"
                  << anchored.StallAngleDeg(0.9) << "\n";
    }
    const double q = 0.5 * Rho * flightSpeed * flightSpeed;

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot open " << outPath << " for writing\n";
        return 1;
    }
    out.precision(9);

    out << "{\n\"meta\":{\"handoff\":\"" << contract.SchemaVersion << "\",\"Vinf\":" << flightSpeed
        << ",\"rho\":" << Rho << ",\"span\":" << contract.Span << ",\"area\":" << ref.Area
        << ",\"chord\":" << ref.Chord << ",\"aspectRatio\":" << aspectRatio
        << ",\"viternaCDmax\":" << viternaCDmax << ",\"trimEta\":" << cleanBuilder.TrimEta()
        << ",\"wingStrips\":" << wing.size() << ",\"bodyPanels\":" << body.size()
        << ",\"ductPanels\":" << duct.size() << ",\"wingLEx\":" << wingLeadingEdgeX
        << ",\"wingPlaneZ\":" << wingPlaneZ << ",\"rootChord\":" << rootChord << ",\"refPoint\":";
    WriteVec3(out, fc.RefPoint);
    out << ",\"sectionModel\":\"" << modelName << "\",\"model\":{\"ClMax\":" << analytic.ClMax
        << ",\"Cd0\":" << analytic.Cd0 << ",\"KCd\":" << analytic.KCd
        << ",\"DeepStallStartDeg\":" << analytic.DeepStallStartDeg
        << ",\"DeepStallEndDeg\":" << analytic.DeepStallEndDeg
        << ",\"PlateNormal\":" << analytic.PlateNormal << "},\n \"alphas\":[";
    const std::vector<double> alphas = BuildAlphaGrid();
    for (std::size_t i = 0; i < alphas.size(); ++i) out << (i ? "," : "") << alphas[i];
    out << "],\"betas\":[";
    for (std::size_t i = 0; i < std::size(Betas); ++i) out << (i ? "," : "") << Betas[i];
    out << "]},\n\"conditions\":[\n";

    const auto start = std::chrono::steady_clock::now();
    bool firstCondition = true;
    for (const double betaDeg : Betas) {
        if (betaFiltered && std::fabs(betaDeg - betaOnly) > 1e-9) continue;
        S::ViscousCouplingOptions options; // fresh per column: first alpha cold-starts
        options.Relaxation = relaxation;
        options.AndersonDepth = andersonDepth;
        options.MaxIterations = maxIterations;
        for (const double alphaDeg : alphas) {
            fc.alphaDeg = alphaDeg;
            fc.betaDeg = betaDeg;

            const S::ViscousCoupledResult res =
                S::SolveViscousCoupled(wing, strips, fc, ref, trail, model, options, sources);
            options.InitialGamma = res.Base.gamma; // continuation up the alpha column

            // Inviscid references on both lattices: the carry/clean gap
            // prices the omitted carry-through, and the clean one is the
            // "no section model" control for the same lattice.
            const S::SolveResult carryInviscid = S::SolveWithSystem(preparedCarry, fc, ref);
            const S::SolveResult cleanInviscid = S::SolveWithSystem(preparedClean, fc, ref);

            // Total force and moment in solver body axes (x aft, y right,
            // z up), the strips' circulatory + profile forces plus the
            // body's pressure integral.
            S::Vec3 forceWing(0, 0, 0);
            for (const S::StripState& strip : res.Strips) forceWing = forceWing + strip.Force;
            const S::Vec3 force = forceWing + res.SourceForce;
            // Wing-only moment alongside the total: the near-field pressure
            // integral over the closed source bodies is the integral the
            // Trefftz work already caught inventing drag on a coupled
            // configuration (TODO.md 3a), so the plate-convergence metrics
            // must be computable without it.
            const S::Vec3 momentWing = res.InducedMoment + res.ProfileMoment + res.SectionMoment;
            const S::Vec3 moment = momentWing + res.SourceMoment;

            const double qS = q * ref.Area;
            const double CN = force.z / qS;             // chord-plane normal
            const double CC = force.x / qS;             // chordwise, +aft
            const double CSide = force.y / qS;
            const double CL = res.Base.L / qS;
            const double CD = res.Base.Di / qS;         // induced + profile (no CD0 buildup)
            const double CY = res.Base.Y / qS;
            const double Cm = moment.y / (qS * ref.Chord);
            const double Croll = moment.x / (qS * ref.Span);
            const double Cn = moment.z / (qS * ref.Span);

            // Metric #1: total-force angle off the chord-plane normal, and
            // the attached-flow reference (the wind lift direction's own
            // angle off +z -- what the angle reads while suction is full).
            const double forceAngleDeg =
                Math::RadToDeg(std::atan2(std::hypot(force.x, force.y), force.z));
            const S::Vec3 dragDir = S::FreestreamVelocity(fc) * (1.0 / flightSpeed);
            const S::Vec3 liftDir = Cross(dragDir, S::Vec3(0, 1, 0)).Normalized();
            const double attachedAngleDeg =
                Math::RadToDeg(std::acos(std::clamp(liftDir.z, -1.0, 1.0)));

            // Metric #2: the collapse variable.
            const double sinSigma = std::sin(Math::DegToRad(alphaDeg)) * std::cos(Math::DegToRad(betaDeg));
            const double sigmaDeg = Math::RadToDeg(std::asin(std::clamp(sinSigma, -1.0, 1.0)));

            // Metric #3: centre of pressure along x from the pitching
            // moment about RefPoint, forces applied in the wing plane:
            //   My = (z_wing - z_ref) Fx - (x_cp - x_ref) Fz.
            const bool xcpValid = std::fabs(force.z) > 1e-9;
            const double xcp =
                xcpValid
                    ? fc.RefPoint.x + ((wingPlaneZ - fc.RefPoint.z) * force.x - moment.y) / force.z
                    : 0.0;

            // Coupling diagnostics: strips outside the residual contract,
            // and strips the section model already treats as deep-stalled.
            int reversedStrips = 0, deepStallStrips = 0;
            for (std::size_t i = 0; i < res.Strips.size(); ++i) {
                if (std::fabs(res.Strips[i].alphaEffDeg) > options.ResidualIncidenceLimitDeg)
                    ++reversedStrips;
                if (std::fabs(res.Strips[i].alphaEffDeg - strips[i].Alpha0Deg) >
                    analytic.DeepStallStartDeg)
                    ++deepStallStrips;
            }

            std::cout << "alpha=" << alphaDeg << " beta=" << betaDeg << ": CL=" << CL
                      << " CN=" << CN << " CD=" << CD << " xcp=" << (xcpValid ? xcp : 0.0)
                      << " iters=" << res.Iterations << (res.Converged ? "" : " (cycle-mean)")
                      << " res=" << res.MaxResidual << " deepStall=" << deepStallStrips << "/"
                      << res.Strips.size() << "\n";

            if (!firstCondition) out << ",\n";
            firstCondition = false;
            out << R"( {"alphaDeg":)" << alphaDeg << R"(,"betaDeg":)" << betaDeg
                << R"(,"sigmaDeg":)" << sigmaDeg << R"(,"CN":)" << CN << R"(,"CC":)" << CC
                << R"(,"CSide":)" << CSide << R"(,"CL":)" << CL << R"(,"CD":)" << CD
                << R"(,"CY":)" << CY << R"(,"Cm":)" << Cm << R"(,"Croll":)" << Croll
                << R"(,"Cn":)" << Cn << R"(,"forceAngleDeg":)" << forceAngleDeg
                << R"(,"attachedAngleDeg":)" << attachedAngleDeg << R"(,"xcpValid":)"
                << (xcpValid ? "true" : "false") << R"(,"xcp":)" << xcp
                << R"(,"CLCarryInviscid":)" << carryInviscid.CL << R"(,"CLCleanInviscid":)"
                << cleanInviscid.CL << R"(,"force":)";
            WriteVec3(out, force);
            out << R"(,"moment":)";
            WriteVec3(out, moment);
            out << R"(,"forceWing":)";
            WriteVec3(out, forceWing);
            out << R"(,"momentWing":)";
            WriteVec3(out, momentWing);
            out << R"(,"converged":)" << (res.Converged ? "true" : "false")
                << R"(,"iterations":)" << res.Iterations << R"(,"maxResidual":)"
                << res.MaxResidual << R"(,"reversedStrips":)" << reversedStrips
                << R"(,"deepStallStrips":)" << deepStallStrips << ",\n  \"surfaces\":{";
            bool firstSurface = true;
            for (const auto& [name, lift] : res.Base.LiftBySurface) {
                if (!firstSurface) out << ',';
                firstSurface = false;
                const auto dragIt = res.Base.DragBySurface.find(name);
                out << '"' << name << R"(":{"lift":)" << lift << R"(,"drag":)"
                    << ((dragIt != res.Base.DragBySurface.end()) ? dragIt->second : 0.0) << '}';
            }
            out << "},\n  \"strips\":[";
            for (std::size_t i = 0; i < res.Strips.size(); ++i) {
                const double signedEta = ((wing[i].A + wing[i].B) * 0.5).y / halfSpan;
                if (i) out << ',';
                out << '[' << signedEta << ',' << res.Strips[i].alphaEffDeg << ','
                    << res.Strips[i].cl << ',' << res.Strips[i].cd << ','
                    << res.Strips[i].Residual << ',' << res.Strips[i].cm << ']';
            }
            out << "]}";
        }
    }
    out << "\n]\n}\n";

    const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - start);
    std::cout << "wrote " << outPath << " (" << alphas.size() * std::size(Betas)
              << " conditions in " << elapsed.count() << " s)\n";
    return 0;
}
