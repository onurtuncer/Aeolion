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
//   AILERON INCREMENTS come from the same coupled solve with the flap
//   carried in the SECTION (StripSection::FlapChordFraction), which is
//   what makes them computable at all: a hinge cannot be represented on
//   a single chordwise row, and before that fix a deflection through
//   this path was exactly zero at every attitude, silently. Neutral and
//   deflected are solved in the SAME pass, each warm-started up its own
//   alpha column, because differencing two limit-cycle means taken from
//   different continuation histories would add noise to the increment
//   that has nothing to do with the control surface.
//
// Parasite drag is NOT here -- it is its own driver
// (aeolion_parasite_drag), since no solve produces it.
//
// TWO ONE-SIDED SWEEPS, both mirrored by the assembler on symmetry
// arguments that this driver also checks numerically rather than
// assuming. Beta: CY, Cl, Cn odd, CX, CZ, Cm even. Aileron: the
// configuration is mirror-symmetric about xz and mirroring maps +delta
// onto -delta, so the same odd/even split holds in deflection; the sweep
// solves one negative deflection to verify it.
//
// Usage:
//   aeolion_aero_map <handoff.json> <out.json> [Vinf] [relaxation]
//                    [maxIterations] [beta|all] [all|baseline|aileron]

#include "Aeolion/Geometry/CstSurface.h"
#include "Aeolion/Geometry/FlapEffectiveness.h"
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
                                              const std::vector<Geometry::AirfoilSection>& sections,
                                              const Geometry::ControlSurface* aileron = nullptr,
                                              double aileronDeg = 0.0) {
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

        // The aileron, carried in the SECTION rather than the panel
        // geometry: a hinge cannot be represented on a single chordwise
        // row, but to thin-airfoil theory a deflected flap is a camber
        // change, so it shifts the zero-lift angle the section model is
        // posed against (StripSection::FlapChordFraction). ANTISYMMETRIC:
        // the right semi-span takes +delta and the left -delta, which is
        // what makes it an aileron rather than a flaperon.
        if (aileron && aileronDeg != 0.0 && strip.Eta >= aileron->EtaStart &&
            strip.Eta <= aileron->EtaEnd) {
            strip.FlapChordFraction = aileron->ChordFraction;
            strip.FlapDeflectionDeg = (mid.y >= 0.0) ? aileronDeg : -aileronDeg;
        }
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
    // Which blocks to compute. The baseline map is expensive and rarely
    // needs regenerating alongside the aileron sweep, so they are
    // selectable: "all" | "baseline" | "aileron".
    const std::string blocks = (argc > 7) ? argv[7] : "all";
    const bool wantBaseline = (blocks == "all" || blocks == "baseline");
    const bool wantAileron = (blocks == "all" || blocks == "aileron");
    if (!wantBaseline && !wantAileron) {
        std::cerr << "unknown block selector '" << blocks << "' (all | baseline | aileron)\n";
        return 1;
    }

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

    // --- rate derivatives ----------------------------------------------------
    // TWO SETS, because neither alone covers the envelope honestly.
    //
    // The INVISCID set (central differences on the prepared system) is the
    // path Part I validated and TestBodyAxes pins against the textbook
    // Cl_p = -0.45. It is exact where the flow is attached and meaningless
    // past separation, since the lattice contains no stall.
    //
    // The COUPLED set differences the Level-2 solve, so the section model's
    // post-stall lift slope enters. That matters more than it sounds: past
    // stall dcl/dalpha goes NEGATIVE, which can drive Cl_p POSITIVE -- roll
    // ANTI-damping, i.e. autorotation, the mechanism of a spin. A table that
    // tapered the attached value to zero would miss that; one that clamped
    // at its last attached value (which is what an alphaRateBp ending at 20
    // silently does) would grant a simulator full attached roll damping at
    // 90 degrees. Both are worse than measuring.
    //
    // Post-stall the coupled solve returns a limit-cycle mean, so a
    // derivative differenced across it is only meaningful if the signal
    // exceeds the cycle's own width. The rate step is therefore enlarged
    // past stall, and each row carries the cycle fluctuation of its own
    // perturbed solves so the assembler can tell a resolved derivative from
    // one lost in the cycle.
    out << "\"rates\":[\n";
    bool firstRate = true;
    for (const double alphaDeg : BuildAlphaGrid()) {
        S::FreestreamConditions base = fc;
        base.alphaDeg = alphaDeg;
        base.betaDeg = 0.0;

        const bool attached = alphaDeg <= RateDerivativeMaxAlphaDeg + 1e-9;
        S::BodyAxisRateDerivatives inv;
        if (attached) inv = S::ComputeBodyAxisRateDerivatives(prepared, base, ref);

        // Coupled roll damping, the derivative whose SIGN carries the
        // physics. A bigger step past stall so the difference clears the
        // cycle width; the attached range keeps the small step so the two
        // sets are comparable there.
        const double step = attached ? 0.05 : 0.40;
        const double q = 0.5 * Rho * flightSpeed * flightSpeed;
        const double reduce = ref.Span / (2.0 * flightSpeed);

        const auto solveAtRate = [&](double p) {
            S::FreestreamConditions fcp = base;
            fcp.p = p;
            S::ViscousCouplingOptions opts = coupling;
            return S::SolveViscousCoupled(wing, strips, fcp, ref, trail, model, opts, sources);
        };
        const auto rp = solveAtRate(+step);
        const auto rm = solveAtRate(-step);
        const S::BodyAxisCoefficients wp = S::BodyAxisFromCoupled(rp, q, ref);
        const S::BodyAxisCoefficients wm = S::BodyAxisFromCoupled(rm, q, ref);
        // FRAME, and the trap this very measurement walked into. Cl and p
        // BOTH flip under the solver->contract rotation, so the DERIVATIVE
        // is invariant -- but only if both sides are in the same frame.
        // Here the moment is already FRD (BodyAxisFromCoupled) while the
        // rate was set on FreestreamConditions in the SOLVER frame, so
        // exactly one flip is outstanding and the quotient needs negating.
        // Caught because the attached range must reproduce the inviscid
        // Cl_p = -0.45 and instead read +0.54: right magnitude, wrong sign,
        // which is the same failure Solver/BodyAxes.h was written about.
        const double clpCoupled = -(wp.Cl - wm.Cl) / (2.0 * step * reduce);
        const double cnpCoupled = -(wp.Cn - wm.Cn) / (2.0 * step * reduce);
        const double fluct = std::max(rp.CycleFluctuation(), rm.CycleFluctuation());

        if (!firstRate) out << ",\n";
        firstRate = false;
        out << R"( {"alphaDeg":)" << alphaDeg << R"(,"attached":)" << (attached ? "true" : "false")
            << R"(,"CZq":)" << inv.CZq << R"(,"Cmq":)" << inv.Cmq << R"(,"Clp":)" << inv.Clp
            << R"(,"Cnp":)" << inv.Cnp << R"(,"CYp":)" << inv.CYp << R"(,"Clr":)" << inv.Clr
            << R"(,"Cnr":)" << inv.Cnr << R"(,"CYr":)" << inv.CYr
            << R"(,"ClpCoupled":)" << clpCoupled << R"(,"CnpCoupled":)" << cnpCoupled
            << R"(,"rateStep":)" << step << R"(,"cycleFluctuation":)" << fluct
            << R"(,"converged":)" << ((rp.Converged && rm.Converged) ? "true" : "false") << '}';
        out.flush();
        std::cout << "rates alpha=" << alphaDeg << "  Clp_inv=" << (attached ? inv.Clp : 0.0)
                  << "  Clp_coupled=" << clpCoupled << "  fluct=" << fluct
                  << ((rp.Converged && rm.Converged) ? "" : "  (cycle-mean)") << std::endl;
    }
    out << "\n],\n\"baseline\":[\n";

    std::cout << "\nalpha  beta      CX        CZ        Cm        Cl        Cn   iters\n";

    bool firstRow = true;
    int unconverged = 0;
    for (const double betaDeg : Betas) {
        if (!wantBaseline) break;
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

    out << "\n]";

    // --- aileron increments, beta = 0 ----------------------------------------
    // ONE-SIDED in deflection. The configuration is mirror-symmetric about
    // its xz-plane, and mirroring maps an antisymmetric command of +delta
    // onto one of -delta while flipping the lateral wrench. So
    //
    //     dCY, dCl, dCn  are ODD in delta_a
    //     dCX, dCZ, dCm  are EVEN
    //
    // and the assembler mirrors, exactly as it does for sideslip. That is
    // an argument, not a measurement, so the sweep also solves one
    // NEGATIVE deflection and checks it -- recorded in "aileronMirror".
    //
    // Neutral and deflected are solved in the SAME pass, each warm-started
    // up its own alpha column. Differencing two limit-cycle means computed
    // from different continuation histories would add noise to the
    // increment that has nothing to do with the control surface.
    if (wantAileron) {
        out << ",\n\"aileron\":[\n";
        const Geometry::ControlSurface* aileron = nullptr;
        for (const Geometry::ControlSurface& cs : contract.ControlSurfaces)
            if (cs.Name == "aileron") { aileron = &cs; break; }
        if (!aileron) {
            std::cerr << "contract states no surface named 'aileron'\n";
            return 1;
        }
        const double tau = Geometry::FlapEffectiveness(1.0 - aileron->ChordFraction);
        std::cout << "\naileron: chord fraction " << aileron->ChordFraction << ", eta "
                  << aileron->EtaStart << "-" << aileron->EtaEnd
                  << ", thin-airfoil tau = " << tau << "\n"
                  << "alpha  delta     dCX        dCZ        dCl        dCn    iters\n";

        // Positive deflections only; plus one negative, for the mirror check.
        const std::vector<double> deltas = {5.0, 10.0, 20.0};
        const double mirrorProbeDeg = -10.0;

        struct Chain {
            double Delta;
            std::vector<S::StripSection> Strips;
            std::vector<double> Warm;
        };
        std::vector<Chain> chains;
        chains.push_back({0.0, StripsFromPanels(wing, halfSpan, contract.AirfoilSections), {}});
        for (const double d : deltas)
            chains.push_back(
                {d, StripsFromPanels(wing, halfSpan, contract.AirfoilSections, aileron, d), {}});
        chains.push_back({mirrorProbeDeg,
                          StripsFromPanels(wing, halfSpan, contract.AirfoilSections, aileron,
                                           mirrorProbeDeg),
                          {}});

        const double q = 0.5 * Rho * flightSpeed * flightSpeed;
        bool firstAil = true;
        fc.betaDeg = 0.0;
        for (const double alphaDeg : BuildAlphaGrid()) {
            fc.alphaDeg = alphaDeg;
            S::BodyAxisCoefficients neutral;
            for (Chain& chain : chains) {
                S::ViscousCouplingOptions opts = coupling;
                opts.InitialGamma = chain.Warm;
                const auto res = S::SolveViscousCoupled(wing, chain.Strips, fc, ref, trail, model,
                                                        opts, sources);
                chain.Warm = res.Base.gamma;
                const S::BodyAxisCoefficients w = S::BodyAxisFromCoupled(res, q, ref);
                if (chain.Delta == 0.0) {
                    neutral = w;
                    continue;
                }
                if (!res.Converged) ++unconverged;
                if (!firstAil) out << ",\n";
                firstAil = false;
                out << R"( {"alphaDeg":)" << alphaDeg << R"(,"deltaDeg":)" << chain.Delta
                    << R"(,"dCX":)" << (w.CX - neutral.CX) << R"(,"dCY":)" << (w.CY - neutral.CY)
                    << R"(,"dCZ":)" << (w.CZ - neutral.CZ) << R"(,"dCl":)" << (w.Cl - neutral.Cl)
                    << R"(,"dCm":)" << (w.Cm - neutral.Cm) << R"(,"dCn":)" << (w.Cn - neutral.Cn)
                    << R"(,"converged":)" << (res.Converged ? "true" : "false")
                    << R"(,"iterations":)" << res.Iterations << R"(,"residual":)"
                    << res.MaxResidual << '}';
                out.flush();

                std::cout << alphaDeg << "\t" << chain.Delta << "\t" << (w.CX - neutral.CX) << "\t"
                          << (w.CZ - neutral.CZ) << "\t" << (w.Cl - neutral.Cl) << "\t"
                          << (w.Cn - neutral.Cn) << "\t" << res.Iterations
                          << (res.Converged ? "" : "  (cycle-mean)") << std::endl;
            }
        }
        out << "\n],\n\"aileronMirrorProbeDeg\":" << mirrorProbeDeg
            << ",\n\"aileronTau\":" << tau
            << ",\n\"aileronEtaStart\":" << aileron->EtaStart
            << ",\n\"aileronChordFraction\":" << aileron->ChordFraction;
    }

    out << "\n}\n";
    std::cout << "\nwrote " << outPath << " (" << unconverged
              << " conditions exported as cycle means)\n";
    return 0;
}
