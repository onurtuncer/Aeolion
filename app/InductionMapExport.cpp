// InductionMapExport.cpp -- the coupling* half of the DAVE-ML flight
// model (models/README.md): what the aft fan does to the AIRFRAME, as
// increments from the power-off tables.
//
// --- why this driver exists, and why it was long declared blocked ---------
// The model's airframe tables are power-off. On this configuration that is
// a real omission rather than a conservative one: the ducted fan sits AFT
// of the wing -- duct leading edge about a fifth of a chord behind the
// wing trailing edge -- so it does not blow the wing. What it does is
// induce a favourable pressure gradient UPSTREAM of itself, over the wing,
// which delays separation. And it does so over precisely the span that
// sheds first: duct outer radius to semi-span is 0.21, against a worst
// separation station at |2y/b| ~ 0.15-0.23.
//
// Solver::SlipstreamField cannot compute this. It is a momentum-theory
// wake and returns identically zero ahead of the disk by construction, so
// pointed at a wing that sits upstream it reports exactly zero
// interaction -- a property of the model, not of the aircraft, and the
// most dangerous kind of answer because it is confidently null.
//
// THE BLOCKER IS STALE. Solver/DiskInduction.h implements the missing
// half: a uniformly loaded actuator disk is exactly equivalent to a
// semi-infinite cylindrical vortex sheet, which has a computable field
// everywhere including upstream, and an exact closed form on the axis to
// check against. It is pinned by TestDiskInduction and already used by
// the attachment sweep. Nothing further was needed to generate these
// tables; the specification simply had not caught up with the solver.
//
// --- what is computed ----------------------------------------------------
// For each (alpha, Tc) the fan's thrust follows from the interaction
// index itself, T = Tc * qbar * S, momentum theory gives the induced
// velocity at the disk, and the vortex cylinder gives the field. The
// coupled solve then runs with that field in its externalField hook and
// the result is differenced against the SAME solve with the fan off.
//
// Neutral and powered are solved in one pass, each warm-started up its own
// alpha column, for the reason the aileron sweep gives: differencing two
// limit-cycle means taken from different continuation histories injects
// noise into the increment that has nothing to do with the fan.
//
// --- the annulus ---------------------------------------------------------
// The fan is an annulus around the tail boom, not a disk, so the hub
// radius is the body's own radius at the duct station. A uniformly loaded
// annulus sheds at BOTH edges -- outer at +gamma_t, inner at -gamma_t --
// which DiskInduction.h represents as two superposed cylinders.
//
// --- declared limits -----------------------------------------------------
// Uniform disk loading (real loading tapers at both edges); no swirl,
// which is correct here because the rotor's axial wake vorticity lives
// DOWNSTREAM of the disk and contributes nothing ahead of it; and the
// q-normalized Tc degenerates toward hover, so the tables are declared
// valid for V >= 10 m/s.
//
// Usage:
//   aeolion_induction_map <handoff.json> <out.json> [Vinf] [relaxation]
//                         [maxIterations]

#include "Aeolion/Geometry/CstSurface.h"
#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/BodyAxes.h"
#include "Aeolion/Solver/DiskInduction.h"
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

constexpr double FlightSpeed = 25.0;
constexpr double Rho = 1.225;
constexpr double TrailSpans = 50.0;
constexpr int BodySectors = 16;

// The model's tcBp. Logarithmically spaced because momentum theory makes
// the induced-velocity ratio go as sqrt(Tc), so equal RATIOS of Tc are
// roughly equal increments of effect.
const double ThrustCoefficients[] = {0.5, 1.0, 2.0, 4.0, 8.0};

std::vector<double> BuildAlphaGrid() {
    std::vector<double> alphas;
    for (double a = -4.0; a <= 26.0 + 1e-9; a += 2.0) alphas.push_back(a);
    for (const double a : {30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0})
        alphas.push_back(a);
    return alphas;
}

// True chord frame, not the cambered panel axes -- the camber double-count
// trap ViscousCoupling.h warns about and this repo has hit twice.
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

/**
 * The fan as an annular actuator disk, in SOLVER axes. The contract states
 * the duct centre in its own frame, so the ingest flip applies here as it
 * does to the moment reference point.
 */
S::ActuatorDisk MakeFanDisk(const Geometry::HandoffContract& contract, double thrust, double rho,
                            double axialSpeed) {
    S::ActuatorDisk disk;
    if (!(thrust > 0.0) || !contract.Duct.IsStated || !(contract.Propulsion.DiskRadius > 0.0))
        return disk;

    const double diskContractX = contract.Duct.Center.x;
    disk.Center = S::Vec3(-diskContractX, contract.Duct.Center.y, -contract.Duct.Center.z);
    disk.Axis = S::Vec3(1.0, 0.0, 0.0); // downstream is solver +x (aft)
    disk.Radius = contract.Propulsion.DiskRadius;
    // Annular: the tail boom occupies the middle, so the hub radius is the
    // body's own radius at the duct station.
    disk.HubRadius = Geometry::RadiusAt(contract.Body, diskContractX);
    if (!(disk.HubRadius < disk.Radius)) disk.HubRadius = 0.0;
    disk.InducedVelocity = S::InducedVelocityFromThrust(thrust, rho, disk.Area(), axialSpeed);
    return disk;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: aeolion_induction_map <handoff.json> <out.json> [Vinf]"
                     " [relaxation] [maxIterations]\n";
        return 1;
    }
    const std::string handoffPath = argv[1];
    const std::string outPath = argv[2];
    const double flightSpeed = (argc > 3) ? std::atof(argv[3]) : FlightSpeed;
    const double relaxation = (argc > 4) ? std::atof(argv[4]) : 0.05;
    const int maxIterations = (argc > 5) ? std::atoi(argv[5]) : 1000;

    Geometry::HandoffContract contract;
    try {
        contract = Geometry::LoadHandoff(handoffPath);
    } catch (const std::exception& error) {
        std::cerr << "cannot load " << handoffPath << ": " << error.what() << "\n";
        return 1;
    }
    if (!contract.Duct.IsStated || !(contract.Propulsion.DiskRadius > 0.0)) {
        std::cerr << "contract states no duct or disk radius; nothing to induce with\n";
        return 1;
    }
    contract.Mesh.ChordwisePanels = 1;

    PB::LatticeOptions options;
    options.BodyCircumferentialPanels = BodySectors;
    options.CarryThroughLift = false;
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
    fc.betaDeg = 0.0;
    if (contract.MomentReferencePointStated)
        fc.RefPoint = S::Vec3(-contract.MomentReferencePoint.x, contract.MomentReferencePoint.y,
                              -contract.MomentReferencePoint.z);

    S::PostStallSectionModel anchored;
    anchored.AspectRatio = aspectRatio;
    {
        auto tables =
            S::BuildSeparationTables(prepared, wing, strips, contract.AirfoilSections, fc, ref);
        std::vector<double> etas;
        etas.reserve(strips.size());
        for (const S::StripSection& s : strips) etas.push_back(s.Eta);
        anchored.SeparationPoint = S::MakeSeparationFunction(std::move(tables), std::move(etas));
    }
    const S::SectionModel model = anchored;

    S::ViscousCouplingOptions coupling;
    coupling.Relaxation = relaxation;
    coupling.AndersonDepth = 0;
    coupling.MaxIterations = maxIterations;

    const double q = 0.5 * Rho * flightSpeed * flightSpeed;
    const double qS = q * ref.Area;

    // Report the geometry that makes this a real effect rather than a
    // rounding error, so the run's own log states its premise.
    {
        const S::ActuatorDisk probe = MakeFanDisk(contract, qS, Rho, flightSpeed);
        const double ductLEx = -(contract.Duct.Center.x + 0.5 * contract.Duct.Chord);
        std::cout << "fan: R=" << probe.Radius << " hub=" << probe.HubRadius
                  << " centre(solver x)=" << probe.Center.x
                  << "  R/semi-span=" << probe.Radius / halfSpan
                  << "  duct LE at solver x=" << ductLEx << "\n";
    }

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot open " << outPath << " for writing\n";
        return 1;
    }
    out.precision(9);
    out << "{\n\"meta\":{\"designId\":\"" << contract.DesignId << "\",\"schema\":\""
        << contract.SchemaVersion << "\",\"Vinf\":" << flightSpeed << ",\"rho\":" << Rho
        << ",\"area\":" << ref.Area << ",\"span\":" << ref.Span << ",\"chord\":" << ref.Chord
        << ",\"diskRadius\":" << contract.Propulsion.DiskRadius
        << ",\"model\":\"semi-infinite vortex cylinder (Solver/DiskInduction.h)\""
        << ",\"validityMinSpeedMps\":10"
        << ",\"excludes\":\"uniform disk loading; no swirl (the rotor's axial wake"
           " vorticity is downstream of the disk and contributes nothing ahead of"
           " it)\"},\n\"rows\":[\n";

    std::cout << "\nalpha   Tc      vi     dCX         dCZ         dCm      iters\n";

    struct Chain {
        double Tc;
        std::vector<double> Warm;
    };
    std::vector<Chain> chains;
    chains.push_back({0.0, {}});
    for (const double tc : ThrustCoefficients) chains.push_back({tc, {}});

    bool first = true;
    int unconverged = 0;
    for (const double alphaDeg : BuildAlphaGrid()) {
        fc.alphaDeg = alphaDeg;
        S::BodyAxisCoefficients off;
        for (Chain& chain : chains) {
            const double thrust = chain.Tc * qS;
            const S::ActuatorDisk disk = MakeFanDisk(contract, thrust, Rho, flightSpeed);
            std::function<S::Vec3(const S::Vec3&)> field;
            if (chain.Tc > 0.0 && disk.Valid()) field = S::DiskInductionField(disk);

            S::ViscousCouplingOptions opts = coupling;
            opts.InitialGamma = chain.Warm;
            const auto start = std::chrono::steady_clock::now();
            const auto res = S::SolveViscousCoupled(wing, strips, fc, ref, trail, model, opts,
                                                    sources, field);
            const double seconds =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            chain.Warm = res.Base.gamma;
            const S::BodyAxisCoefficients w = S::BodyAxisFromCoupled(res, q, ref);

            if (chain.Tc == 0.0) {
                off = w;
                continue;
            }
            if (!res.Converged) ++unconverged;
            if (!first) out << ",\n";
            first = false;
            out << R"( {"alphaDeg":)" << alphaDeg << R"(,"Tc":)" << chain.Tc
                << R"(,"thrustN":)" << thrust << R"(,"viMps":)" << disk.InducedVelocity
                << R"(,"dCX":)" << (w.CX - off.CX) << R"(,"dCY":)" << (w.CY - off.CY)
                << R"(,"dCZ":)" << (w.CZ - off.CZ) << R"(,"dCl":)" << (w.Cl - off.Cl)
                << R"(,"dCm":)" << (w.Cm - off.Cm) << R"(,"dCn":)" << (w.Cn - off.Cn)
                << R"(,"converged":)" << (res.Converged ? "true" : "false")
                << R"(,"iterations":)" << res.Iterations << R"(,"residual":)"
                << res.MaxResidual
                << R"(,"cycleFluctuation":)" << res.CycleFluctuation()
                << R"(,"seconds":)" << seconds << '}';
            out.flush();

            std::cout << alphaDeg << "\t" << chain.Tc << "\t" << disk.InducedVelocity << "\t"
                      << (w.CX - off.CX) << "\t" << (w.CZ - off.CZ) << "\t" << (w.Cm - off.Cm)
                      << "\t" << res.Iterations << (res.Converged ? "" : "  (cycle-mean)")
                      << std::endl;
        }
    }

    out << "\n]}\n";
    std::cout << "\nwrote " << outPath << " (" << unconverged
              << " conditions exported as cycle means)\n";
    return 0;
}
