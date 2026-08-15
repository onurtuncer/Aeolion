// AileronEffectivenessExport.cpp -- the aeroDCl(alpha) question of the
// DAVE-ML model (models/README.md open item): does aileron roll authority
// decay past separation the way it must, or does the tabulation hand a
// simulator attached-flow roll power at 40 degrees incidence?
//
// The driver answers it by solving the SAME antisymmetric deflection two
// ways at every attitude and reporting both:
//
//   1. INVISCID -- the plain lattice (Solve/SolveWithSystem). The aileron
//      enters as rotated panel normals, so this is the geometric control
//      effect with no stall anywhere in it. It should NOT decay: a linear
//      lattice knows nothing about separation, and any decay here would be
//      induced/geometric rather than viscous.
//   2. LEVEL-2 COUPLED -- the sectional-feedback solve, the path that
//      actually generates the model's force tables.
//
// The comparison is the measurement, and it tests a specific structural
// worry rather than just producing a curve. StripSection::Alpha0Deg is
// read from the section's CST camber line
// (Geometry::SectionZeroLiftAngleDeg) and knows nothing about a deflected
// flap, while the coupling poses the section model in incidence measured
// FROM that zero-lift line (ViscousCoupling.h: alphaEffDeg -
// strip.Alpha0Deg) and drives the lattice's cl onto the section's. If the
// section model cannot see the flap, the fixed point may erase the very
// control effect the lattice geometry carries -- and it would do so at
// EVERY incidence, not merely past stall. That is a different and more
// serious defect than a missing post-stall decay, so the sweep starts in
// the attached range where the answer is unambiguous.
//
// The analytic section polar is used deliberately: the structural
// question -- whether a deflection reaches the section model at all -- is
// independent of which polar is mounted, and the analytic blend carries
// its own stall and deep-stall behaviour, so a genuine decay would still
// show. The anchored model's separation tables would add cost without
// changing what is being asked.
//
// ---------------------------------------------------------------------
// MEASURED (2026-08-15), and the answer is not the one the question
// expected. The Level-2 coupled path CANNOT REPRESENT AN AILERON AT ALL:
//
//   alpha   dCl (8 rows)   dCl (1 row, inviscid)   dCl (coupled)
//       0      0.009254                        0               0
//       6      0.009361                        0               0
//      18      0.008959                        0               0
//      60      0.003182                        0               0
//
// Both single-row columns are EXACTLY zero at every attitude, and the
// cause is a collision of two contracts that are each individually
// reasonable. PanelBuilder's MinRowsToResolveHinge = 2: a hinge line
// needs at least two chordwise rows to be split across, or
// ChordwiseRowBounds silently returns the undivided strip and the
// deflection is never applied. SolveViscousCoupled requires EXACTLY ONE
// Weissinger row per strip. The two cannot both be satisfied, so any
// aeroDC* table generated from the coupled path is identically zero --
// a flight model with no roll control, and nothing in the pipeline
// flags it.
//
// The inviscid multi-row column is a real control effect but is not a
// substitute: its decay from 0.00936 at alpha = 6 to 0.00318 at 60 is
// purely geometric (freestream projection and induced effects in a
// lattice that contains no stall). It passes straight through the
// configuration's stall near 18 degrees with no break at all, so a
// table built from it would hand a simulator most of its attached-flow
// roll authority deep into stall -- the departure-recovery failure the
// DAVE-ML report warns about, in the opposite direction.
//
// Neither available source is usable as it stands. The fix is
// solver-side, not driver-side: a deflected strip must carry the flap in
// its SECTION description, i.e. StripSection::Alpha0Deg shifted by the
// thin-airfoil flap increment for its hinge position and deflection (and
// ideally a flap-shifted stall angle), so the section model the coupling
// drives the lattice onto actually knows the flap is down. That is a
// behaviour change to a tested module and wants its own branch.
// ---------------------------------------------------------------------
//
// Usage:
//   aeolion_aileron_effectiveness <handoff.json> <out.json> [deltaDeg] [Vinf]

#include "Aeolion/Geometry/CstSurface.h"
#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

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

constexpr double Rho = 1.225;
constexpr double TrailSpans = 50.0;
constexpr int BodySectors = 16; // the measured value; see AttachmentSweepExport.cpp

// Attached through deep stall, on the model's own alpha breakpoints.
std::vector<double> BuildAlphaGrid() {
    std::vector<double> alphas;
    for (double a = 0.0; a <= 26.0 + 1e-9; a += 2.0) alphas.push_back(a);
    for (const double a : {30.0, 40.0, 50.0, 60.0}) alphas.push_back(a);
    return alphas;
}

// The true chord frame, NOT the cambered panel axes -- the camber
// double-count trap documented in ViscousCoupling.h and measured twice in
// this repo. The handoff wing is rectangular, unswept and untwisted, so
// its chord frame is the solver frame.
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

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: aeolion_aileron_effectiveness <handoff.json> <out.json>"
                     " [deltaDeg] [Vinf]\n";
        return 1;
    }
    const std::string handoffPath = argv[1];
    const std::string outPath = argv[2];
    const double deltaDeg = (argc > 3) ? std::atof(argv[3]) : 10.0;
    const double flightSpeed = (argc > 4) ? std::atof(argv[4]) : 25.0;

    Geometry::HandoffContract contract;
    try {
        contract = Geometry::LoadHandoff(handoffPath);
    } catch (const std::exception& error) {
        std::cerr << "cannot load " << handoffPath << ": " << error.what() << "\n";
        return 1;
    }
    const int nativeRows = contract.Mesh.ChordwisePanels;

    std::size_t aileronIndex = contract.ControlSurfaces.size();
    for (std::size_t i = 0; i < contract.ControlSurfaces.size(); ++i)
        if (contract.ControlSurfaces[i].Name == "aileron") { aileronIndex = i; break; }
    if (aileronIndex >= contract.ControlSurfaces.size()) {
        std::cerr << "contract states no surface named 'aileron'\n";
        return 1;
    }

    PB::LatticeOptions options;
    options.BodyCircumferentialPanels = BodySectors;
    options.CarryThroughLift = false; // every strip midpoint outside the body

    // Two lattices of the SAME wing at different chordwise resolution. The
    // native one (the contract's own row count) can resolve the hinge; the
    // single-row one is what the Level-2 coupling requires, and
    // PanelBuilder's MinRowsToResolveHinge = 2 means a hinge is simply not
    // split there. Building both is the measurement.
    Geometry::HandoffContract single = contract;
    single.Mesh.ChordwisePanels = 1;

    PB::LatticeBuilder nativeBuilder(contract, options);
    PB::LatticeBuilder singleBuilder(single, options);

    const auto body = singleBuilder.BuildBody();
    const auto duct = singleBuilder.BuildDuct();
    std::vector<Lattice::SourcePanel> sources = body;
    sources.insert(sources.end(), duct.begin(), duct.end());

    const auto nativeNeutral = nativeBuilder.Build();
    nativeBuilder.Deflect(PB::Antisymmetric(aileronIndex, deltaDeg));
    const auto nativeDeflected = nativeBuilder.Build();
    nativeBuilder.ClearDeflections();

    const auto wingNeutral = singleBuilder.Build();
    singleBuilder.Deflect(PB::Antisymmetric(aileronIndex, deltaDeg));
    const auto wingDeflected = singleBuilder.Build();
    singleBuilder.ClearDeflections();

    const double halfSpan = 0.5 * contract.Span;
    const auto stripsNeutral = StripsFromPanels(wingNeutral, halfSpan, contract.AirfoilSections);
    const auto stripsDeflected = StripsFromPanels(wingDeflected, halfSpan, contract.AirfoilSections);

    const double trail = TrailSpans * contract.Span;
    const auto preparedNeutral = S::Prepare(S::PanelSystem{wingNeutral, sources}, trail);
    const auto preparedDeflected = S::Prepare(S::PanelSystem{wingDeflected, sources}, trail);
    const auto preparedNativeNeutral = S::Prepare(S::PanelSystem{nativeNeutral, sources}, trail);
    const auto preparedNativeDeflected = S::Prepare(S::PanelSystem{nativeDeflected, sources}, trail);

    S::ReferenceGeometry ref;
    ref.Area = singleBuilder.GrossPlanformArea();
    ref.Span = contract.Span;
    ref.Chord = ref.Area / contract.Span;

    S::FreestreamConditions fc;
    fc.Vinf = flightSpeed;
    fc.rho = Rho;
    if (contract.MomentReferencePointStated)
        fc.RefPoint = S::Vec3(-contract.MomentReferencePoint.x, contract.MomentReferencePoint.y,
                              -contract.MomentReferencePoint.z);

    // Plain damped iteration: the Anderson-accelerated fixed point lands on
    // spurious sawtooth equilibria on this many tightly packed strips (the
    // measured Phase-0 finding).
    S::ViscousCouplingOptions coupling;
    coupling.Relaxation = 0.05;
    coupling.AndersonDepth = 0;
    coupling.MaxIterations = 1000;

    const S::AnalyticSectionModel model{};

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot open " << outPath << " for writing\n";
        return 1;
    }
    out.precision(9);
    out << "{\n\"meta\":{\"deltaDeg\":" << deltaDeg << ",\"Vinf\":" << flightSpeed
        << ",\"rho\":" << Rho << ",\"area\":" << ref.Area << ",\"span\":" << ref.Span
        << ",\"nativeChordwisePanels\":" << nativeRows << ",\"minRowsToResolveHinge\":"
        << PB::MinRowsToResolveHinge << ",\"sectionModel\":\"analytic\",\"aileronEtaStart\":"
        << contract.ControlSurfaces[aileronIndex].EtaStart << ",\"aileronChordFraction\":"
        << contract.ControlSurfaces[aileronIndex].ChordFraction
        << ",\"note\":\"dClInviscid is the geometric control effect; dClCoupled is what"
           " the Level-2 path that generates the model's tables actually produces\"},\n"
        << "\"rows\":[\n";

    std::cout << "delta = " << deltaDeg << " deg;  native rows = " << nativeRows
              << " (hinge resolved), single row = 1 (the coupling's contract)\n"
              << "alpha   dCl_native   dCl_1row_inv   dCl_coupled   coupled/native\n";

    bool first = true;
    for (const double alphaDeg : BuildAlphaGrid()) {
        fc.alphaDeg = alphaDeg;

        const auto natNeutral = S::SolveWithSystem(preparedNativeNeutral, fc, ref);
        const auto natDeflected = S::SolveWithSystem(preparedNativeDeflected, fc, ref);
        const double dClNative = natDeflected.Croll - natNeutral.Croll;

        const auto invNeutral = S::SolveWithSystem(preparedNeutral, fc, ref);
        const auto invDeflected = S::SolveWithSystem(preparedDeflected, fc, ref);
        const double dClInviscid = invDeflected.Croll - invNeutral.Croll;

        const auto cplNeutral = S::SolveViscousCoupled(wingNeutral, stripsNeutral, fc, ref, trail,
                                                       model, coupling, sources);
        const auto cplDeflected = S::SolveViscousCoupled(wingDeflected, stripsDeflected, fc, ref,
                                                         trail, model, coupling, sources);
        const double dClCoupled = cplDeflected.Base.Croll - cplNeutral.Base.Croll;

        // Retention against the resolved-hinge lattice: what fraction of the
        // real geometric control effect each single-row path still carries.
        const double ratio = (std::fabs(dClNative) > 1e-12) ? dClCoupled / dClNative : 0.0;
        const double ratioInviscid =
            (std::fabs(dClNative) > 1e-12) ? dClInviscid / dClNative : 0.0;

        if (!first) out << ",\n";
        first = false;
        out << R"( {"alphaDeg":)" << alphaDeg << R"(,"dClNative":)" << dClNative
            << R"(,"dClInviscid":)" << dClInviscid << R"(,"dClCoupled":)" << dClCoupled
            << R"(,"ratio":)" << ratio << R"(,"ratioInviscid":)" << ratioInviscid
            << R"(,"ClInviscidNeutral":)" << invNeutral.Croll << R"(,"ClCoupledNeutral":)"
            << cplNeutral.Base.Croll << R"(,"CLInviscid":)" << invNeutral.CL
            << R"(,"CLCoupled":)" << cplNeutral.Base.CL << R"(,"coupledConverged":)"
            << ((cplNeutral.Converged && cplDeflected.Converged) ? "true" : "false")
            << R"(,"coupledResidual":)"
            << std::max(cplNeutral.MaxResidual, cplDeflected.MaxResidual) << '}';
        out.flush();

        std::cout << alphaDeg << "\t" << dClNative << "\t" << dClInviscid << "\t" << dClCoupled
                  << "\t" << ratio
                  << (cplNeutral.Converged && cplDeflected.Converged ? "" : "  (cycle-mean)")
                  << std::endl;
    }

    out << "\n]}\n";
    std::cout << "wrote " << outPath << "\n";
    return 0;
}
