// ParticleCrosscheckExport.cpp -- the post-separation study's 3-D unsteady
// cross-check (Part II, "The unsteady cross-check"): the particle-wake solve of
// Solver/ParticleWake.h run at the declared handful of attitudes,
// alpha = 30/60/90 x beta = 0/15, on the geometry handoff's wing.
//
// Two deliberate reductions, stated here because they shape the numbers:
//
//   - The lattice is COARSENED spanwise (SpanwisePanelsPerSection = 2,
//     one Weissinger row) before building: a particle wake costs N^2 and
//     the cross-check's question -- what do the quasi-steady cycle means
//     leave out in mean and fluctuation -- is posed at configuration
//     level, not at the map's spanwise resolution. The comparison target
//     is the map's WING-ONLY forces at the same attitudes.
//   - No body or duct sources, and the separation decision uses the
//     incidence-threshold fallback: at 30 degrees and beyond every strip
//     is far past its boundary, so the f-tables would decide nothing
//     the threshold does not.
//
// Output: JSON rows [alphaDeg, betaDeg, meanCL, meanCD, meanCY, meanCN,
// rmsCL, rmsCN, circulationCL, maxParticles], consumed by
// papers/journal-of-aircraft-poststall/figures/render-crosscheck-figure.py.

#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/Geometry/CstSurface.h"
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/ParticleWake.h"

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
constexpr int CoarseSpanwisePerSection = 1;

const double Alphas[] = {30.0, 60.0, 90.0};
const double Betas[] = {0.0, 15.0};

std::vector<S::StripSection> StripsFromPanels(const std::vector<S::Panel>& panels, double halfSpan,
                                              const std::vector<Geometry::AirfoilSection>& sections) {
    // True-chord frames, as everywhere since the camber double-count.
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
    const std::string handoffPath =
        (argc > 1) ? argv[1] : "tests/Data/AeolionGeometryHandoff-1.8.0.json";
    const std::string outPath = (argc > 2) ? argv[2] : "particle-crosscheck.json";
    // Averaging duration [convective times]: the paper's long-averaging
    // runs pass a larger value here; the default matches the solver's.
    const double duration = (argc > 3) ? std::atof(argv[3]) : S::PwDefaultDuration;

    Geometry::HandoffContract contract;
    try {
        contract = Geometry::LoadHandoff(handoffPath);
    } catch (const std::exception& error) {
        std::cerr << "cannot load " << handoffPath << ": " << error.what() << "\n";
        return 1;
    }
    contract.Mesh.ChordwisePanels = 1;
    contract.Mesh.SpanwisePanelsPerSection = CoarseSpanwisePerSection;

    PB::LatticeOptions options;
    options.CarryThroughLift = false;
    PB::LatticeBuilder builder(contract, options);
    const auto wing = builder.Build();
    const double halfSpan = 0.5 * contract.Span;
    const auto strips = StripsFromPanels(wing, halfSpan, contract.AirfoilSections);

    S::ReferenceGeometry ref;
    ref.Area = builder.GrossPlanformArea();
    ref.Span = contract.Span;
    ref.Chord = ref.Area / contract.Span;

    S::FreestreamConditions fc;
    fc.Vinf = FlightSpeed;
    fc.rho = Rho;

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot open " << outPath << "\n";
        return 1;
    }
    out.precision(9);
    out << "{\n\"meta\":{\"strips\":" << wing.size() << ",\"area\":" << ref.Area
        << ",\"span\":" << contract.Span << ",\"Vinf\":" << FlightSpeed << "},\n\"cross-check\":[\n";

    std::cerr << "particle cross-check: " << wing.size() << " strips, duration " << duration
              << " convective times\n"; // cerr: unbuffered progress under redirection
    bool first = true;
    for (const double alphaDeg : Alphas) {
        for (const double betaDeg : Betas) {
            fc.alphaDeg = alphaDeg;
            fc.betaDeg = betaDeg;
            S::ParticleWakeOptions wakeOptions;
            wakeOptions.Duration = duration;
            const S::ParticleWakeResult run =
                S::SolveParticleWake(wing, strips, fc, ref, wakeOptions);
            std::cerr << "alpha=" << alphaDeg << " beta=" << betaDeg << ": mean CL "
                      << run.MeanCL << " CD " << run.MeanCD << " CN " << run.MeanCN
                      << "  rms CL " << run.RmsCL << " CN " << run.RmsCN << "  particles "
                      << run.MaxParticles << (run.Valid ? "" : "  [INVALID]") << std::endl;
            if (!run.Valid) continue;
            if (!first) out << ",\n";
            first = false;
            out << " [" << alphaDeg << ',' << betaDeg << ',' << run.MeanCL << ',' << run.MeanCD
                << ',' << run.MeanCY << ',' << run.MeanCN << ',' << run.RmsCL << ','
                << run.RmsCN << ',' << run.MeanCirculationCL << ',' << run.MaxParticles << ']';
            out.flush();
        }
    }
    out << "\n]\n}\n";
    std::cout << "wrote " << outPath << "\n";
    return 0;
}
