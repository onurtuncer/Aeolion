// ParticleCrosscheckExport.cpp -- the post-separation study's 3-D unsteady
// cross-check (Part II, "The unsteady cross-check"): the particle-wake solve of
// Solver/ParticleWake.h run at the declared handful of attitudes,
// alpha = 30/60/90 x beta = 0/15, on the geometry handoff's wing.
//
// Three deliberate reductions, stated here because they shape the numbers:
//
//   - The wing is FULL SPAN and contiguous, built here from the contract's
//     own planform law (chord over eta; this wing is rectangular, unswept
//     and untwisted -- Part I's application geometry) rather than from the
//     trimmed builder lattice. MEASURED reason: the trimmed wing's two
//     inner tip-vortex streams face each other across the body gap with no
//     body there to occupy it, and their close-range dynamics dominated
//     the runs -- RMS growing with the averaging window instead of
//     converging (13.5 on an O(1) mean at alpha = 30, duration 30).
//     Absent the body, span continuity is the physical statement, and it
//     is also what the map's carry-through solve asserts.
//   - The lattice is COARSE (12 strips): a particle wake costs N^2 and
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
#include "Aeolion/Solver/ParticleWake.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace Aeolion;
namespace S = Aeolion::Solver;

namespace {

constexpr double FlightSpeed = 25.0;
constexpr double Rho = 1.225;
constexpr int FullSpanStrips = 12;

const double Alphas[] = {30.0, 60.0, 90.0};
const double Betas[] = {0.0, 15.0};

} // namespace

int main(int argc, char** argv) {
    const std::string handoffPath =
        (argc > 1) ? argv[1] : "tests/Data/AeolionGeometryHandoff-1.8.0.json";
    const std::string outPath = (argc > 2) ? argv[2] : "particle-crosscheck.json";
    // Averaging duration [convective times]: the paper's long-averaging
    // runs pass a larger value here; the default matches the solver's.
    const double duration = (argc > 3) ? std::atof(argv[3]) : S::PwDefaultDuration;
    // Pilot mode: argv[4]/argv[5] restrict the matrix to one attitude, and
    // argv[6] overrides the eddy-viscosity coefficient -- the Phase-A
    // convergence pilots sweep it (0 recovers the inviscid tier).
    const bool pilot = argc > 5;
    const double alphaOnly = pilot ? std::atof(argv[4]) : 0.0;
    const double betaOnly = pilot ? std::atof(argv[5]) : 0.0;
    const double nuCoeff =
        (argc > 6) ? std::atof(argv[6]) : S::PwTurbulentViscosityCoeff;

    Geometry::HandoffContract contract;
    try {
        contract = Geometry::LoadHandoff(handoffPath);
    } catch (const std::exception& error) {
        std::cerr << "cannot load " << handoffPath << ": " << error.what() << "\n";
        return 1;
    }
    // Full-span contiguous lattice from the contract's planform law (see
    // the header note): chord interpolated over |eta| from the stations,
    // straight quarter-chord line, true-chord frames.
    const double halfSpan = 0.5 * contract.Span;
    const auto chordAt = [&](double eta) {
        const auto& st = contract.Stations;
        if (st.empty()) return 0.0;
        if (eta <= st.front().Eta) return st.front().Chord;
        for (std::size_t k = 1; k < st.size(); ++k)
            if (eta <= st[k].Eta) {
                const double f = (eta - st[k - 1].Eta) /
                                 std::max(st[k].Eta - st[k - 1].Eta, 1e-12);
                return st[k - 1].Chord + f * (st[k].Chord - st[k - 1].Chord);
            }
        return st.back().Chord;
    };

    std::vector<S::Panel> wing;
    std::vector<S::StripSection> strips;
    double grossArea = 0.0;
    const double width = contract.Span / FullSpanStrips;
    for (int i = 0; i < FullSpanStrips; ++i) {
        const double y0 = -halfSpan + width * i;
        const double etaMid = std::fabs(y0 + 0.5 * width) / halfSpan;
        const double chord = chordAt(etaMid);

        S::Panel p;
        p.A = S::Vec3(0.25 * chord, y0, 0.0);
        p.B = S::Vec3(0.25 * chord, y0 + width, 0.0);
        p.ControlPoint = S::Vec3(0.75 * chord, y0 + 0.5 * width, 0.0);
        p.Normal = S::Vec3(0.0, 0.0, 1.0);
        p.TrailDirA = p.TrailDirB = S::Vec3(1.0, 0.0, 0.0);
        p.PlanformArea = chord * width;
        p.Area = p.PlanformArea;
        p.SpanwiseWidth = width;
        p.Surface = "wing";
        wing.push_back(p);

        S::StripSection strip;
        strip.ChordDir = S::Vec3(1.0, 0.0, 0.0);
        strip.LiftDir = S::Vec3(0.0, 0.0, 1.0);
        strip.Chord = chord;
        strip.Width = width;
        strip.Eta = etaMid;
        strip.Alpha0Deg = Geometry::SectionZeroLiftAngleDeg(contract.AirfoilSections, etaMid);
        strips.push_back(strip);
        grossArea += chord * width;
    }

    S::ReferenceGeometry ref;
    ref.Area = grossArea;
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
    out << "{\n\"meta\":{\"strips\":" << wing.size() << ",\"fullSpan\":true,\"area\":" << ref.Area
        << ",\"span\":" << contract.Span << ",\"Vinf\":" << FlightSpeed << "},\n\"cross-check\":[\n";

    std::cerr << "particle cross-check: " << wing.size() << " strips, duration " << duration
              << " convective times\n"; // cerr: unbuffered progress under redirection
    bool first = true;
    for (const double alphaDeg : Alphas) {
        for (const double betaDeg : Betas) {
            if (pilot && (alphaDeg != alphaOnly || betaDeg != betaOnly)) continue;
            fc.alphaDeg = alphaDeg;
            fc.betaDeg = betaDeg;
            S::ParticleWakeOptions wakeOptions;
            wakeOptions.Duration = duration;
            wakeOptions.TurbulentViscosityCoeff = nuCoeff;
            const S::ParticleWakeResult run =
                S::SolveParticleWake(wing, strips, fc, ref, wakeOptions);
            std::cerr << "alpha=" << alphaDeg << " beta=" << betaDeg << " nu=" << nuCoeff
                      << ": mean CL " << run.MeanCL << " CD " << run.MeanCD << " CN "
                      << run.MeanCN << " +-" << run.MeanCN_CI << " (" << run.Batches
                      << " batches)  rms CL " << run.RmsCL << " CN " << run.RmsCN
                      << "  particles " << run.MaxParticles << " merged " << run.Merged
                      << (run.Valid ? "" : "  [INVALID]") << std::endl;
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
