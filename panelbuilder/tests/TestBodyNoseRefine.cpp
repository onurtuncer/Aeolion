// LatticeOptions::BodyNoseRefineStations: axial mesh resolution as a consumer
// choice, the counterpart to BodyCircumferentialPanels.
//
// The property that matters is NOT that refinement adds panels -- that is
// obvious and would be true of a wrong implementation too. It is that
// refinement changes the DISCRETIZATION and not the GEOMETRY. Added stations
// are evaluated on the contract's own piecewise-linear radius law, so the
// surface a consumer gets back is the same surface, more finely cut. An
// implementation that resampled onto a smooth curve, or that dropped original
// breakpoints while adding new ones, would cut corners off the shape and pass
// any test that only counted panels.
//
// The check is on panel CORNERS, not centroids. Corners are placed at
// (x, RadiusAt(x)) exactly, so they must lie on the law to round-off. A flat
// quad's centroid does NOT -- it chords the ring azimuthally and the surface
// axially -- so a centroid test measures the panelling, not the geometry, and
// reports large deviations at the end caps where rings degenerate.
#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/PanelBuilder/PanelBuilder.h"

#include <cmath>
#include <iostream>
#include <string>

static int failures = 0;
#define CHECK(cond, msg) \
    do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace PB = Aeolion::PanelBuilder;
namespace G = Aeolion::Geometry;

int main() {
    const G::HandoffContract contract =
        G::LoadHandoff(std::string(AEOLION_TEST_DATA_DIR) + "/AeolionGeometryHandoff-1.8.0.json");
    CHECK(contract.Body.IsPresent(), "fixture must carry a body");

    PB::LatticeOptions plain;
    plain.BodyCircumferentialPanels = 16;
    PB::LatticeOptions refined = plain;
    refined.BodyNoseRefineStations = 20;
    refined.BodyNoseRefineFraction = 0.10;

    const auto bodyPlain = PB::LatticeBuilder(contract, plain).BuildBody();
    const auto bodyRefined = PB::LatticeBuilder(contract, refined).BuildBody();

    // 1. It refines.
    CHECK(bodyRefined.size() > bodyPlain.size(), "refinement added no panels");

    // 2. It does not move the surface. Every corner must sit ON the contract's
    //    radius law, with one legitimate exception: the body is CAPPED at both
    //    ends, and a cap is panelled as concentric annuli, so its corners sit
    //    at every intermediate radius from the axis out to the rim -- all at
    //    ONE axial station. RadiusAt is single-valued in x and cannot describe
    //    them, but they are surface points all the same. So the rule is: a
    //    corner is on the law, or it is a cap point -- at an end station, no
    //    further out than the rim there.
    //
    //    Measured rather than assumed. Off-cap corners match the law to
    //    ~2e-18; the only departures are at contract x = -0.491, the base
    //    station, at radii 0.0406/4 apart -- the annuli. The unrefined mesh
    //    shows exactly the same thing, which is what proves the departure
    //    belongs to the cap and not to the refinement.
    const double xNose = contract.Body.Stations.front().x;
    const double xTail = contract.Body.Stations.back().x;
    const auto offLaw = [&](const auto& mesh) {
        double worst = 0.0;
        for (const auto& p : mesh)
            for (const auto& c : p.Corners) {
                const double xc = -c.x;
                const double r = std::sqrt(c.y * c.y + c.z * c.z);
                const double rLaw = G::RadiusAt(contract.Body, xc);
                const bool onCap = (std::fabs(xc - xNose) < 1e-9 || std::fabs(xc - xTail) < 1e-9) &&
                                   r <= rLaw + 1e-9;
                if (onCap) continue;
                worst = std::max(worst, std::fabs(r - rLaw));
            }
        return worst;
    };
    const double worstPlain = offLaw(bodyPlain);
    const double worst = offLaw(bodyRefined);
    CHECK(worst < 1e-9, "refined panel corners left the radius law by " << worst << " m");
    CHECK(worstPlain < 1e-9, "unrefined corners already off the law by " << worstPlain << " m");
    CHECK(worst <= worstPlain + 1e-12,
          "refinement made the surface worse: " << worstPlain << " -> " << worst);

    // 3. Zero means untouched, so the default cannot change any existing
    //    consumer's mesh.
    PB::LatticeOptions off = plain;
    off.BodyNoseRefineStations = 0;
    CHECK(PB::LatticeBuilder(contract, off).BuildBody().size() == bodyPlain.size(),
          "zero refinement must be a no-op");

    if (failures == 0)
        std::cout << "TestBodyNoseRefine: " << bodyPlain.size() << " -> " << bodyRefined.size()
                  << " panels, corners on the radius law to " << worst << " m -- OK\n";
    return failures == 0 ? 0 : 1;
}
