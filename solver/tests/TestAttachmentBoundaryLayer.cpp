// TestAttachmentBoundaryLayer.cpp -- validates the separation march against
// closed-form results, not against itself.
//
//   - Flat plate. With no pressure gradient Thwaites must reproduce
//     Blasius' momentum thickness, theta = 0.664 s/sqrt(Re_s); Thwaites'
//     own constant gives 0.671, a known 1% high. The layer must never
//     separate.
//
//   - Howarth's linearly retarded flow, Ue = U0 (1 - s/L). This is THE
//     textbook test of a Thwaites march: the exact answer separates at
//     s/L = 0.120 and Thwaites' correlation is known to place it at
//     0.123. Getting this right is what makes the separation locations
//     this module reports mean anything. The march treats that crossing as
//     a bubble and transitions there, so the checked quantity is the
//     BUBBLE location -- the same lambda = -0.09 crossing under the name
//     the module reports it by.
//
//   - The stagnation start. A run beginning at a real attachment point
//     (Ue ~ a s) must produce theta_0 = sqrt(0.075 nu / a) on its own,
//     with no seed -- the same number StagnationMomentumThickness computes
//     from the Hiemenz similarity solution by a completely different
//     route. The two agreeing is what says the march is starting where it
//     claims to.
//
//   - A favourable gradient must NOT separate, which is the control that
//     stops the Howarth test from being passed by anything that reports
//     separation eagerly.
#include "Aeolion/Solver/AttachmentBoundaryLayer.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace S = Aeolion::Solver;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

// A synthetic run with a prescribed edge-velocity law. `psi` is carried
// along as the chordwise station so the reported separation psi can be
// checked too; here it is simply the arc length, which makes psi and s the
// same number and keeps the expected values readable.
S::SurfaceRun MakeRun(double length, int stations, const std::function<double(double)>& ue) {
    S::SurfaceRun run;
    for (int i = 1; i <= stations; ++i) {
        const double s = length * i / stations;
        run.S.push_back(s);
        run.Ue.push_back(ue(s));
        run.Psi.push_back(s);
        run.Zeta.push_back(0.0);
    }
    return run;
}

void TestFlatPlateRecoversBlasius() {
    // Re low enough that Michel does not trip before the far end, so the
    // whole run is the laminar branch being tested.
    const double Re = 2.0e4;
    const S::SurfaceRun run = MakeRun(1.0, 400, [](double) { return 1.0; });
    const S::SurfaceMarch march = S::MarchSurfaceRun(run, Re);

    CHECK(march.Valid, "a 400-station run must march");
    CHECK(!march.Separated(),
          "a flat plate has no adverse gradient and must not separate; got separation at psi=" +
              std::to_string(march.SeparationPsi));

    // theta(s) = sqrt(0.45 s / Re) for constant Ue. At s = 1 that is
    // 0.671/sqrt(Re), against Blasius' 0.664/sqrt(Re).
    const double thetaEnd = march.Theta.back();
    const double thwaites = std::sqrt(0.45 / Re);
    const double blasius = 0.664 / std::sqrt(Re);
    CHECK(std::fabs(thetaEnd - thwaites) < 0.01 * thwaites,
          "Thwaites on a flat plate must give theta = sqrt(0.45 s/Re), got " +
              std::to_string(thetaEnd) + " against " + std::to_string(thwaites));
    CHECK(std::fabs(thetaEnd - blasius) < 0.02 * blasius,
          "and must therefore sit within ~1% of Blasius, got " + std::to_string(thetaEnd) +
              " against " + std::to_string(blasius));

    // Zero pressure gradient means lambda = 0, so Thwaites' H is 2.61.
    CHECK(std::fabs(march.ShapeFactor.back() - 2.61) < 0.02,
          "a zero-gradient laminar layer has H = 2.61, got " +
              std::to_string(march.ShapeFactor.back()));
}

void TestHowarthFlowSeparatesAtTheKnownStation() {
    // Ue = 1 - s. Thwaites places separation at s = 0.123 (exact: 0.120).
    const double Re = 1.0e5;
    const S::SurfaceRun run = MakeRun(0.5, 2000, [](double s) { return 1.0 - s; });
    const S::SurfaceMarch march = S::MarchSurfaceRun(run, Re);

    CHECK(march.Valid, "the Howarth run must march");
    CHECK(march.BubbleFormed,
          "Howarth's retarded flow must reach laminar separation before Michel's criterion");
    CHECK(std::fabs(march.TransitionArc - 0.123) < 0.006,
          "Thwaites must place Howarth laminar separation at s/L = 0.123, got " +
              std::to_string(march.TransitionArc));
    // psi tracks s in this fixture, so the reported chordwise station must
    // agree with the arc length -- the interpolation must not drift.
    CHECK(std::fabs(march.TransitionPsi - march.TransitionArc) < 1e-6,
          "the interpolated bubble psi must track its arc length");

    // A layer that separates laminarly in a gradient this adverse cannot
    // then survive it turbulent either: the verdict must be separation.
    CHECK(march.Mode == S::SeparationMode::TurbulentSeparation,
          "a turbulent layer reattaching into Howarth's gradient must separate again");
    CHECK(march.SeparationArc > march.TransitionArc,
          "turbulent separation must lie downstream of the bubble that triggered it");

    // Refining the mesh must not move the answer: a crossing quantized to
    // stations would jump here.
    const S::SurfaceRun fine = MakeRun(0.5, 8000, [](double s) { return 1.0 - s; });
    const S::SurfaceMarch refined = S::MarchSurfaceRun(fine, Re);
    CHECK(std::fabs(refined.TransitionArc - march.TransitionArc) < 1e-3,
          "the bubble location must be mesh-converged: 2000 stations gave " +
              std::to_string(march.TransitionArc) + ", 8000 gave " +
              std::to_string(refined.TransitionArc));
}

void TestStagnationStartRecoversHiemenzTheta() {
    // Ue = a s near an attachment point. Thwaites' integral from zero must
    // produce theta_0 = sqrt(0.075/(Re a)) with no seed supplied.
    const double Re = 1.0e5;
    for (const double a : {0.5, 2.0, 8.0}) {
        const S::SurfaceRun run = MakeRun(0.02, 400, [a](double s) { return a * s; });
        const S::SurfaceMarch march = S::MarchSurfaceRun(run, Re);
        CHECK(march.Valid, "the stagnation run must march");

        const double expected = std::sqrt(0.075 / (Re * a));
        CHECK(std::fabs(march.StartingTheta - expected) < 0.02 * expected,
              "a march from a stagnation point must supply its own theta_0 = sqrt(0.075/(Re a)); at a=" +
                  std::to_string(a) + " got " + std::to_string(march.StartingTheta) + " against " +
                  std::to_string(expected));

        // And it must stay there: in a constant-strain region theta is
        // constant, which is the Hiemenz result.
        CHECK(std::fabs(march.Theta.back() - expected) < 0.03 * expected,
              "theta must stay at the Hiemenz value through a constant-strain region, got " +
                  std::to_string(march.Theta.back()));
    }
}

void TestFavourableGradientStaysAttached() {
    // The control for the Howarth test: an ACCELERATING flow over the same
    // length and Reynolds number must not separate.
    const double Re = 1.0e5;
    const S::SurfaceRun run = MakeRun(0.5, 2000, [](double s) { return 1.0 + s; });
    const S::SurfaceMarch march = S::MarchSurfaceRun(run, Re);

    CHECK(march.Valid, "the favourable-gradient run must march");
    CHECK(!march.BubbleFormed,
          "an accelerating flow cannot reach laminar separation; got a bubble at s=" +
              std::to_string(march.TransitionArc));
    CHECK(!march.Separated(),
          "and must reach the trailing edge attached; got separation at s=" +
              std::to_string(march.SeparationArc));
}

void TestEmptyAndDegenerateRunsAreDeclined() {
    CHECK(!S::MarchSurfaceRun(S::SurfaceRun{}, 1e5).Valid, "an empty run must be declined");

    const S::SurfaceRun tiny = MakeRun(1.0, 2, [](double) { return 1.0; });
    CHECK(!S::MarchSurfaceRun(tiny, 1e5).Valid,
          "a run below MinMarchStations cannot resolve a gradient and must be declined");

    const S::SurfaceRun run = MakeRun(1.0, 100, [](double) { return 1.0; });
    CHECK(!S::MarchSurfaceRun(run, 0.0).Valid, "a zero Reynolds number must be declined");
}

} // namespace

int main() {
    TestFlatPlateRecoversBlasius();
    TestHowarthFlowSeparatesAtTheKnownStation();
    TestStagnationStartRecoversHiemenzTheta();
    TestFavourableGradientStaysAttached();
    TestEmptyAndDegenerateRunsAreDeclined();

    if (failures == 0) std::cout << "PASS: TestAttachmentBoundaryLayer\n";
    return failures == 0 ? 0 : 1;
}
