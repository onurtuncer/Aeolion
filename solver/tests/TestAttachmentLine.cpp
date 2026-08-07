// TestAttachmentLine.cpp -- validates the wing attachment line, and in
// particular that it responds to sideslip the way a swept wing physically
// must.
//
//   - Zero incidence. A symmetric section on a flat lattice carries no
//     circulation at alpha = 0, so there is no induced flow and the
//     attachment line must sit exactly on the leading edge at every strip.
//
//   - Positive incidence. The attachment line moves onto the LOWER surface
//     at every strip, and its 3-D position must lie on the section surface
//     just aft of and below the leading edge.
//
//   - Geometric sweep. With no sideslip the effective sweep -- the angle
//     between the local flow and the leading-edge-normal plane -- must
//     recover the wing's geometric sweep, and must be the same magnitude on
//     both wings.
//
//   - Sideslip, the case the whole leading-edge-normal formulation exists
//     for. On a SWEPT wing at sideslip the effective sweep rises on one
//     wing and falls on the other, so the attachment-line Reynolds number
//     -- and hence the wing's susceptibility to leading-edge contamination
//     -- becomes asymmetric. A streamwise formulation returns a symmetric
//     answer here, so this is the assertion that distinguishes the two.
//
//   - Unswept control. The same sideslip on an UNSWEPT wing must NOT produce
//     that asymmetry: the flow runs along the span equally on both sides.
//     Without this, the sideslip test above would also pass for code that
//     merely keyed off the sign of y.
#include "Aeolion/Solver/AttachmentLine.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using Aeolion::Geometry::AirfoilSection;
using Aeolion::Math::Vec3;
namespace S = Aeolion::Solver;
namespace M = Aeolion::Math;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

struct Wing {
    std::vector<S::Panel> Panels;
    std::vector<S::StripSection> Strips;
    double Span = 0.0, Chord = 0.0, SweepRad = 0.0;
};

// A flat, untwisted, constant-chord wing with a swept leading edge: one
// Weissinger row per strip, ordered left tip to right tip.
Wing BuildSweptWing(double span, double chord, double sweepDeg, int strips) {
    Wing wing;
    wing.Span = span;
    wing.Chord = chord;
    wing.SweepRad = M::DegToRad(sweepDeg);
    const double tanSweep = std::tan(wing.SweepRad);
    const double halfSpan = 0.5 * span;

    for (int i = 0; i < strips; ++i) {
        const double yA = -halfSpan + span * i / strips;
        const double yB = -halfSpan + span * (i + 1) / strips;
        const double yC = 0.5 * (yA + yB);

        // Leading edge sweeps aft with |y|; the bound vortex sits a quarter
        // chord behind it.
        const auto leadingEdgeX = [&](double y) { return std::fabs(y) * tanSweep; };

        S::Panel panel;
        panel.A = Vec3(leadingEdgeX(yA) + 0.25 * chord, yA, 0.0);
        panel.B = Vec3(leadingEdgeX(yB) + 0.25 * chord, yB, 0.0);
        panel.ControlPoint = Vec3(leadingEdgeX(yC) + 0.75 * chord, yC, 0.0);
        panel.Normal = Vec3(0, 0, 1);
        panel.TrailDirA = Vec3(1, 0, 0);
        panel.TrailDirB = Vec3(1, 0, 0);
        panel.SpanwiseWidth = yB - yA;
        panel.PlanformArea = chord * (yB - yA);
        panel.Area = panel.PlanformArea;
        panel.Surface = "wing";
        wing.Panels.push_back(panel);

        S::StripSection strip;
        strip.ChordDir = Vec3(1, 0, 0);
        strip.LiftDir = Vec3(0, 0, 1);
        strip.Chord = chord;
        strip.Width = yB - yA;
        strip.Eta = std::fabs(yC) / halfSpan;
        wing.Strips.push_back(strip);
    }
    return wing;
}

std::vector<AirfoilSection> SymmetricSections(double scale) {
    AirfoilSection section;
    section.Eta = 0.0;
    section.CoefficientsUpper = {scale, scale, scale, scale};
    section.CoefficientsLower = {-scale, -scale, -scale, -scale};
    return {section};
}

// Solve a wing and hand back its attachment line. The PreparedSystem is
// created inside, so the FlowField cannot outlive it -- everything the
// caller needs is extracted before returning.
S::AttachmentLine SolveWing(const Wing& wing, const std::vector<AirfoilSection>& sections,
                            double alphaDeg, double betaDeg) {
    const S::PanelSystem system{wing.Panels, {}};
    const S::PreparedSystem prepared = S::Prepare(system, 50.0 * wing.Span);

    S::FreestreamConditions fc;
    fc.Vinf = 40.0;
    fc.alphaDeg = alphaDeg;
    fc.betaDeg = betaDeg;

    S::ReferenceGeometry ref;
    ref.Area = wing.Span * wing.Chord;
    ref.Span = wing.Span;
    ref.Chord = wing.Chord;

    const S::SolveResult result = S::SolveWithSystem(prepared, fc, ref);
    const S::FlowField field = S::MakeFlowField(prepared, fc, result.gamma, result.sigma);
    return S::ComputeAttachmentLine(field, wing.Panels, wing.Strips, sections);
}

// --------------------------------------------------------------- tests ----

void TestZeroIncidenceSitsOnTheLeadingEdge() {
    const Wing wing = BuildSweptWing(8.0, 1.0, 0.0, 24);
    const S::AttachmentLine line = SolveWing(wing, SymmetricSections(0.17), 0.0, 0.0);

    CHECK(line.Stations.size() == wing.Panels.size(), "one attachment station per strip");
    int found = 0;
    for (const S::AttachmentStation& station : line.Stations) {
        if (!station.Found) continue;
        ++found;
        CHECK(std::fabs(station.StagnationOffset) < 1e-3,
              "a symmetric section at zero incidence attaches ON the leading edge, got offset " +
                  std::to_string(station.StagnationOffset));
        CHECK(std::fabs(station.AlphaNormalDeg) < 1e-6,
              "an uncirculated flat lattice induces nothing at alpha=0, got alpha_n=" +
                  std::to_string(station.AlphaNormalDeg));
    }
    CHECK(found == static_cast<int>(line.Stations.size()), "every strip must resolve an attachment point");
}

void TestPositiveIncidenceMovesOntoTheLowerSurface() {
    const Wing wing = BuildSweptWing(8.0, 1.0, 0.0, 24);
    const S::AttachmentLine line = SolveWing(wing, SymmetricSections(0.17), 6.0, 0.0);

    for (const S::AttachmentStation& station : line.Stations) {
        CHECK(station.Found, "every strip must resolve an attachment point at 6 deg");
        if (!station.Found) continue;

        CHECK(station.StagnationOffset > 0.0,
              "at positive incidence the attachment line lies on the LOWER surface, got offset " +
                  std::to_string(station.StagnationOffset));
        CHECK(station.AttachmentPoint.z < 0.0,
              "the attachment point must sit below the chord plane, got z=" +
                  std::to_string(station.AttachmentPoint.z));
        CHECK(station.AttachmentPoint.x > station.LeadingEdge.x,
              "the attachment point must sit aft of the leading edge");

        // Induced downwash reduces the local incidence below the geometric
        // 6 degrees, but never reverses it on a wing at positive alpha.
        CHECK(station.AlphaNormalDeg > 0.0 && station.AlphaNormalDeg < 6.0,
              "local incidence must be positive but below geometric, got " +
                  std::to_string(station.AlphaNormalDeg));

        CHECK(station.StrainRate > 0.0, "the strain rate must be positive");
        CHECK(station.MomentumThickness > 0.0 && station.MomentumThickness < 1e-3,
              "the starting momentum thickness must be small but nonzero, got " +
                  std::to_string(station.MomentumThickness));
    }
}

void TestEffectiveSweepRecoversGeometricSweep() {
    const double sweepDeg = 30.0;
    const Wing wing = BuildSweptWing(8.0, 1.0, sweepDeg, 24);
    const S::AttachmentLine line = SolveWing(wing, SymmetricSections(0.17), 4.0, 0.0);

    // Skip the tip strips, whose one-sided leading-edge difference and strong
    // induced flow both blur the comparison. The ROOT strips are not skipped:
    // a swept wing's leading edge genuinely breaks at the centreline, and the
    // point of detecting that is to get the sweep right there anyway.
    int kinked = 0;
    for (std::size_t i = 2; i + 2 < line.Stations.size(); ++i) {
        const S::AttachmentStation& station = line.Stations[i];
        if (!station.Found) continue;
        if (station.AtKink) ++kinked;
        const double effective = M::RadToDeg(station.EffectiveSweepRad);
        CHECK(std::fabs(effective - sweepDeg) < 3.0,
              "with no sideslip the effective sweep must recover the geometric sweep, got " +
                  std::to_string(effective) + " deg at strip " + std::to_string(i) +
                  (station.AtKink ? " (at kink)" : ""));
    }
    // Exactly the two strips flanking the centreline see the root break.
    CHECK(kinked == 2, "a swept wing's root kink must be detected at exactly two strips, got " +
                           std::to_string(kinked));

    // An unswept wing has no break to find.
    const Wing straight = BuildSweptWing(8.0, 1.0, 0.0, 24);
    const S::AttachmentLine straightLine = SolveWing(straight, SymmetricSections(0.17), 4.0, 0.0);
    for (const S::AttachmentStation& station : straightLine.Stations)
        CHECK(!station.AtKink, "a straight leading edge must not be reported as kinked");

    // And it must be symmetric between the two wings.
    const std::size_t last = line.Stations.size() - 1;
    for (std::size_t i = 2; i + 2 < line.Stations.size() / 2; ++i) {
        const S::AttachmentStation& left = line.Stations[i];
        const S::AttachmentStation& right = line.Stations[last - i];
        if (!left.Found || !right.Found) continue;
        CHECK(std::fabs(left.EffectiveSweepRad - right.EffectiveSweepRad) < M::DegToRad(1.0),
              "at zero sideslip the two wings must see the same effective sweep");
        CHECK(std::fabs(left.AttachmentLineReynolds - right.AttachmentLineReynolds) <
                  0.05 * std::max(left.AttachmentLineReynolds, 1.0),
              "at zero sideslip the attachment-line Reynolds number must be symmetric");
    }
}

// The load-bearing test: sideslip must break the left/right symmetry of a
// SWEPT wing's attachment line.
void TestSideslipBreaksSweptWingSymmetry() {
    const double sweepDeg = 30.0, betaDeg = 10.0;
    const Wing wing = BuildSweptWing(8.0, 1.0, sweepDeg, 24);
    const S::AttachmentLine line = SolveWing(wing, SymmetricSections(0.17), 2.0, betaDeg);

    const std::size_t mid = line.Stations.size() / 2;
    const S::AttachmentStation& left = line.Stations[mid / 2];              // left wing
    const S::AttachmentStation& right = line.Stations[mid + mid / 2];       // right wing
    CHECK(left.Found && right.Found, "both wings must resolve an attachment point at sideslip");
    if (!left.Found || !right.Found) return;

    // Sideslip toward +y adds to the right wing's spanwise flow (its leading
    // edge runs outboard AND aft) and subtracts from the left's.
    const double leftSweep = M::RadToDeg(left.EffectiveSweepRad);
    const double rightSweep = M::RadToDeg(right.EffectiveSweepRad);
    CHECK(rightSweep > leftSweep + 8.0,
          "positive sideslip must sweep the RIGHT wing further from the flow: left " +
              std::to_string(leftSweep) + " deg vs right " + std::to_string(rightSweep) + " deg");

    // Roughly sweep +/- beta, which is the infinite-swept-wing reading.
    CHECK(std::fabs(rightSweep - (sweepDeg + betaDeg)) < 4.0,
          "the windward wing's effective sweep should be near sweep + beta, got " +
              std::to_string(rightSweep));
    CHECK(std::fabs(leftSweep - (sweepDeg - betaDeg)) < 4.0,
          "the leeward wing's effective sweep should be near sweep - beta, got " +
              std::to_string(leftSweep));

    // More spanwise flow means a higher attachment-line Reynolds number, so
    // the two wings are NOT equally close to leading-edge contamination.
    CHECK(right.AttachmentLineReynolds > 1.3 * left.AttachmentLineReynolds,
          "the more swept wing must carry the higher attachment-line Reynolds number: left " +
              std::to_string(left.AttachmentLineReynolds) + " vs right " +
              std::to_string(right.AttachmentLineReynolds));

    // Reversing sideslip must mirror the asymmetry exactly.
    const S::AttachmentLine mirrored = SolveWing(wing, SymmetricSections(0.17), 2.0, -betaDeg);
    const S::AttachmentStation& mirroredLeft = mirrored.Stations[mid / 2];
    const S::AttachmentStation& mirroredRight = mirrored.Stations[mid + mid / 2];
    CHECK(std::fabs(mirroredLeft.EffectiveSweepRad - right.EffectiveSweepRad) < M::DegToRad(1.0),
          "reversing sideslip must swap the two wings' effective sweeps");
    CHECK(std::fabs(mirroredRight.EffectiveSweepRad - left.EffectiveSweepRad) < M::DegToRad(1.0),
          "reversing sideslip must swap the two wings' effective sweeps");
}

// The control: an UNSWEPT wing at the same sideslip must stay symmetric.
void TestSideslipLeavesUnsweptWingSymmetric() {
    const double betaDeg = 10.0;
    const Wing wing = BuildSweptWing(8.0, 1.0, 0.0, 24);
    const S::AttachmentLine line = SolveWing(wing, SymmetricSections(0.17), 2.0, betaDeg);

    const std::size_t mid = line.Stations.size() / 2;
    const S::AttachmentStation& left = line.Stations[mid / 2];
    const S::AttachmentStation& right = line.Stations[mid + mid / 2];
    CHECK(left.Found && right.Found, "both wings must resolve an attachment point");
    if (!left.Found || !right.Found) return;

    CHECK(std::fabs(left.EffectiveSweepRad - right.EffectiveSweepRad) < M::DegToRad(1.0),
          "an unswept wing sees the same spanwise flow on both sides at sideslip: left " +
              std::to_string(M::RadToDeg(left.EffectiveSweepRad)) + " vs right " +
              std::to_string(M::RadToDeg(right.EffectiveSweepRad)));

    // And that shared effective sweep is just the sideslip angle itself.
    CHECK(std::fabs(M::RadToDeg(left.EffectiveSweepRad) - betaDeg) < 2.0,
          "on an unswept wing the effective sweep IS the sideslip, got " +
              std::to_string(M::RadToDeg(left.EffectiveSweepRad)));
}

void TestAttachmentLineStateThresholds() {
    // Rbar rises with speed and with sweep, so a fast, strongly swept wing
    // should cross into contamination while a slow unswept one does not.
    const Wing swept = BuildSweptWing(12.0, 1.5, 45.0, 24);
    const S::AttachmentLine hot = SolveWing(swept, SymmetricSections(0.17), 3.0, 0.0);

    double worst = 0.0;
    for (const S::AttachmentStation& station : hot.Stations)
        if (station.Found) worst = std::max(worst, station.AttachmentLineReynolds);
    CHECK(worst > 0.0, "a swept wing must produce a nonzero attachment-line Reynolds number");

    // Whatever the number, the reported state must agree with the thresholds.
    for (const S::AttachmentStation& station : hot.Stations) {
        if (!station.Found) continue;
        const double rbar = station.AttachmentLineReynolds;
        const S::AttachmentLineState expected =
            (rbar >= S::AttachmentLineTransitionReynolds)      ? S::AttachmentLineState::Turbulent
            : (rbar >= S::AttachmentLineContaminationReynolds) ? S::AttachmentLineState::Contaminated
                                                               : S::AttachmentLineState::Laminar;
        CHECK(station.State == expected, "the reported state must follow the thresholds at Rbar=" +
                                             std::to_string(rbar));
    }
}

void TestNoSectionDataYieldsNoAttachmentPoint() {
    const Wing wing = BuildSweptWing(8.0, 1.0, 0.0, 12);
    const S::AttachmentLine line = SolveWing(wing, {}, 5.0, 0.0);
    for (const S::AttachmentStation& station : line.Stations)
        CHECK(!station.Found,
              "with no section data there is no thickness, so no stagnation point may be claimed");
}

} // namespace

int main() {
    TestZeroIncidenceSitsOnTheLeadingEdge();
    TestPositiveIncidenceMovesOntoTheLowerSurface();
    TestEffectiveSweepRecoversGeometricSweep();
    TestSideslipBreaksSweptWingSymmetry();
    TestSideslipLeavesUnsweptWingSymmetric();
    TestAttachmentLineStateThresholds();
    TestNoSectionDataYieldsNoAttachmentPoint();

    if (failures == 0) std::cout << "PASS: TestAttachmentLine\n";
    return failures == 0 ? 0 : 1;
}
