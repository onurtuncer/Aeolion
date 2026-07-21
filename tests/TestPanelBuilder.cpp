// TestPanelBuilder.cpp -- validates that the lattice built from a geometry
// handoff is a genuine cambered, chordwise-discretized surface:
//
//   - Reduction: with one chordwise row and a symmetric (zero-camber)
//     section, the builder must reproduce the trusted
//     ToWingParams()+BuildWing() lattice essentially to floating point.
//     Two variants isolate twist and dihedral, because the two builders'
//     section frames agree exactly when only one of the pair is nonzero
//     (see PanelBuilder.h's "knowing deviations" note).
//   - Camber physics: a positively-cambered wing MUST lift at zero
//     geometric incidence, and a symmetric one must not. This is the whole
//     point of consuming the CST sections, and it is the check a flat
//     lattice structurally cannot pass.
//   - Chordwise discretization: requesting M rows must produce M panels per
//     spanwise strip, still report ONE station per strip, conserve total
//     planform area, and converge in CL as M increases.
//   - The on-disk JSON fixture must build and solve.
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Geometry/CstSurface.h"
#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/Solver/Solver.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using Aeolion::Geometry::AirfoilSection;
using Aeolion::Geometry::HandoffContract;
using Aeolion::Geometry::PlanformStation;
using Aeolion::Solver::Panel;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

namespace PB = Aeolion::PanelBuilder;

// Thin wrappers so the assertions below stay about aerodynamics rather than
// about builder plumbing.
std::vector<Panel> BuildLattice(const HandoffContract& contract) {
    return PB::LatticeBuilder(contract).Build();
}
std::vector<Panel> BuildLattice(const HandoffContract& contract, const PB::ControlDeflection& deflection) {
    return PB::LatticeBuilder(contract).Deflect(deflection).Build();
}
std::vector<Panel> BuildLattice(const HandoffContract& contract, const PB::LatticeOptions& options) {
    return PB::LatticeBuilder(contract, options).Build();
}

// A section whose upper and lower coefficients are exact negatives, so the
// mean line is identically zero: a symmetric airfoil, the only case that can
// be compared against the flat BuildWing() lattice.
AirfoilSection SymmetricSection(double eta) {
    AirfoilSection section;
    section.Eta = eta;
    section.CoefficientsUpper = {0.20, 0.24, 0.22, 0.18};
    section.CoefficientsLower = {-0.20, -0.24, -0.22, -0.18};
    return section;
}

// Upper and lower offset in the same direction, so the mean line bows
// upward: positive camber, which must lift at zero incidence.
AirfoilSection CamberedSection(double eta, double camberScale) {
    AirfoilSection section;
    section.Eta = eta;
    section.CoefficientsUpper = {0.20 + camberScale, 0.24 + camberScale, 0.22 + camberScale, 0.18 + camberScale};
    section.CoefficientsLower = {-0.20 + camberScale, -0.24 + camberScale, -0.22 + camberScale, -0.18 + camberScale};
    return section;
}

void CheckPanelsClose(const Panel& a, const Panel& b, double tol, const std::string& where) {
    CHECK(std::fabs(a.A.x - b.A.x) < tol && std::fabs(a.A.y - b.A.y) < tol && std::fabs(a.A.z - b.A.z) < tol,
          where + ": A mismatch");
    CHECK(std::fabs(a.B.x - b.B.x) < tol && std::fabs(a.B.y - b.B.y) < tol && std::fabs(a.B.z - b.B.z) < tol,
          where + ": B mismatch");
    CHECK(std::fabs(a.ControlPoint.x - b.ControlPoint.x) < tol &&
              std::fabs(a.ControlPoint.y - b.ControlPoint.y) < tol &&
              std::fabs(a.ControlPoint.z - b.ControlPoint.z) < tol,
          where + ": ControlPoint mismatch");
    CHECK(std::fabs(a.Normal.x - b.Normal.x) < tol && std::fabs(a.Normal.y - b.Normal.y) < tol &&
              std::fabs(a.Normal.z - b.Normal.z) < tol,
          where + ": Normal mismatch");
    CHECK(std::fabs(a.Area - b.Area) < tol, where + ": Area mismatch");
    CHECK(std::fabs(a.PlanformArea - b.PlanformArea) < tol, where + ": PlanformArea mismatch");
    CHECK(std::fabs(a.SpanwiseWidth - b.SpanwiseWidth) < tol, where + ": SpanwiseWidth mismatch");
}

// A single-trapezoid contract: exactly what ToWingParams() accepts.
HandoffContract MakeTrapezoidContract(double twistTipDeg, double dihedralDeg) {
    HandoffContract contract;
    contract.Span = 3.2;
    contract.Mesh.SpanwisePanelsPerSection = 14;
    contract.Mesh.ChordwisePanels = 1;

    PlanformStation root;
    root.Eta = 0.0; root.Chord = 0.5; root.TwistDeg = 0.0;
    root.SweepQuarterChordDeg = 12.0; root.DihedralDeg = dihedralDeg;

    PlanformStation tip;
    tip.Eta = 1.0; tip.Chord = 0.3; tip.TwistDeg = twistTipDeg;
    tip.SweepQuarterChordDeg = 12.0; tip.DihedralDeg = dihedralDeg;

    contract.Stations = {root, tip};
    contract.AirfoilSections = {SymmetricSection(0.0), SymmetricSection(1.0)};
    return contract;
}

// A rectangular wing, uniform section -- the cleanest case for isolating
// camber's effect from planform effects.
HandoffContract MakeRectangularContract(int chordwisePanels, double camberScale) {
    HandoffContract contract;
    contract.Span = 4.0;
    contract.Mesh.SpanwisePanelsPerSection = 12;
    contract.Mesh.ChordwisePanels = chordwisePanels;

    PlanformStation root;
    root.Eta = 0.0; root.Chord = 0.5; root.TwistDeg = 0.0;
    root.SweepQuarterChordDeg = 0.0; root.DihedralDeg = 0.0;
    PlanformStation tip = root;
    tip.Eta = 1.0;

    contract.Stations = {root, tip};
    contract.AirfoilSections = {CamberedSection(0.0, camberScale), CamberedSection(1.0, camberScale)};
    return contract;
}

// A cranked, three-station planform with sweep AND dihedral both changing
// between segments -- the case ToWingParams() explicitly refuses.
HandoffContract MakeCrankedContract() {
    HandoffContract contract;
    contract.Span = 4.0;
    contract.Mesh.SpanwisePanelsPerSection = 10;
    contract.Mesh.ChordwisePanels = 4;

    PlanformStation root;
    root.Eta = 0.0; root.Chord = 0.6; root.TwistDeg = 0.0;
    root.SweepQuarterChordDeg = 5.0; root.DihedralDeg = 0.0;

    PlanformStation kink;
    kink.Eta = 0.4; kink.Chord = 0.45; kink.TwistDeg = -1.0;
    kink.SweepQuarterChordDeg = 5.0; kink.DihedralDeg = 2.0;

    PlanformStation tip;
    tip.Eta = 1.0; tip.Chord = 0.2; tip.TwistDeg = -4.0;
    tip.SweepQuarterChordDeg = 25.0; tip.DihedralDeg = 6.0;

    contract.Stations = {root, kink, tip};
    contract.AirfoilSections = {CamberedSection(0.0, 0.05), CamberedSection(0.55, 0.03),
                                CamberedSection(1.0, 0.0)};
    return contract;
}

Aeolion::Solver::SolveResult SolveLattice(const std::vector<Panel>& panels, double span, double alphaDeg) {
    // Reference area is the planform projection, per the coefficient
    // convention (Solver::Panel).
    double area = 0.0;
    for (const auto& p : panels) area += p.PlanformArea;
    Aeolion::Solver::ReferenceGeometry ref;
    ref.Area = area; ref.Span = span; ref.Chord = area / span;
    Aeolion::Solver::FreestreamConditions fc;
    fc.Vinf = 25.0; fc.alphaDeg = alphaDeg; fc.rho = 1.225;
    return Aeolion::Solver::Solve(panels, fc, ref, 50.0 * span);
}

// --- CST evaluation, checked independently of any lattice ------------------
void TestCstCamberEvaluation() {
    const AirfoilSection symmetric = SymmetricSection(0.0);
    for (double psi : {0.1, 0.25, 0.5, 0.75, 0.9}) {
        CHECK(std::fabs(Aeolion::Geometry::SectionCamber(symmetric, psi)) < 1e-15,
              "symmetric section must have zero camber at psi=" + std::to_string(psi));
        CHECK(std::fabs(Aeolion::Geometry::SectionCamberSlope(symmetric, psi)) < 1e-12,
              "symmetric section must have zero camber slope at psi=" + std::to_string(psi));
    }

    const AirfoilSection cambered = CamberedSection(0.0, 0.08);
    CHECK(Aeolion::Geometry::SectionCamber(cambered, 0.5) > 0.0, "positive camber scale must bow the mean line up");
    // The sharp-TE fit drives both surfaces, and so the mean line, to zero
    // at the trailing edge (see CstSurface.h).
    CHECK(std::fabs(Aeolion::Geometry::SectionCamber(cambered, 1.0)) < 1e-9,
          "camber must vanish at the trailing edge");

    // Spanwise interpolation must land on the endpoints and stay bracketed.
    const std::vector<AirfoilSection> sections = {CamberedSection(0.0, 0.10), CamberedSection(1.0, 0.0)};
    const double atRoot = Aeolion::Geometry::CamberAt(sections, 0.0, 0.5);
    const double atTip = Aeolion::Geometry::CamberAt(sections, 1.0, 0.5);
    const double atMid = Aeolion::Geometry::CamberAt(sections, 0.5, 0.5);
    CHECK(std::fabs(atRoot - Aeolion::Geometry::SectionCamber(sections[0], 0.5)) < 1e-15,
          "eta=0 must evaluate the root section exactly");
    CHECK(std::fabs(atTip - Aeolion::Geometry::SectionCamber(sections[1], 0.5)) < 1e-15,
          "eta=1 must evaluate the tip section exactly");
    CHECK(atMid < atRoot && atMid > atTip, "mid-span camber must lie between the two sections");
}

// --- reduction to the trusted flat builder ---------------------------------
void CheckReducesToBuildWing(const HandoffContract& contract, const std::string& label) {
    const auto wing = Aeolion::Geometry::ToWingParams(contract);
    std::vector<Panel> expected = Aeolion::Solver::BuildWing(wing);
    for (auto& p : expected) p.Surface = "wing";

    // ToWingParams() leaves WingParams::CosineSpacing false, so the
    // comparison must ask for uniform explicitly -- the builder's own
    // default is cosine.
    const Aeolion::PanelBuilder::LatticeOptions uniform{Aeolion::PanelBuilder::Spacing::Uniform,
                                                        Aeolion::PanelBuilder::Spacing::Uniform};
    const std::vector<Panel> actual = BuildLattice(contract, uniform);

    CHECK(actual.size() == expected.size(),
          label + ": panel count " + std::to_string(actual.size()) + " != " + std::to_string(expected.size()));
    if (actual.size() != expected.size()) return;
    for (std::size_t i = 0; i < actual.size(); ++i)
        CheckPanelsClose(actual[i], expected[i], 1e-9, label + " panel[" + std::to_string(i) + "]");
}

void TestReducesToFlatBuilder() {
    // Twist only, then dihedral only: the two frames coincide exactly in
    // each case, so this is a strict equality check, not a tolerance dodge.
    CheckReducesToBuildWing(MakeTrapezoidContract(-3.5, 0.0), "twist-only");
    CheckReducesToBuildWing(MakeTrapezoidContract(0.0, 5.0), "dihedral-only");
}

// --- the actual point: camber generates lift at zero incidence -------------
void TestCamberLiftsAtZeroIncidence() {
    const HandoffContract symmetric = MakeRectangularContract(4, 0.0);
    const HandoffContract cambered = MakeRectangularContract(4, 0.06);
    const HandoffContract moreCambered = MakeRectangularContract(4, 0.12);

    const auto symmetricResult = SolveLattice(BuildLattice(symmetric), symmetric.Span, 0.0);
    const auto camberedResult = SolveLattice(BuildLattice(cambered), cambered.Span, 0.0);
    const auto moreResult = SolveLattice(BuildLattice(moreCambered), moreCambered.Span, 0.0);

    std::cout << "alpha=0: CL symmetric=" << symmetricResult.CL << "  cambered=" << camberedResult.CL
              << "  more cambered=" << moreResult.CL << "\n";

    CHECK(std::fabs(symmetricResult.CL) < 1e-6,
          "a symmetric section at zero incidence must produce ~zero lift, got " + std::to_string(symmetricResult.CL));
    CHECK(camberedResult.CL > 0.05,
          "a cambered section MUST lift at zero incidence -- a flat lattice cannot do this; got " +
              std::to_string(camberedResult.CL));
    CHECK(moreResult.CL > camberedResult.CL, "more camber must produce more lift at zero incidence");

    // Zero-lift angle must go negative for a positively-cambered section --
    // the classic thin-airfoil signature of camber.
    const auto camberedAtNegative =
        SolveLattice(BuildLattice(cambered), cambered.Span, -4.0);
    CHECK(camberedAtNegative.CL < camberedResult.CL, "CL must still increase with alpha for a cambered wing");
    CHECK(camberedAtNegative.CL < 0.0, "a moderately cambered wing should be below zero lift at alpha=-4deg");
}

// --- chordwise discretization ---------------------------------------------
void TestChordwiseRows() {
    const int rows = 5;
    const HandoffContract contract = MakeRectangularContract(rows, 0.04);
    const std::vector<Panel> panels = BuildLattice(contract);

    const std::size_t strips = 2 * static_cast<std::size_t>(contract.Mesh.SpanwisePanelsPerSection) *
                               (contract.Stations.size() - 1);
    CHECK(panels.size() == strips * static_cast<std::size_t>(rows),
          "expected " + std::to_string(strips * rows) + " panels, got " + std::to_string(panels.size()));

    // Every strip must hold exactly `rows` panels, and they must be distinct
    // chordwise positions (no duplicated row).
    std::map<int, int> perStrip;
    for (const auto& p : panels) {
        CHECK(p.StripIndex >= 0, "every panel from a handoff lattice must carry a StripIndex");
        ++perStrip[p.StripIndex];
    }
    CHECK(perStrip.size() == strips, "expected " + std::to_string(strips) + " strips, got " +
                                         std::to_string(perStrip.size()));
    for (const auto& [strip, count] : perStrip)
        CHECK(count == rows, "strip " + std::to_string(strip) + " has " + std::to_string(count) +
                                 " panels, expected " + std::to_string(rows));

    // The solver must collapse each chordwise stack into ONE station.
    const auto result = SolveLattice(panels, contract.Span, 3.0);
    CHECK(result.Stations.size() == strips,
          "solver should report one station per strip (" + std::to_string(strips) + "), got " +
              std::to_string(result.Stations.size()));
    CHECK(result.gamma.size() == panels.size(), "gamma must stay per-panel");

    // Planform area is conserved regardless of how finely the chord is cut.
    // (Summing Area instead would pick up the camber arc excess -- that is
    // the projection, and only the projection, being checked here.)
    double totalPlanform = 0.0;
    for (const auto& p : panels) totalPlanform += p.PlanformArea;
    const double expectedArea = contract.Stations.front().Chord * contract.Span;
    CHECK(std::fabs(totalPlanform - expectedArea) < 1e-9 * expectedArea,
          "total planform area " + std::to_string(totalPlanform) + " != " + std::to_string(expectedArea));

    // CL must converge as chordwise resolution increases, not wander.
    double previous = 0.0;
    for (int m : {1, 2, 4, 8}) {
        const HandoffContract refined = MakeRectangularContract(m, 0.04);
        const auto refinedResult =
            SolveLattice(BuildLattice(refined), refined.Span, 3.0);
        std::cout << "chordwise rows=" << m << "  CL=" << refinedResult.CL << "\n";
        CHECK(std::isfinite(refinedResult.CL), "CL must be finite at chordwise rows=" + std::to_string(m));
        if (m > 1)
            CHECK(std::fabs(refinedResult.CL - previous) < 0.05,
                  "CL jumped too much refining chordwise rows to " + std::to_string(m));
        previous = refinedResult.CL;
    }
}

void TestCrankedPlanform() {
    const HandoffContract contract = MakeCrankedContract();
    const std::vector<Panel> panels = BuildLattice(contract);

    double totalPlanform = 0.0;
    for (const auto& p : panels) totalPlanform += p.PlanformArea;

    // Independent area check: trapezoid per segment, doubled for both
    // semi-spans -- shares no code path with the builder's own chord math.
    double expectedHalfArea = 0.0;
    const double halfSpan = contract.Span * 0.5;
    for (std::size_t i = 0; i + 1 < contract.Stations.size(); ++i) {
        const auto& a = contract.Stations[i];
        const auto& b = contract.Stations[i + 1];
        expectedHalfArea += 0.5 * (a.Chord + b.Chord) * (b.Eta - a.Eta) * halfSpan;
    }
    const double expectedArea = 2.0 * expectedHalfArea;
    CHECK(std::fabs(totalPlanform - expectedArea) < 1e-9 * expectedArea,
          "cranked total planform area " + std::to_string(totalPlanform) + " != " + std::to_string(expectedArea));

    for (std::size_t i = 1; i < panels.size(); ++i)
        CHECK(panels[i].A.y >= panels[i - 1].A.y, "panels not sorted by A.y at " + std::to_string(i));

    const auto result = SolveLattice(panels, contract.Span, 4.0);
    CHECK(std::isfinite(result.CL) && std::isfinite(result.CDi), "cranked solve produced non-finite CL/CDi");
    CHECK(result.CL > 0.0, "cranked wing CL should be positive at alpha=4deg");
    CHECK(result.CDi > 0.0, "cranked wing CDi should be positive");
}

void TestJsonFixture() {
    const std::string path = std::string(AEOLION_TEST_DATA_DIR) + "/AeolionGeometryHandoff-1.0.0.json";
    std::ifstream input(path);
    CHECK(static_cast<bool>(input), "could not open fixture " + path);
    if (!input) return;

    nlohmann::json root;
    input >> root;
    const HandoffContract contract = Aeolion::Geometry::ParseHandoff(root);
    CHECK(contract.Mesh.ChordwisePanels > 1, "fixture should request more than one chordwise row");

    const std::vector<Panel> panels = BuildLattice(contract);
    const std::size_t strips = 2 * static_cast<std::size_t>(contract.Mesh.SpanwisePanelsPerSection) *
                               (contract.Stations.size() - 1);
    CHECK(panels.size() == strips * static_cast<std::size_t>(contract.Mesh.ChordwisePanels),
          "fixture panel count mismatch");

    double totalPlanform = 0.0;
    double totalSurface = 0.0;
    for (const auto& p : panels) { totalPlanform += p.PlanformArea; totalSurface += p.Area; }
    const double expectedArea = contract.Stations.front().Chord * contract.Span; // rectangular planform
    CHECK(std::fabs(totalPlanform - expectedArea) < 1e-6 * expectedArea, "fixture planform area mismatch");
    // The fixture's sections are genuinely cambered, so its wetted surface
    // must come out larger than its projection.
    CHECK(totalSurface > totalPlanform, "fixture cambered surface area should exceed its planform area");

    const auto result = SolveLattice(panels, contract.Span, 5.0);
    std::cout << "fixture: panels=" << panels.size() << "  stations=" << result.Stations.size()
              << "  CL=" << result.CL << "  CDi=" << result.CDi << "\n";
    CHECK(std::isfinite(result.CL) && result.CL > 0.0, "fixture CL should be positive and finite at alpha=5deg");
}

// --- true surface area vs planform projection ------------------------------
void TestCamberedAreaExceedsPlanform() {
    // Flat, no dihedral: the curved surface IS its projection, exactly.
    const HandoffContract flat = MakeRectangularContract(6, 0.0);
    for (const auto& p : BuildLattice(flat))
        CHECK(std::fabs(p.Area - p.PlanformArea) < 1e-12 * p.PlanformArea,
              "a flat un-dihedralled panel's area must equal its planform projection");

    // Cambered: every panel's true surface must be strictly longer.
    const HandoffContract cambered = MakeRectangularContract(6, 0.10);
    double totalArea = 0.0;
    double totalPlanform = 0.0;
    for (const auto& p : BuildLattice(cambered)) {
        CHECK(p.Area > p.PlanformArea, "a cambered panel's arc area must exceed its projection");
        totalArea += p.Area;
        totalPlanform += p.PlanformArea;
    }
    std::cout << "cambered surface area=" << totalArea << "  planform=" << totalPlanform
              << "  ratio=" << totalArea / totalPlanform << "\n";
    CHECK(totalArea > totalPlanform, "total surface area must exceed total planform area");

    // Planform area is unchanged by camber -- it is a projection.
    const double expectedPlanform = cambered.Stations.front().Chord * cambered.Span;
    CHECK(std::fabs(totalPlanform - expectedPlanform) < 1e-9 * expectedPlanform,
          "camber must not change the planform area");

    // Dihedral tilts the surface out of the reference plane by exactly
    // 1/cos(dihedral), independent of camber.
    const HandoffContract dihedral = MakeTrapezoidContract(0.0, 20.0);
    for (const auto& p : BuildLattice(dihedral)) {
        const double expected = p.PlanformArea / std::cos(Aeolion::Math::DegToRad(20.0));
        CHECK(std::fabs(p.Area - expected) < 1e-9 * expected,
              "dihedral must scale a flat panel's area by 1/cos(dihedral)");
    }
}

// --- control surface deflection --------------------------------------------
// A trailing-edge flap over the outer span, hinged about the spanwise axis.
HandoffContract MakeFlappedContract(Aeolion::Geometry::ControlSurfaceBinding binding) {
    HandoffContract contract = MakeRectangularContract(6, 0.0); // symmetric section: isolate the flap
    Aeolion::Geometry::ControlSurface flap;
    flap.Name = "aileron";
    flap.ChordFraction = 0.25;
    flap.EtaStart = 0.4;
    flap.EtaEnd = 1.0;
    flap.HingeAxis = Aeolion::Math::Vec3(0.0, 1.0, 0.0);
    flap.Binding = binding;
    contract.ControlSurfaces = {flap};
    return contract;
}

void TestControlSurfaceDeflection() {
    using Aeolion::Geometry::ControlSurfaceBinding;
    const HandoffContract contract = MakeFlappedContract(ControlSurfaceBinding::Wing);

    const auto neutral = BuildLattice(contract);
    const auto flapDown =
        BuildLattice(contract, Aeolion::PanelBuilder::Symmetric(0, 10.0));
    const auto flapUp =
        BuildLattice(contract, Aeolion::PanelBuilder::Symmetric(0, -10.0));
    const auto aileron =
        BuildLattice(contract, Aeolion::PanelBuilder::Antisymmetric(0, 10.0));

    // Deflection must not change the lattice topology -- only its shape.
    CHECK(flapDown.size() == neutral.size(), "deflection must not change panel count");

    const auto neutralResult = SolveLattice(neutral, contract.Span, 0.0);
    const auto downResult = SolveLattice(flapDown, contract.Span, 0.0);
    const auto upResult = SolveLattice(flapUp, contract.Span, 0.0);
    const auto aileronResult = SolveLattice(aileron, contract.Span, 0.0);

    std::cout << "flap: neutral CL=" << neutralResult.CL << "  +10deg CL=" << downResult.CL
              << "  -10deg CL=" << upResult.CL << "\n"
              << "aileron +/-10deg: CL=" << aileronResult.CL << "  Croll=" << aileronResult.Croll << "\n";

    // Symmetric section, zero incidence, no deflection: no lift.
    CHECK(std::fabs(neutralResult.CL) < 1e-6, "undeflected symmetric wing must not lift at alpha=0");

    // Positive (trailing-edge down) deflection must generate lift, and the
    // response must be odd in the deflection angle.
    CHECK(downResult.CL > 0.1, "trailing-edge-down flap must generate lift, got " + std::to_string(downResult.CL));
    CHECK(upResult.CL < -0.1, "trailing-edge-up flap must generate downforce, got " + std::to_string(upResult.CL));
    CHECK(std::fabs(downResult.CL + upResult.CL) < 1e-6, "flap lift response must be odd in deflection angle");

    // A flap deflection is symmetric: it must NOT roll the aircraft.
    CHECK(std::fabs(downResult.Croll) < 1e-6, "a symmetric flap deflection must not produce rolling moment");

    // An aileron deflection is antisymmetric: it must roll and must not lift.
    CHECK(std::fabs(aileronResult.CL) < 1e-6, "an antisymmetric aileron deflection must not change total lift");
    CHECK(std::fabs(aileronResult.Croll) > 1e-3,
          "an aileron deflection MUST produce a rolling moment, got " + std::to_string(aileronResult.Croll));

    // Rigid rotation preserves arc length, so the true surface area is
    // deflection-invariant even though the projection is not.
    double neutralArea = 0.0, deflectedArea = 0.0;
    for (const auto& p : neutral) neutralArea += p.Area;
    for (const auto& p : flapDown) deflectedArea += p.Area;
    CHECK(std::fabs(neutralArea - deflectedArea) < 1e-9 * neutralArea,
          "a rigid hinge rotation must preserve total surface area");
}

void TestDuctJetSurfacesAreNeverLatticed() {
    using Aeolion::Geometry::ControlSurfaceBinding;
    // The identical surface, bound to the duct jet instead of the wing:
    // ControlSurface.h is emphatic that it must not enter the wing lattice.
    const HandoffContract ductJet = MakeFlappedContract(ControlSurfaceBinding::DuctJet);

    const auto neutral = BuildLattice(ductJet);
    const auto commanded =
        BuildLattice(ductJet, Aeolion::PanelBuilder::Symmetric(0, 20.0));

    CHECK(neutral.size() == commanded.size(), "duct-jet surface must not change the wing lattice size");
    for (std::size_t i = 0; i < neutral.size(); ++i)
        CheckPanelsClose(commanded[i], neutral[i], 1e-12,
                         "duct-jet deflection must leave wing panel " + std::to_string(i) + " untouched");

    const auto result = SolveLattice(commanded, ductJet.Span, 0.0);
    CHECK(std::fabs(result.CL) < 1e-6, "commanding a duct-jet vane must not lift the wing lattice");
}

// The hinge line must fall on a panel edge, not inside a panel, so the
// control surface's chord extent is exactly what the contract asked for.
void TestHingeLandsOnPanelEdge() {
    using Aeolion::Geometry::ControlSurfaceBinding;
    const HandoffContract contract = MakeFlappedContract(ControlSurfaceBinding::Wing);
    const double hingeChordFraction = contract.ControlSurfaces.front().ChordFraction;

    const auto neutral = BuildLattice(contract);
    const auto deflected =
        BuildLattice(contract, Aeolion::PanelBuilder::Symmetric(0, 15.0));

    // Panels that moved are exactly the aft ones; their planform areas must
    // sum to the commanded chord fraction of the covered strips' area.
    double movedPlanform = 0.0;
    double coveredPlanform = 0.0;
    const double halfSpan = contract.Span * 0.5;
    for (std::size_t i = 0; i < neutral.size(); ++i) {
        const double stripEta = std::fabs(0.5 * (neutral[i].A.y + neutral[i].B.y)) / halfSpan;
        const bool covered = stripEta >= contract.ControlSurfaces.front().EtaStart &&
                             stripEta <= contract.ControlSurfaces.front().EtaEnd;
        if (!covered) continue;
        coveredPlanform += neutral[i].PlanformArea;

        const double motion = (deflected[i].ControlPoint - neutral[i].ControlPoint).Norm();
        if (motion > 1e-9) movedPlanform += neutral[i].PlanformArea;
    }
    CHECK(coveredPlanform > 0.0, "the flap band should cover some strips");
    const double movedFraction = movedPlanform / coveredPlanform;
    std::cout << "hinge: moved chord fraction=" << movedFraction << " (requested " << hingeChordFraction << ")\n";
    CHECK(std::fabs(movedFraction - hingeChordFraction) < 1e-9,
          "the deflected portion must be exactly the requested chord fraction, got " +
              std::to_string(movedFraction));
}

// --- the mesh must land on control surface edges ---------------------------
// Without a breakpoint there, a surface's span snaps to whatever panel
// boundary happens to be nearest, making its control authority a property
// of the mesh rather than of the contract.
void TestControlSurfaceBandIsExact() {
    using Aeolion::Geometry::ControlSurfaceBinding;
    HandoffContract contract = MakeRectangularContract(6, 0.0);
    // Band edges deliberately chosen NOT to fall on any boundary of the
    // underlying 12-panel division, so the mesh has to insert them.
    Aeolion::Geometry::ControlSurface flap;
    flap.Name = "aileron";
    flap.ChordFraction = 0.25;
    flap.EtaStart = 0.37;
    flap.EtaEnd = 0.91;
    flap.HingeAxis = Aeolion::Math::Vec3(0.0, 1.0, 0.0);
    flap.Binding = ControlSurfaceBinding::Wing;
    contract.ControlSurfaces = {flap};

    const auto neutral = BuildLattice(contract);
    const auto moved =
        BuildLattice(contract, Aeolion::PanelBuilder::Antisymmetric(0, 10.0));

    const double halfSpan = contract.Span * 0.5;
    double lo = 2.0, hi = -1.0;
    for (std::size_t i = 0; i < neutral.size(); ++i) {
        if ((moved[i].ControlPoint - neutral[i].ControlPoint).Norm() < 1e-12) continue;
        const double etaA = std::fabs(neutral[i].A.y) / halfSpan;
        const double etaB = std::fabs(neutral[i].B.y) / halfSpan;
        lo = std::min({lo, etaA, etaB});
        hi = std::max({hi, etaA, etaB});
    }
    std::cout << "band: requested [" << flap.EtaStart << ", " << flap.EtaEnd << "]  meshed [" << lo << ", " << hi
              << "]\n";
    CHECK(std::fabs(lo - flap.EtaStart) < 1e-9,
          "the deflected band's inboard edge must sit exactly at eta_start, got " + std::to_string(lo));
    CHECK(std::fabs(hi - flap.EtaEnd) < 1e-9,
          "the deflected band's outboard edge must sit exactly at eta_end, got " + std::to_string(hi));

    // Subdividing at the band edges must redistribute the section's panel
    // budget, not inflate it.
    const auto plain = BuildLattice(MakeRectangularContract(6, 0.0));
    CHECK(neutral.size() == plain.size(),
          "adding a control surface must not change panel count (" + std::to_string(neutral.size()) + " vs " +
              std::to_string(plain.size()) + ")");
}

// Cosine must actually cluster, and must leave the wing's total geometry
// alone -- it redistributes panels, it does not add or remove area.
void TestCosineSpacingClusters() {
    using Aeolion::PanelBuilder::LatticeOptions;
    using Aeolion::PanelBuilder::Spacing;
    const HandoffContract contract = MakeRectangularContract(1, 0.0);

    const auto uniform = BuildLattice(contract, {Spacing::Uniform, Spacing::Uniform});
    const auto cosine = BuildLattice(contract, {Spacing::Cosine, Spacing::Uniform});
    CHECK(uniform.size() == cosine.size(), "spacing must not change panel count");

    const auto widthSpread = [](const std::vector<Panel>& panels) {
        double lo = 1e30, hi = -1e30;
        for (const auto& p : panels) { lo = std::min(lo, p.SpanwiseWidth); hi = std::max(hi, p.SpanwiseWidth); }
        return hi / lo;
    };
    CHECK(std::fabs(widthSpread(uniform) - 1.0) < 1e-9, "uniform spacing must give equal-width panels");
    CHECK(widthSpread(cosine) > 3.0,
          "cosine spacing should cluster strongly at the tips, width ratio was " +
              std::to_string(widthSpread(cosine)));

    double uniformArea = 0.0, cosineArea = 0.0;
    for (const auto& p : uniform) uniformArea += p.PlanformArea;
    for (const auto& p : cosine) cosineArea += p.PlanformArea;
    CHECK(std::fabs(uniformArea - cosineArea) < 1e-9 * uniformArea,
          "spacing must not change total planform area");
}

} // namespace

int main() {
    TestCstCamberEvaluation();
    TestReducesToFlatBuilder();
    TestCamberLiftsAtZeroIncidence();
    TestChordwiseRows();
    TestCrankedPlanform();
    TestCamberedAreaExceedsPlanform();
    TestControlSurfaceDeflection();
    TestDuctJetSurfacesAreNeverLatticed();
    TestHingeLandsOnPanelEdge();
    TestControlSurfaceBandIsExact();
    TestCosineSpacingClusters();
    TestJsonFixture();

    if (failures == 0) { std::cout << "PASS: TestPanelBuilder\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestPanelBuilder\n";
    return 1;
}
