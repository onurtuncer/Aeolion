// PanelBuilder/PanelBuilder.h
//
// LatticeBuilder turns a parsed geometry handoff into the solver's lattice
// -- station by station, on the CAMBER SURFACE the handoff's CST sections
// describe, discretized both spanwise and chordwise, with wing-bound
// control surfaces deflected about their stated hinge lines. This is the
// "station-resolved lattice builder" that HandoffContract.h's own comment
// calls out as what supersedes Geometry::ToWingParams(), which can only
// reduce a design to Solver::WingParams's single flat trapezoid.
//
// Construction does the geometry-only work (breakpoints, spanwise station
// march) once; deflections are applied per Build(), so sweeping a control
// surface through many angles does not redo the planform march:
//
//     LatticeBuilder builder(contract);
//     auto neutral = builder.Build();
//     auto rolled  = builder.Deflect(Antisymmetric(aileron, 10.0)).Build();
//
// --- planform: how the quarter-chord line is traced ------------------------
// Geometry::PlanformStation.h states the rule this relies on: Eta, Chord,
// TwistDeg, SweepQuarterChordDeg, and DihedralDeg are all linearly
// interpolated between adjacent stations. Chord and twist are therefore
// exact by direct linear interpolation. Sweep and dihedral are angles, and
// a linearly-varying angle has no closed-form quarter-chord curve (the
// curve integrates tan(angle) over span), so the quarter-chord position is
// accumulated outward from the root one panel boundary at a time using each
// sub-interval's midpoint angle. That accumulation is EXACT whenever
// sweep/dihedral are constant across a station segment (the tangent is then
// genuinely constant and the sum telescopes to the closed form, regardless
// of spanwise_panels_per_section) and converges with panel count otherwise.
//
// --- camber: what makes this not a flat lattice ----------------------------
// Every lattice point is placed on the mean line of the local CST section
// (Geometry/CstSurface.h), interpolated across span between the handoff's
// airfoil_sections. The panel's flow-tangency normal is taken from the
// camber SLOPE at its own three-quarter-chord point, not merely from twist,
// so a cambered section develops lift at zero geometric incidence -- the
// physical effect a flat lattice structurally cannot represent.
//
// Panel::Area is the area of that curved surface, from the true camber-line
// arc length and the dihedral-tilted spanwise extent. Panel::PlanformArea
// is its projection, which is what every force coefficient normalizes by
// (see Solver::Panel for why the two must not be confused).
//
// --- chordwise discretization ----------------------------------------------
// mesh_topology.chordwise_panels rows are laid along each section. Row m
// spans some [psi_a, psi_b]; its bound vortex sits at that row's own quarter
// chord and its control point at that row's three-quarter chord. This is the
// classic horseshoe-per-panel arrangement (Katz & Plotkin), which the
// existing solver core already handles unmodified: it operates on an
// arbitrary panel list and never assumed one row. Every panel in a chordwise
// stack shares a StripIndex, so the solver reports the stack as ONE spanwise
// station with the summed section circulation.
//
// --- spacing and breakpoints ------------------------------------------------
// Panels are distributed by LatticeOptions (uniform on both axes by
// default -- see Spacing for the measurements behind that choice). The
// distribution is PIECEWISE between breakpoints the mesh must land on:
//   spanwise -- every planform station eta, plus both eta edges of every
//               wing-bound control surface;
//   chordwise -- 0, 1, and any hinge line crossing the strip.
// A section's panel budget is redistributed among its sub-intervals in
// proportion to their widths rather than added to, so panel count does not
// grow with the number of control surfaces.
//
// Landing on those breakpoints is a correctness matter, not just an
// accuracy one: without them a control surface's extent snaps to whatever
// panel boundary is nearest, and its span -- hence its control authority --
// becomes a property of the mesh instead of the contract.
//
// --- control surfaces -------------------------------------------------------
// Where a wing-bound control surface covers a strip, the chordwise rows are
// redistributed so that a panel EDGE falls exactly on the hinge line rather
// than a panel straddling it -- a straddling panel would have to be either
// wholly deflected or wholly not, making the hinge position an artifact of
// the mesh. Rows are then split between the fore and aft portions in
// proportion to their chord fractions.
//
// Panels aft of the hinge are rigidly rotated about the hinge line: both the
// lattice geometry AND the flow-tangency normal, so the deflection is a real
// change in shape rather than a linearized normal-only approximation. Since
// the rotation is rigid it preserves camber-line arc length, so a deflected
// panel's Area is unchanged -- only its projection moves.
//
// Duct-jet-bound surfaces are skipped entirely. ControlSurface.h is
// emphatic that they are not lifting surfaces and must never enter the wing
// lattice; their eta band measures a duct-exit radius fraction, a different
// quantity that merely shares the name.
//
// Hinge axes arrive in the contract's own aetherion_body_frd frame and are
// converted here, at ingest, exactly as ADR-0016 and ControlSurface.h
// prescribe -- rotating the axis vector and keeping the right-hand rule
// about it, never re-deriving a sign convention from prose.
//
// With one chordwise row, a symmetric section and no control surfaces, all
// of this reproduces Solver::BuildWing() exactly, which TestPanelBuilder
// asserts. Two knowing deviations from that builder, both toward being more
// geometrically consistent:
//   - The section frame here is orthonormal: the chordwise axis is twisted
//     about global y, and the surface-normal axis comes from crossing it
//     with the dihedral-tilted spanwise axis. BuildWing() instead rotates
//     the normal by twist and dihedral in sequence while leaving the
//     chordwise offset untwisted by dihedral. The two agree exactly when
//     only one of twist/dihedral is nonzero, and differ by O(sin(twist) *
//     (1 - cos(dihedral))) when both are.
//   - Spanwise panel edges follow the interpolated planform, so a cranked
//     wing gets a genuine kink rather than a fitted straight taper.

#pragma once

#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/Lattice/Panel.h"

#include <cstddef>
#include <utility>
#include <vector>

namespace Aeolion::PanelBuilder {

// --- module constants ------------------------------------------------------
// A chordwise row count below this cannot resolve a hinge line at all --
// there is no interior edge to place it on -- so a covered strip keeps its
// uniform distribution and the surface simply does not articulate.
inline constexpr int MinRowsToResolveHinge = 2;

// Two breakpoints closer than this (in eta) are the same breakpoint. Without
// a merge tolerance, a control surface edge that lands a hair off a planform
// station would spawn a sliver interval with its own panels.
inline constexpr double BreakpointMergeTolerance = 1e-9;

// Eta 0 and 1 are always planform stations, so a control surface edge at
// either is already a breakpoint and must not be inserted twice.
inline constexpr double EtaInteriorMin = 0.0;
inline constexpr double EtaInteriorMax = 1.0;

// Fewest panels any sub-interval may receive when a section's budget is
// split at a breakpoint.
inline constexpr int MinPanelsPerSubInterval = 1;

// --- commanded control deflection ------------------------------------------
// Positive by the right-hand rule about the surface's own hinge axis (for
// the usual spanwise axis that is trailing-edge down).
//
// The two semi-spans are commanded separately because ONE contract entry
// describes the band on BOTH of them -- its eta is a semi-span fraction --
// and the interesting control modes differ in exactly this sign. Equal
// angles give a flap/elevator; opposite angles give an aileron. There is no
// safe default that covers both, so both are stated.
//
// Surfaces are addressed by INDEX into HandoffContract::ControlSurfaces,
// never by name: ControlSurface.h warns that names are not unique (a
// four-vane duct emits four entries all called "vane").
struct ControlDeflection {
    std::size_t SurfaceIndex = 0;
    double RightAngleDeg = 0.0; // +y semi-span
    double LeftAngleDeg = 0.0;  // -y semi-span
};

[[nodiscard]] inline ControlDeflection Symmetric(std::size_t surfaceIndex, double angleDeg) {
    return {surfaceIndex, angleDeg, angleDeg};
}

[[nodiscard]] inline ControlDeflection Antisymmetric(std::size_t surfaceIndex, double angleDeg) {
    return {surfaceIndex, angleDeg, -angleDeg};
}

// --- spacing ----------------------------------------------------------------
// How panels are distributed within one interval between breakpoints. The
// handoff schema has no spacing field, so this is the consumer's choice,
// not the producer's.
//
// Cosine clusters panels toward both ends of an interval, which is where
// the loading gradient is steep. It is the textbook choice and it is NOT
// the default here, for a measured reason.
//
// MEASURED WARNING. Cosine clustering shrinks the end panels like 1/N^2
// while the mean panel is 1/N, so the panel size ratio grows without bound
// as the mesh refines -- 57:1 at 90 spanwise panels. This lattice does not
// tolerate that. Against a converged uniform reference (rectangular AR=8,
// cambered, alpha=3deg, 4320 panels), where the trusted answer is
// CL = 0.361482 at an LU pivot ratio of 0.90:
//
//     spanwise cosine only    CL = 0.363477   pivot 2.1e-02   (+0.6%)
//     chordwise cosine only   CL = 0.375924   pivot 1.9e-01   (+4.0%)
//     both                    CL = 0.383217   pivot 3.4e-03   (+6.0%)
//
// The error is not a different discretization converging elsewhere: it
// grows with refinement, tracks the collapse in pivot ratio, and at higher
// clustering the solve fails outright (a variant clustering only outboard
// returned CL = -4.3). Chordwise is the worse offender despite the milder
// conditioning hit, because the 1/4-3/4 collocation rule is what makes a
// uniform chordwise division reproduce thin-airfoil theory exactly; on a
// non-uniform division that property is lost.
//
// Cosine is kept because it is a modest win at low panel counts, but it is
// opt-in and should not be used on fine meshes without checking
// LUFactorization::MinPivotRatio.
enum class Spacing { Uniform, Cosine };

struct LatticeOptions {
    Spacing Spanwise = Spacing::Uniform;
    Spacing Chordwise = Spacing::Uniform;
};

// --- the builder ------------------------------------------------------------
class LatticeBuilder {
public:
    // The contract is copied: a builder outlives the expression that made
    // it in every intended use (LatticeBuilder(LoadHandoff(path))), and the
    // copy is negligible beside the O(N^3) solve it feeds.
    explicit LatticeBuilder(Geometry::HandoffContract contract, LatticeOptions options = {});

    // Command a control surface. Repeat calls for the same surface replace
    // the previous command, so a sweep can just re-Deflect and re-Build.
    LatticeBuilder& Deflect(const ControlDeflection& deflection);
    LatticeBuilder& ClearDeflections();

    [[nodiscard]] std::vector<Solver::Panel> Build() const;

    // The spanwise panel boundaries this builder will use, in eta. Exposed
    // because "did the mesh land on the control surface edge?" is a
    // question worth being able to ask without rebuilding the lattice.
    [[nodiscard]] const std::vector<double>& BoundaryEtas() const { return m_BoundaryEtas; }
    [[nodiscard]] const Geometry::HandoffContract& Contract() const { return m_Contract; }

private:
    // One accumulated station on the (implicit right) semi-span: the
    // chord-line quarter-chord point plus everything needed to place lattice
    // points on the section there. All of these are functions of |y| only,
    // so a single semi-span pass is mirrored for both sides.
    struct SpanStation {
        double Eta = 0.0;
        double y = 0.0;
        double QuarterChordX = 0.0;
        double QuarterChordZ = 0.0;
        double Chord = 0.0;
        double TwistDeg = 0.0;
        double DihedralDeg = 0.0;
        double SweepQuarterChordDeg = 0.0; // only needed to seed the next accumulation step
    };

    // The section's local orthonormal frame at a station. ChordDir runs aft
    // along the (twisted) chord line; UpDir is perpendicular to it in the
    // section plane, tilted by dihedral -- the direction camber is measured
    // along, and the panel normal when the camber slope is zero.
    struct SectionFrame {
        Solver::Vec3 ChordDir;
        Solver::Vec3 UpDir;
    };

    // The hinge a given strip articulates about, already resolved to angle.
    // Inactive when no wing-bound control surface covers the strip.
    struct HingeSpec {
        bool Active = false;
        double PsiHinge = 1.0;
        Solver::Vec3 Axis{0, 1, 0};
        double RightAngleRad = 0.0;
        double LeftAngleRad = 0.0;
    };

    // -- pure geometry helpers, independent of any particular contract --
    [[nodiscard]] static double SpacingFraction(Spacing spacing, int index, int count);
    [[nodiscard]] static std::vector<int> AllocateProportional(const std::vector<double>& widths, int total,
                                                               int minimum);
    [[nodiscard]] static Geometry::PlanformStation InterpolatePlanform(
        const std::vector<Geometry::PlanformStation>& stations, double eta);
    [[nodiscard]] static Solver::Vec3 HingeAxisToSolverFrame(const Math::Vec3& contractAxis);
    [[nodiscard]] static Solver::Vec3 MirrorHingeAxis(const Solver::Vec3& axis);
    [[nodiscard]] static SectionFrame FrameAt(const SpanStation& station, double sign);
    [[nodiscard]] static SpanStation MidStation(const SpanStation& a, const SpanStation& b);

    // -- construction-time work, done once --
    [[nodiscard]] std::vector<double> ComputeBreakpoints() const;
    [[nodiscard]] std::vector<double> ComputeBoundaryEtas() const;
    [[nodiscard]] std::vector<SpanStation> ComputeSemiSpan() const;

    // -- per-Build work --
    [[nodiscard]] HingeSpec HingeForStrip(double eta) const;
    [[nodiscard]] std::vector<std::pair<double, double>> ChordwiseRowBounds(const HingeSpec& hinge) const;
    [[nodiscard]] Solver::Vec3 SurfacePoint(const SpanStation& station, double sign, double psi) const;
    void EmitStrip(const SpanStation& inner, const SpanStation& outer, double sign, int stripIndex,
                   std::vector<Solver::Panel>& panels) const;

    Geometry::HandoffContract m_Contract;
    LatticeOptions m_Options;
    std::vector<ControlDeflection> m_Deflections;
    std::vector<double> m_BoundaryEtas;
    std::vector<SpanStation> m_SemiSpan;
};

} // namespace Aeolion::PanelBuilder
