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
#include "Aeolion/Lattice/SourcePanel.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <cstddef>
#include <functional>
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

// Fewest azimuthal divisions that still describe a closed surface of
// revolution rather than a sliver.
inline constexpr int MinBodySectors = 3;

// Panels smaller than these are the degenerate ones a collapsing ring
// produces at a pointed nose or closed tail; they carry no surface and
// would only put near-zero rows into the influence matrix.
inline constexpr double BodyDegenerateLength = 1e-12;
inline constexpr double BodyDegenerateArea = 1e-18;

// Radial subdivisions of the base cap. Potential flow is singular at a
// sharp base rim, and a single fan of full-radius wedges resolves it badly
// enough that a closed body visibly stops satisfying d'Alembert.
inline constexpr int BodyBaseRadialRings = 4;

// A cut beyond this fraction of semi-span would leave almost no wing; a
// body that wide means the placement or the radius law is wrong, and
// returning an empty lattice would be a confusing way to say so.
inline constexpr double MaxTrimEta = 0.9;

inline constexpr double Two3rds = 2.0 / 3.0; // radial centroid of a triangle spanning the axis to the rim

// Surface tags. The base is tagged apart from the rest of the body because
// it is the face that will carry the duct efflux rather than a solid wall.
inline constexpr const char* BodySurfaceName = "fuselage";
inline constexpr const char* BodyBaseSurfaceName = "fuselage_base";

// Fewest azimuthal divisions that still describe a closed duct ring rather
// than a sliver -- same reasoning as MinBodySectors.
inline constexpr int MinDuctSectors = 3;

// Surface tags for the duct's four faces: the outer wall, the inner bore
// wall the slipstream passes through, and the two annular end caps.
inline constexpr const char* DuctOuterSurfaceName = "duct_outer";
inline constexpr const char* DuctInnerSurfaceName = "duct_inner";
inline constexpr const char* DuctLeadingCapSurfaceName = "duct_leading_cap";
inline constexpr const char* DuctTrailingCapSurfaceName = "duct_trailing_cap";

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
/** A commanded control-surface deflection, addressed by index into HandoffContract::ControlSurfaces. */
struct ControlDeflection {
    std::size_t SurfaceIndex = 0;
    double RightAngleDeg = 0.0; ///< +y semi-span.
    double LeftAngleDeg = 0.0;  ///< -y semi-span.
};

/** Equal deflection on both semi-spans (flap/elevator mode). */
[[nodiscard]] inline ControlDeflection Symmetric(std::size_t surfaceIndex, double angleDeg) {
    return {surfaceIndex, angleDeg, angleDeg};
}

/** Opposite deflection on each semi-span (aileron mode). */
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
/** How panels are distributed within one interval between breakpoints. */
enum class Spacing { Uniform, Cosine };

/** Discretization and body-handling options for LatticeBuilder. */
struct LatticeOptions {
    Spacing Spanwise = Spacing::Uniform;
    Spacing Chordwise = Spacing::Uniform;

    /**
     * Whether to panel the fuselage at all. A wing-only study is a
     * legitimate thing to want, and the body roughly doubles the unknown
     * count -- which the dense solve pays for cubically.
     */
    bool IncludeBody = true;

    /**
     * Cut the wing lattice off where it enters the fuselage. A lifting
     * panel whose control point is inside a solid body is not a meaningful
     * boundary condition, so the default is to trim.
     *
     * The cut is at the body radius evaluated at the ROOT LEADING EDGE
     * station -- a single straight cut rather than one following the body's
     * curvature down the chord. That is the conventional wing-body
     * treatment and it is what the placement anchor names, so the two agree
     * by construction. Requires a stated placement; without one there is no
     * defined relationship to cut against and nothing is trimmed.
     */
    bool TrimWingAtBody = true;

    /**
     * Carry the wing's bound circulation THROUGH the fuselage.
     *
     * Trimming alone loses lift that the real aircraft keeps: a source
     * distribution cannot carry circulation, and lift is circulation, so a
     * sources-only body recovers only about three quarters of the trimmed
     * root load. The classical fix (Pitts, Nielsen & Kaattari) is a
     * carry-through vortex spanning the body at the strength of the
     * innermost exposed strip.
     *
     * Implemented by EXTENDING that strip's bound segment inboard to the
     * centreline rather than adding an unknown: the strength is then equal
     * to its neighbour's by construction instead of by a constraint row, so
     * the system stays one block of tangency equations. Tangency is still
     * enforced only on the exposed wing, since the extended part of the
     * segment runs inside the body where a boundary condition would be
     * meaningless.
     */
    bool CarryThroughLift = true;

    /**
     * Circumferential divisions around the body. The handoff schema states
     * an axial station list but no azimuthal resolution, so this is the
     * consumer's choice. The body is axisymmetric, but the FLOW around it
     * is not once there is incidence or sideslip, which is the whole reason
     * to resolve it circumferentially rather than treat it as rings.
     */
    int BodyCircumferentialPanels = 16;

    /**
     * Whether to panel the duct at all. Cheap relative to the fuselage, but
     * still not free, and a contract with no `duct` block builds an empty
     * list regardless of this flag.
     */
    bool IncludeDuct = true;

    /** Circumferential divisions around the duct, same role as BodyCircumferentialPanels. */
    int DuctCircumferentialPanels = 16;

    /**
     * Axial divisions along the duct's chord, on the two cylindrical walls
     * only (the end caps are a single radial band -- see BuildDuct()). The
     * duct schema states a single chord, not a station list like the
     * fuselage, so how finely to resolve that chord axially is this
     * consumer's choice rather than the producer's.
     */
    int DuctAxialPanels = 4;
};

// --- base efflux ------------------------------------------------------------
// The base that BuildBody() caps is a duct exit, not a wall: fluid leaves
// through it. This turns those panels from a solid cap into an outflow
// boundary carrying the propeller jet.
//
// WHAT IS PRESCRIBED. Lattice::SourcePanel::PrescribedNormalVelocity is the
// TOTAL normal velocity the boundary condition enforces, so what goes there
// is the jet's exit speed along the outward normal -- the oncoming flow plus
// the slipstream's own induced velocity. The solver then asks the sources
// for the difference between that and what the freestream already supplies,
// which is exactly the excess the jet has to add. A solid wall prescribes
// zero and gets the ordinary tangency condition, so the two cases are the
// same equation with a different right-hand side, not different treatments.
//
// WHY A REFERENCE CONDITION. The exit velocity of a ducted fan is set by the
// propeller's operating point, not by the airframe's attitude, so it is
// baked in once rather than recomputed per alpha. That is an approximation
// and worth stating plainly: the freestream's contribution to the exit
// velocity does vary with incidence (as cos(alpha) at the base), and this
// ignores that variation. It is small next to the jet's own induced velocity
// for any meaningful thrust setting, and it keeps the influence matrix
// flight-condition-independent, which is what makes a derivative sweep cheap.
//
// The slipstream is taken as a plain callable so this does not drag any
// propeller-model dependency into the panel builder; e.g. the external
// BEMT project's SlipstreamField satisfies it after a Vec3 field-wise
// conversion. With no propeller state, simply never call this: the base stays
// the solid cap BuildBody() produced, which is a closed-body answer -- wrong
// in a known direction rather than an unknown one.
//
/**
 * Prescribe the duct exit's outflow (freestream + slipstream induced
 * velocity) on the base panels of body, in place. Returns the number of
 * base panels affected, so a caller can tell the difference between
 * "applied" and "there was no base".
 */
[[nodiscard]] int ApplyBaseEfflux(std::vector<Lattice::SourcePanel>& body,
                                  const std::function<Math::Vec3(const Math::Vec3&)>& slipstream,
                                  const Math::Vec3& referenceFreestream);

// --- propeller blade lattice ------------------------------------------------
// The panel-method propeller: every blade meshed as horseshoe-vortex panels
// (one Weissinger row per radial strip, exactly the Panel vocabulary the
// wing uses) and solved by the SAME Solver::Solve as the airframe -- the
// rotation enters through FreestreamConditions::p, the solver's roll rate
// about +x, whose kinematic-velocity term hands every blade panel its true
// Omega x r tangential onset flow. Quasi-steady and potential-flow: the
// torque a solve reports is induced torque, with no profile-drag
// contribution.
//
// Conventions: the prop spins at +Omega about +x (right-hand rule), the
// freestream arrives along +x (solver body axes, x aft), so thrust comes
// out along -x: Thrust = -SolveResult::Di, shaft torque = -SolveResult::Mx
// (with RefPoint at the hub), power = torque * Omega. Blades sit on the
// CST camber surface when prop.Sections carries one (chords wrapped on
// their radius cylinders -- see the .cpp), flat chord lines otherwise.
//
// The wake direction is why this function needs the operating point, not
// just geometry: each trailing leg leaves along the LOCAL kinematic
// velocity at its root (axial inflow + Omega x r tangential sweep, with
// the momentum-theory hover inflow as the axial floor), the linearized
// helix. Trailing purely axially instead is not a small error at hover --
// the local wind there is almost entirely tangential, the legs would
// leave near-perpendicular to the flow, and the lattice's
// quarter/three-quarter-chord geometry stops meaning anything (in
// practice the circulation solve saw-tooths radially and diverges with
// refinement). Rebuild the lattice whenever rpm or inflow change; it is
// cheap.
//
// ONE chordwise row per strip, deliberately (Weissinger): the chords are
// wrapped on their radius cylinders, and a straight trailing leg leaves a
// curved chord immediately -- the surface falls away from its tangent
// line by s^2/2r, so a forward row's legs would graze the aft rows'
// control points at fractions of a millimetre and the chordwise
// circulation saw-tooths (found empirically). Camber still enters exactly
// where thin-airfoil theory wants it, through the mean-line slope at the
// 3/4-chord boundary condition; resolving the chordwise LOADING
// distribution needs a vortex-ring wake or finite-core filaments, which
// belongs with the future viscous/boundary-layer work.
/**
 * Mesh a metric propeller into horseshoe-vortex panels: one Weissinger
 * row per radial strip (station pair), per blade, on the CST camber
 * surface when prop.Sections carries one, with trailing legs along the
 * local relative wind of the stated operating point (axialSpeed [m/s]
 * inflow along +x, omega [rad/s] about +x). Solve with
 * FreestreamConditions{ Vinf = axialSpeed, p = omega, RefPoint = hub
 * center (the origin) }.
 */
[[nodiscard]] std::vector<Lattice::Panel> BuildPropellerLattice(const Geometry::Propeller& prop,
                                                                double axialSpeed, double omega);

/**
 * The per-strip section frames and section data for the Level-2 viscous
 * coupling (Solver::SolveViscousCoupled), aligned one-to-one with
 * BuildPropellerLattice's panels: chord/lift directions at mid-chord of
 * each radial strip, the local chord and width, and the thin-airfoil
 * zero-lift angle of the strip's CST camber line -- the section model
 * must carry the camber the coupling supersedes in the lattice's own
 * boundary condition. Pure geometry: independent of the operating point.
 */
[[nodiscard]] std::vector<Solver::StripSection> BuildPropellerStrips(const Geometry::Propeller& prop);

// --- the builder ------------------------------------------------------------
/**
 * Builds the solver's lattice from a parsed geometry handoff, on the
 * camber surface the handoff's CST sections describe, with wing-bound
 * control surfaces deflected about their stated hinge lines. See the file
 * header for the full method.
 */
class LatticeBuilder {
public:
    /**
     * The contract is copied: a builder outlives the expression that made
     * it in every intended use (LatticeBuilder(LoadHandoff(path))), and the
     * copy is negligible beside the O(N^3) solve it feeds.
     */
    explicit LatticeBuilder(Geometry::HandoffContract contract, LatticeOptions options = {});

    /**
     * Command a control surface. Repeat calls for the same surface replace
     * the previous command, so a sweep can just re-Deflect and re-Build.
     */
    LatticeBuilder& Deflect(const ControlDeflection& deflection);
    /** Clear all commanded deflections. */
    LatticeBuilder& ClearDeflections();

    /** Build the wing lattice with the currently commanded deflections applied. */
    [[nodiscard]] std::vector<Solver::Panel> Build() const;

    /**
     * The fuselage as source panels. Empty when the contract carries no
     * body or LatticeOptions::IncludeBody is false, which is exactly the
     * degenerate case the solver's blocked system reduces to.
     */
    [[nodiscard]] std::vector<Lattice::SourcePanel> BuildBody() const;

    /**
     * The duct as source panels: two concentric cylinders (an outer wall and
     * an inner bore wall, straight-walled over the stated chord -- see
     * Geometry::DuctGeometry.h for what the schema does and does not carry),
     * capped front and back by annuli. Disjoint from BuildBody()'s fuselage,
     * with its own surface tags. Empty when the contract carries no duct or
     * LatticeOptions::IncludeDuct is false.
     */
    [[nodiscard]] std::vector<Lattice::SourcePanel> BuildDuct() const;

    /**
     * Semi-span fraction the wing is trimmed at, or 0 when nothing is
     * trimmed. Exposed because "how much of the wing did the body eat?" is
     * a question worth asking without diffing panel counts.
     */
    [[nodiscard]] double TrimEta() const { return m_TrimEta; }

    /**
     * Area of the FULL trapezoidal planform, including the part buried in
     * the fuselage, computed from the contract's stations rather than by
     * summing panels.
     *
     * This is the reference area coefficients must be normalized by, and it
     * must NOT follow the trim: CL is conventionally referred to gross
     * planform area, so letting a trimmed panel sum drive it would rescale
     * CL by whatever fraction was cut -- a silent redefinition of exactly
     * the kind Solver::Panel's Area/PlanformArea split exists to prevent.
     */
    [[nodiscard]] double GrossPlanformArea() const;

    /**
     * The spanwise panel boundaries this builder will use, in eta. Exposed
     * because "did the mesh land on the control surface edge?" is a
     * question worth being able to ask without rebuilding the lattice.
     */
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
    [[nodiscard]] double ComputeTrimEta() const;
    [[nodiscard]] Solver::Vec3 ComputePlacementOffset() const;
    [[nodiscard]] std::vector<double> ComputeBreakpoints() const;
    [[nodiscard]] std::vector<double> ComputeBoundaryEtas() const;
    [[nodiscard]] std::vector<SpanStation> ComputeSemiSpan() const;

    // -- per-Build work --
    [[nodiscard]] HingeSpec HingeForStrip(double eta) const;
    [[nodiscard]] std::vector<std::pair<double, double>> ChordwiseRowBounds(const HingeSpec& hinge) const;
    [[nodiscard]] Solver::Vec3 SurfacePoint(const SpanStation& station, double sign, double psi) const;
    void EmitStrip(const SpanStation& inner, const SpanStation& outer, double sign, int stripIndex,
                   bool carryThrough, std::vector<Solver::Panel>& panels) const;

    Geometry::HandoffContract m_Contract;
    LatticeOptions m_Options;
    Solver::Vec3 m_PlacementOffset{0, 0, 0}; // solver frame, from the contract's root-LE anchor
    double m_TrimEta = 0.0;                  // semi-span fraction cut away at the body surface
    std::vector<ControlDeflection> m_Deflections;
    std::vector<double> m_BoundaryEtas;
    std::vector<SpanStation> m_SemiSpan;
};

} // namespace Aeolion::PanelBuilder
