// Solver/AttachmentLine.h
//
// The wing's attachment line: the spanwise locus of stagnation points where
// the flow divides between the upper and lower surfaces, and the boundary
// layer of each begins.
//
// --- the division of labour ---------------------------------------------------
// The lattice cannot answer this on its own -- it is a zero-thickness camber
// sheet and the flow never stops on it (see SectionPanelMethod.h). What it
// owns is the INDUCED field, and therefore what each strip sees. What the
// section owns is the thickness, and therefore where the flow actually
// stops. This file is the seam: it reads a local velocity per strip out of
// the 3-D solve, poses a 2-D section problem in the right plane, and maps
// the answer back into three dimensions.
//
// --- why the leading-edge-NORMAL plane ----------------------------------------
// The section cuts a handoff carries are streamwise, taken at a span
// fraction. That is the wrong plane to pose an attachment-line problem in,
// and the error is not a small one once there is sweep or sideslip.
//
// Infinite-swept-wing theory says the flow near a swept leading edge
// separates into two independent parts: a two-dimensional problem in the
// plane NORMAL to the leading edge, and a spanwise velocity along it that
// is simply convected. The normal problem places the attachment line; the
// spanwise velocity does not affect where it is, but it is the whole of
// what decides whether the attachment line stays laminar. Posing the
// section streamwise mixes the two.
//
// Sideslip is what makes this load-bearing rather than academic. The
// EFFECTIVE sweep -- the angle between the local flow and the plane normal
// to the leading edge -- is not the geometric sweep once beta is nonzero,
// and it changes in OPPOSITE directions on the two wings: at positive
// sideslip one wing's leading edge is swept further from the flow and the
// other's is raked toward it. So the attachment line moves asymmetrically,
// its state differs left to right, and a streamwise formulation cannot show
// that at all -- it produces a symmetric answer to an asymmetric question.
// Everything here therefore keys off the leading edge's own direction and
// the local velocity resolved against it.
//
// --- what comes out -------------------------------------------------------------
// Per strip: the attachment point in 3-D, the edge-velocity distribution
// along both surfaces measured from it, the strain rate that starts a
// boundary-layer march, and the attachment-line Reynolds number that decides
// whether that march may start laminar at all.

#pragma once

#include "Aeolion/Geometry/AirfoilSection.h"
#include "Aeolion/Geometry/SectionContour.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Solver/SectionPanelMethod.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace Aeolion::Solver {

// --- attachment-line transition (Poll) -----------------------------------------
// The swept attachment line has its own boundary layer, driven by the
// spanwise velocity and thinned by the normal-plane strain rate. Its state
// is governed by
//
//     Rbar = W_e * eta / nu,     eta = sqrt(nu / (dU_e/ds)),
//
// with W_e the velocity ALONG the leading edge and eta the attachment-line
// length scale. Two thresholds matter and they mean different things:
// below the first, a disturbance introduced at the root (a fuselage
// junction, a trip, a de-icing boot) decays as it travels outboard; above
// it, the disturbance propagates and CONTAMINATES the whole attachment
// line, turning the entire wing's leading edge turbulent regardless of what
// the chordwise flow would have done. The second threshold is where the
// attachment line goes turbulent on its own without any such disturbance.
inline constexpr double AttachmentLineContaminationReynolds = 245.0;
inline constexpr double AttachmentLineTransitionReynolds = 583.0;

/**
 * Momentum thickness of the swept attachment-line boundary layer as a
 * multiple of eta -- the constant of the swept Hiemenz similarity solution.
 * Distinct from the CHORDWISE layer's starting thickness, which is Thwaites'
 * stagnation value (see StagnationMomentumThickness): the two layers are
 * different flows sharing one strain rate.
 */
inline constexpr double AttachmentLineThicknessConstant = 0.404;

/** What state the attachment line is in at a station. */
enum class AttachmentLineState {
    Laminar,      ///< Below contamination: a root disturbance decays going outboard.
    Contaminated, ///< A root disturbance would propagate along the whole leading edge.
    Turbulent,    ///< Turbulent unaided, with no disturbance needed.
};

/** The attachment-line solution at one spanwise strip. */
struct AttachmentStation {
    int Strip = -1;
    double Eta = 0.0;

    // --- the strip's frame ---------------------------------------------------
    Vec3 LeadingEdge{0, 0, 0};    ///< The leading-edge point at this strip, in solver axes.
    Vec3 LeadingEdgeDir{0, 0, 0}; ///< Unit, along the leading edge toward increasing strip index.
    Vec3 ChordDir{0, 0, 0};       ///< Unit, leading edge -> trailing edge.
    Vec3 Normal{0, 0, 0};         ///< Unit, toward the suction side.

    /**
     * This station sits against a break in the leading edge -- a root kink,
     * a crank, a strake junction -- where the leading edge's direction is
     * discontinuous and the one-sided value away from the break was used
     * (see Detail::LeadingEdgeDirection).
     *
     * Worth propagating rather than hiding: infinite-swept-wing theory
     * assumes a locally straight leading edge and has nothing to say at a
     * break, and a break at the root is exactly where the disturbance that
     * contaminates an attachment line comes from.
     */
    bool AtKink = false;

    // --- the leading-edge-normal problem ---------------------------------------
    /** Geometric sweep of the leading edge relative to the section cut [rad]. */
    double SweepRad = 0.0;
    /**
     * Sweep relative to the LOCAL FLOW, asin(|V . e_LE| / |V|) [rad]. This is
     * the one that differs between the two wings at sideslip, and the one
     * the attachment-line Reynolds number is built from.
     */
    double EffectiveSweepRad = 0.0;
    double NormalChord = 0.0;      ///< Chord in the leading-edge-normal plane [m].
    double AlphaNormalDeg = 0.0;   ///< Incidence in that plane, from the chord line.
    double NormalSpeed = 0.0;      ///< |V| resolved into that plane [m/s].
    double SpanwiseSpeed = 0.0;    ///< W_e, along the leading edge [m/s].

    // --- the answer ------------------------------------------------------------
    bool Found = false;
    Vec3 AttachmentPoint{0, 0, 0}; ///< Where the flow divides, in solver axes.

    /**
     * Surface distance of the attachment point from the leading edge, as a
     * fraction of the NORMAL chord. Positive means it lies on the lower
     * (pressure) surface, which is the usual case at positive incidence.
     */
    double StagnationOffset = 0.0;

    double StrainRate = 0.0;         ///< dU_e/ds at the attachment point [1/s].
    double LengthScale = 0.0;        ///< eta = sqrt(nu / StrainRate) [m].
    double MomentumThickness = 0.0;  ///< Chordwise layer's starting theta_0 [m].
    double AttachmentLineTheta = 0.0; ///< The swept attachment line's own theta [m].
    double AttachmentLineReynolds = 0.0; ///< Rbar.
    AttachmentLineState State = AttachmentLineState::Laminar;

    /**
     * The whole 2-D solve: U_e(s) on both surfaces measured from the
     * attachment point, ready for an integral boundary-layer march. Its arc
     * lengths are fractions of NormalChord and its speeds multiples of
     * NormalSpeed.
     */
    SectionSolution Section;
};

/** The attachment line across a lifting surface, one entry per strip. */
struct AttachmentLine {
    std::vector<AttachmentStation> Stations;
};

/** Tuning for an attachment-line computation. */
struct AttachmentLineOptions {
    double KinematicViscosity = SeaLevelKinematicViscosity; ///< [m^2/s]
    int PointsPerSurface = Geometry::DefaultContourPointsPerSurface;
};

/**
 * Below this agreement between the leading-edge segments either side of a
 * station, the leading edge is taken to BREAK there rather than curve.
 * Five degrees: sweep and taper bend a leading edge by far less than that
 * between neighbouring strips, and a planform break bends it by far more.
 */
inline constexpr double LeadingEdgeKinkCosine = 0.9962; // cos(5 degrees)

namespace Detail {

/**
 * The leading-edge point of a strip. A Weissinger row carries its bound
 * vortex on the quarter-chord line, so the leading edge is a quarter chord
 * ahead of the bound segment's midpoint.
 */
[[nodiscard]] inline Vec3 StripLeadingEdge(const Panel& panel, const StripSection& strip) {
    const Vec3 mid = (panel.A + panel.B) * Math::Half;
    return mid - strip.ChordDir * (Math::QuarterChord * strip.Chord);
}

/**
 * The direction of the leading edge at station i, differenced from
 * neighbouring stations -- and refusing to average across a planform break.
 *
 * A swept wing's leading edge is not differentiable at the root: it runs aft
 * going outboard on BOTH sides, so its direction flips discontinuously
 * across the centreline. A plain central difference there returns the mean
 * of the two, which is very nearly unswept, and hands the two innermost
 * strips an effective sweep of about half the true one. That is not a small
 * error in an unimportant place: the root is precisely where a disturbance
 * enters the attachment line and where contamination is seeded.
 *
 * So when the two one-sided differences disagree beyond
 * LeadingEdgeKinkCosine, the station is at a break, and the side to trust is
 * the one that agrees with ITS own outer neighbour -- the side that does not
 * cross the kink. `atKink` reports the situation either way, because a
 * caller reading attachment-line contamination wants to know a station sits
 * against a planform break.
 */
[[nodiscard]] inline Vec3 LeadingEdgeDirection(const std::vector<Vec3>& points, std::size_t i,
                                               bool& atKink) {
    atKink = false;
    const std::size_t n = points.size();
    const auto segment = [&](std::size_t a, std::size_t b) { return (points[b] - points[a]).Normalized(); };

    if (i == 0) return segment(0, 1);
    if (i + 1 >= n) return segment(n - 2, n - 1);

    const Vec3 inboard = segment(i - 1, i);
    const Vec3 outboard = segment(i, i + 1);
    if (Dot(inboard, outboard) >= LeadingEdgeKinkCosine) return (inboard + outboard).Normalized();

    atKink = true;
    double inboardScore = -2.0, outboardScore = -2.0;
    if (i >= 2) inboardScore = Dot(inboard, segment(i - 2, i - 1));
    if (i + 2 < n) outboardScore = Dot(outboard, segment(i + 1, i + 2));
    return (inboardScore >= outboardScore) ? inboard : outboard;
}

} // namespace Detail

/**
 * Locate the attachment line from per-strip local velocities.
 *
 * `panels` and `strips` must align one-to-one, the same contract
 * SolveViscousCoupled states: one Weissinger row per spanwise strip. So must
 * `localVelocity`, which is the velocity each strip sees -- kinematic plus
 * induced, its own bound vortex excluded. That is exactly
 * FlowField::BoundMidpointVelocity, and it is taken as an argument rather
 * than recomputed so this stays usable with a viscous-coupled solve, whose
 * converged strip velocities are its own.
 *
 * Strips must be ordered along the span, because the leading edge's own
 * direction is differenced from neighbouring stations -- there is no other
 * way to know which way a leading edge runs.
 */
[[nodiscard]] inline AttachmentLine ComputeAttachmentLine(
    const std::vector<Panel>& panels, const std::vector<StripSection>& strips,
    const std::vector<Vec3>& localVelocity, const std::vector<Geometry::AirfoilSection>& sections,
    const AttachmentLineOptions& options = {}) {
    AttachmentLine line;
    const std::size_t n = panels.size();
    if (n < 2 || strips.size() != n || localVelocity.size() != n) return line;

    // Leading-edge points first: the leading edge's DIRECTION is a difference
    // between neighbouring stations, so every point has to exist before any
    // direction does.
    std::vector<Vec3> leadingEdge(n);
    for (std::size_t i = 0; i < n; ++i) leadingEdge[i] = Detail::StripLeadingEdge(panels[i], strips[i]);

    line.Stations.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const StripSection& strip = strips[i];
        AttachmentStation station;
        station.Strip = static_cast<int>(i);
        station.Eta = strip.Eta;
        station.LeadingEdge = leadingEdge[i];
        station.ChordDir = strip.ChordDir.Normalized();
        station.Normal = strip.LiftDir.Normalized();

        // Leading-edge direction: central where the leading edge is smooth,
        // one-sided at the tips and across a planform break.
        bool atKink = false;
        const Vec3 along = Detail::LeadingEdgeDirection(leadingEdge, i, atKink);
        station.AtKink = atKink;
        if (!(along.Norm() > Math::Half)) {
            line.Stations.push_back(station);
            continue;
        }
        station.LeadingEdgeDir = along;

        // --- the leading-edge-normal frame -------------------------------------
        const Vec3 eLE = station.LeadingEdgeDir;
        // Geometric sweep: how far the leading edge is raked out of the plane
        // perpendicular to the chord line.
        const double sinSweep = std::clamp(Dot(eLE, station.ChordDir), Math::UnitClampLo, Math::UnitClampHi);
        station.SweepRad = std::asin(sinSweep);
        const double cosSweep = std::sqrt(std::max(1.0 - sinSweep * sinSweep, 0.0));
        if (!(cosSweep > Math::Tiny)) { // a leading edge parallel to the chord has no normal plane
            line.Stations.push_back(station);
            continue;
        }
        station.NormalChord = strip.Chord * cosSweep;

        // Chord direction and normal, resolved into the leading-edge-normal
        // plane. The panel normal is already perpendicular to the surface, so
        // removing its (tiny) spanwise part only cleans up sweep-induced skew.
        Vec3 chordNormalPlane = station.ChordDir - eLE * Dot(station.ChordDir, eLE);
        Vec3 normalInPlane = station.Normal - eLE * Dot(station.Normal, eLE);
        if (!(chordNormalPlane.Norm() > Math::Tiny) || !(normalInPlane.Norm() > Math::Tiny)) {
            line.Stations.push_back(station);
            continue;
        }
        chordNormalPlane = chordNormalPlane.Normalized();
        normalInPlane = normalInPlane.Normalized();

        // --- resolve the local velocity ----------------------------------------
        const Vec3 velocity = localVelocity[i];
        const double speed = velocity.Norm();
        const double spanwise = Dot(velocity, eLE);
        const Vec3 inPlane = velocity - eLE * spanwise;

        station.SpanwiseSpeed = std::fabs(spanwise);
        station.NormalSpeed = inPlane.Norm();
        station.EffectiveSweepRad =
            (speed > Math::Tiny)
                ? std::asin(std::clamp(station.SpanwiseSpeed / speed, 0.0, Math::UnitClampHi))
                : 0.0;
        if (!(station.NormalSpeed > Math::Tiny)) {
            line.Stations.push_back(station);
            continue;
        }

        // Incidence in the normal plane. Positive when the oncoming flow has
        // a component along the surface normal, which is the sense the
        // solver's own freestream carries at positive alpha.
        const double alphaNormal =
            std::atan2(Dot(inPlane, normalInPlane), Dot(inPlane, chordNormalPlane));
        station.AlphaNormalDeg = Math::RadToDeg(alphaNormal);

        // --- the section problem -------------------------------------------------
        const Geometry::SectionContour streamwise =
            Geometry::ContourAt(sections, strip.Eta, options.PointsPerSurface);
        if (!streamwise.Valid()) { // no section data means no thickness, hence no stagnation point
            line.Stations.push_back(station);
            continue;
        }
        const Geometry::SectionContour normalContour =
            Geometry::ToLeadingEdgeNormal(streamwise, station.SweepRad);

        station.Section = SolveSectionContour(normalContour, alphaNormal);
        if (!station.Section.Valid || !station.Section.StagnationFound) {
            line.Stations.push_back(station);
            continue;
        }

        // --- map the answer back into three dimensions -------------------------
        const SectionSolution& section = station.Section;
        station.Found = true;
        station.AttachmentPoint = station.LeadingEdge +
                                  chordNormalPlane * (section.StagnationPsi * station.NormalChord) +
                                  normalInPlane * (section.StagnationZeta * station.NormalChord);

        // Surface offset from the leading edge, signed positive onto the
        // lower surface. The contour's midpoint arc IS the leading edge.
        const double leadingEdgeArc = Math::Half * section.Perimeter;
        station.StagnationOffset = leadingEdgeArc - section.StagnationArc;

        // --- boundary-layer scales ------------------------------------------------
        // The section was nondimensionalized by the NORMAL chord and the
        // NORMAL speed, so the strain rate dimensionalizes with those.
        station.StrainRate = section.StagnationStrain * station.NormalSpeed / station.NormalChord;
        station.MomentumThickness = StagnationMomentumThickness(
            section.StagnationStrain, station.NormalChord, station.NormalSpeed,
            options.KinematicViscosity);

        if (station.StrainRate > Math::Tiny) {
            station.LengthScale = std::sqrt(options.KinematicViscosity / station.StrainRate);
            station.AttachmentLineTheta = AttachmentLineThicknessConstant * station.LengthScale;
            station.AttachmentLineReynolds =
                station.SpanwiseSpeed * station.LengthScale / options.KinematicViscosity;
        }

        station.State = (station.AttachmentLineReynolds >= AttachmentLineTransitionReynolds)
                            ? AttachmentLineState::Turbulent
                        : (station.AttachmentLineReynolds >= AttachmentLineContaminationReynolds)
                            ? AttachmentLineState::Contaminated
                            : AttachmentLineState::Laminar;

        line.Stations.push_back(station);
    }
    return line;
}

/**
 * Convenience overload reading the strip velocities out of a solved field.
 * `panels` must be the same lifting panels the field was solved on, in the
 * same order, since the velocity is taken at each one's own bound-vortex
 * midpoint with its own bound segment excluded.
 */
[[nodiscard]] inline AttachmentLine ComputeAttachmentLine(
    const FlowField& field, const std::vector<Panel>& panels, const std::vector<StripSection>& strips,
    const std::vector<Geometry::AirfoilSection>& sections,
    const AttachmentLineOptions& options = {}) {
    std::vector<Vec3> velocity(panels.size());
    for (std::size_t i = 0; i < panels.size(); ++i)
        velocity[i] = field.BoundMidpointVelocity(static_cast<int>(i));
    return ComputeAttachmentLine(panels, strips, velocity, sections, options);
}

/** The worst attachment-line state anywhere along the line. */
[[nodiscard]] inline AttachmentLineState WorstAttachmentLineState(const AttachmentLine& line) {
    AttachmentLineState worst = AttachmentLineState::Laminar;
    for (const AttachmentStation& station : line.Stations) {
        if (!station.Found) continue;
        if (station.State == AttachmentLineState::Turbulent) return AttachmentLineState::Turbulent;
        if (station.State == AttachmentLineState::Contaminated) worst = AttachmentLineState::Contaminated;
    }
    return worst;
}

} // namespace Aeolion::Solver
