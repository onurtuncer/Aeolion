// Solver/SurfaceFlow.h
//
// The skin-flow field of a source-panelled body: where the flow attaches,
// where it separates, and what path a fluid particle takes over the surface
// between the two. This is the geometry a boundary layer is posed on --
// every integral BL method marches from an attachment point along an
// external streamline -- so it is the prerequisite for coupling one, not
// the coupling itself.
//
// --- why the body can answer this and the lattice cannot -------------------
// A source-panelled body is a REAL closed surface carrying a real velocity
// field, so its stagnation points exist in the solution and are found
// rather than modelled. A vortex lattice is a zero-thickness camber sheet:
// the flow never stops anywhere on it, and its leading edge carries a
// square-root velocity singularity instead of a stagnation point. The
// lattice's attachment line therefore comes from the SECTION problem, which
// needs thickness -- see Solver/AttachmentLine.h. Nothing in this file
// applies to lifting panels, and it does not pretend to.
//
// --- the surface as a structured patch --------------------------------------
// Everything here works in the panelling's own (station, sector) index
// space rather than in any body-fixed coordinate like (x, azimuth). That is
// deliberate. Recovering an azimuth as atan2(z, y) assumes the surface is
// swept about +x through the origin, which is true of the fuselage, false
// of an offset duct, and silently false rather than loudly false -- the
// streamline would just curve wrongly. Lattice::SourcePanel states its own
// (StationIndex, SectorIndex) instead, and this file uses those.
//
// Working in index space needs a metric, which the grid carries:
//
//     a_i = dP/di,  a_j = dP/dj      (covariant basis, per unit index)
//     V_t = u a_i + v a_j            (u, v: contravariant velocity)
//     sqrt(g) = |a_i x a_j|          (area per unit index square)
//
// u and v are solved for once per node from the 2x2 Gram system and stored,
// so tracing is plain interpolation afterwards. They vanish exactly where
// V_t does (the Gram matrix is positive definite), which is what makes the
// critical-point search a search for a common zero of two interpolants.
//
// --- what is invariant and what is not --------------------------------------
// A critical point's TYPE is read from the Jacobian of (u, v) in index
// space, which looks like a coordinate-dependent thing to do and is not: at
// a point where the velocity vanishes, a change of surface coordinates maps
// J to M J M^-1, so the eigenvalues -- and hence node / saddle / focus and
// attachment / separation -- are invariant. The eigenVECTORS are index-space
// directions and are converted through the metric before being reported as
// 3-D directions.
//
// Sign convention: surface streamlines run along V_t. Flow ARRIVES at the
// surface and spreads out at an attachment point, so streamlines emanate
// from it and its eigenvalues have positive real parts. A separation point
// is the reverse. This is the standard Lighthill / Tobak & Peake reading of
// a skin-friction topology, applied to the inviscid surface flow that the
// boundary layer will be built on.

#pragma once

#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Solver/Solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace Aeolion::Solver {

// --- grid validity ----------------------------------------------------------
inline constexpr int MinSurfaceStations = 2; ///< A meridional derivative needs two.
inline constexpr int MinSurfaceSectors = 3;  ///< Fewer than three is not a closed ring.

// --- numerical guards --------------------------------------------------------
inline constexpr double SurfaceMetricFloor = 1e-18;   ///< |a_i x a_j|^2 below this: degenerate cell.
inline constexpr double SurfaceSpeedFloor = 1e-12;    ///< "the surface flow has stopped here"
inline constexpr double CriticalPointTolerance = 1e-10; ///< Newton residual in index units.
inline constexpr int    CriticalPointIterations = 40;

/**
 * How much a critical point must actually SPIRAL before it is called a
 * focus: the rotation rate |Im(lambda)| must reach this fraction of the
 * growth rate |Re(lambda)|.
 *
 * A negative discriminant alone is not enough. An isotropic node -- the
 * stagnation point of a sphere, and very nearly that of any blunt nose --
 * has two equal eigenvalues, which puts it exactly on the node/focus
 * boundary where the discriminant is zero. Panelling noise then pushes it
 * a hair either way, and a bare sign test reports a sphere's forward
 * stagnation point as a spiral. Requiring real rotation keeps the
 * distinction meaningful and reserves "focus" for the vortical separation
 * it is supposed to name.
 */
inline constexpr double FocusSpiralThreshold = 0.25;

// --- streamline integration --------------------------------------------------
inline constexpr double DefaultStreamlineStep = 0.25;  ///< Index units per RK4 step (a quarter cell).
inline constexpr int    DefaultStreamlineSteps = 4000;
/**
 * A trace ends when the local speed drops below this fraction of the
 * patch's mean. The threshold has to be a real fraction rather than a
 * rounding-level epsilon: a streamline running into a separation NODE
 * approaches it exponentially and never arrives, so too small a threshold
 * buys thousands of steps that creep within one panel of the node, inflate
 * the reported arc length, and bunch almost every output point into a
 * region the boundary layer will never be marched through anyway.
 */
inline constexpr double StreamlineStallFraction = 1e-3;
inline constexpr double StreamlineStepOvershoot = 2.0;  ///< A step may exceed the nominal by this much before halving.
inline constexpr int    StreamlineStepHalvings = 24;    ///< Halvings allowed before the step is accepted anyway.

/**
 * How far off a critical point its distinguished streamlines are seeded, in
 * index units. Small enough that the seed still lies in the eigen-direction's
 * linear neighbourhood, large enough that the flow there is well above the
 * panelling's own noise -- seeding a thousandth of a cell out puts the seed
 * inside the faceting error and the line sets off in an arbitrary direction.
 */
inline constexpr double CriticalPointSeedOffset = 0.05;

// ------------------------------------------------------------ one sample ----
/** The solved surface flow at one source panel's control point. */
struct SurfaceSample {
    Vec3 Point{0, 0, 0};      ///< The panel's control point (on the faceted surface).
    Vec3 Normal{0, 0, 0};     ///< Unit outward normal.
    Vec3 Tangential{0, 0, 0}; ///< V - (V.n) n, the surface flow itself [m/s].
    double Speed = 0.0;       ///< |Tangential| [m/s].

    /**
     * Pressure coefficient against the TRANSLATIONAL freestream dynamic
     * pressure. On a solve carrying body rates, an external field, or
     * propeller efflux this is still referenced to Vinf -- the conventional
     * choice, and the one the solver's own force integration uses -- so it
     * is a coefficient, not a local pressure ratio.
     */
    double Cp = 0.0;

    /** Contravariant surface velocity: V_t = Ui * a_i + Vj * a_j [index units/s]. */
    double Ui = 0.0, Vj = 0.0;
    double RootG = 0.0; ///< |a_i x a_j|, surface area per unit index square [m^2].

    /**
     * The basis here spans no area, so the contravariant decomposition has
     * no unique answer and Ui/Vj were left at zero. That zero is NOT a
     * stagnation point, and a critical-point search must skip any cell
     * touching such a node rather than report one -- a collapsed ring on
     * the body axis would otherwise read as an attachment point sitting
     * exactly where the real one often is, which is the worst possible
     * place for a false positive.
     */
    bool Degenerate = false;

    int PanelIndex = -1; ///< Back-reference into the system's source panels.
};

// -------------------------------------------------------- the surface grid ---
/**
 * One source-panelled surface as a structured (station x sector) patch,
 * carrying the solved tangential velocity and its metric at every node.
 *
 * Stations do NOT wrap (the patch has a nose end and a tail end); sectors
 * do, because a body of revolution closes on itself azimuthally.
 */
struct SurfaceGrid {
    std::string Surface;
    int Stations = 0;
    int Sectors = 0;
    std::vector<SurfaceSample> Samples; ///< Row-major, [station * Sectors + sector].
    double MeanSpeed = 0.0;             ///< Patch-average |V_t|, the scale a stall test is measured against.

    [[nodiscard]] bool Valid() const {
        return Stations >= MinSurfaceStations && Sectors >= MinSurfaceSectors &&
               Samples.size() == static_cast<std::size_t>(Stations) * static_cast<std::size_t>(Sectors);
    }

    /** Sector index wrapped into [0, Sectors); station clamped to the patch. */
    [[nodiscard]] int Wrap(int sector) const {
        const int m = sector % Sectors;
        return (m < 0) ? m + Sectors : m;
    }
    [[nodiscard]] int Clamp(int station) const { return std::clamp(station, 0, Stations - 1); }

    [[nodiscard]] const SurfaceSample& At(int station, int sector) const {
        return Samples[static_cast<std::size_t>(Clamp(station)) * static_cast<std::size_t>(Sectors) +
                       static_cast<std::size_t>(Wrap(sector))];
    }
};

namespace Detail {

/**
 * Covariant basis at a node by central differences in index space, one-sided
 * at the station ends. Sectors wrap, so their difference is always central.
 */
inline void SurfaceBasis(const SurfaceGrid& grid, int i, int j, Vec3& ai, Vec3& aj) {
    if (i <= 0) {
        ai = grid.At(1, j).Point - grid.At(0, j).Point;
    } else if (i >= grid.Stations - 1) {
        ai = grid.At(grid.Stations - 1, j).Point - grid.At(grid.Stations - 2, j).Point;
    } else {
        ai = (grid.At(i + 1, j).Point - grid.At(i - 1, j).Point) * Math::Half;
    }
    aj = (grid.At(i, j + 1).Point - grid.At(i, j - 1).Point) * Math::Half;
}

/**
 * Contravariant components of a tangential vector against a covariant basis:
 * the 2x2 Gram solve. Returns false on a degenerate cell (a ring collapsed
 * onto the axis), where the basis spans no area and the decomposition has
 * no unique answer.
 */
[[nodiscard]] inline bool Contravariant(const Vec3& ai, const Vec3& aj, const Vec3& v, double& u,
                                        double& w) {
    const double gii = Dot(ai, ai), gij = Dot(ai, aj), gjj = Dot(aj, aj);
    const double det = gii * gjj - gij * gij;
    if (!(std::fabs(det) > SurfaceMetricFloor)) {
        u = w = 0.0;
        return false;
    }
    const double bi = Dot(ai, v), bj = Dot(aj, v);
    u = (bi * gjj - bj * gij) / det;
    w = (bj * gii - bi * gij) / det;
    return true;
}

} // namespace Detail

/**
 * Assemble the named surface of a solved system into a structured grid.
 *
 * Returns an invalid grid (Valid() == false) rather than guessing when the
 * panelling does not state a complete (StationIndex, SectorIndex) topology:
 * a hole, a duplicate, or an unset index. The fuselage base cap is exactly
 * such a surface -- it is a flat disc stacked in rings at a single station,
 * so its indices are not a unique key and it is declined here by design.
 */
[[nodiscard]] inline SurfaceGrid BuildSurfaceGrid(const FlowField& field, const std::string& surface) {
    SurfaceGrid grid;
    grid.Surface = surface;
    if (!field.System) return grid;
    const std::vector<SourcePanel>& sources = field.System->Sources;

    // Extent of the index space, and a rejection of anything unindexed.
    int maxStation = -1, maxSector = -1;
    int members = 0;
    for (const SourcePanel& panel : sources) {
        if (panel.Surface != surface) continue;
        if (panel.StationIndex < 0 || panel.SectorIndex < 0) return grid; // unstructured: decline
        maxStation = std::max(maxStation, panel.StationIndex);
        maxSector = std::max(maxSector, panel.SectorIndex);
        ++members;
    }
    if (members == 0) return grid;

    const int stations = maxStation + 1;
    const int sectors = maxSector + 1;
    if (stations < MinSurfaceStations || sectors < MinSurfaceSectors) return grid;
    if (members != stations * sectors) return grid; // a hole or a duplicate

    grid.Stations = stations;
    grid.Sectors = sectors;
    grid.Samples.assign(static_cast<std::size_t>(stations) * static_cast<std::size_t>(sectors), {});

    const double vRef = field.Vinf.Norm();
    std::vector<bool> filled(grid.Samples.size(), false);

    for (std::size_t k = 0; k < sources.size(); ++k) {
        const SourcePanel& panel = sources[k];
        if (panel.Surface != surface) continue;
        const std::size_t slot = static_cast<std::size_t>(panel.StationIndex) *
                                     static_cast<std::size_t>(sectors) +
                                 static_cast<std::size_t>(panel.SectorIndex);
        if (filled[slot]) { // duplicate key: the topology is a lie, decline the whole patch
            return SurfaceGrid{};
        }
        filled[slot] = true;

        const Vec3 velocity = field.SourceSurfaceVelocity(static_cast<int>(k));
        SurfaceSample& sample = grid.Samples[slot];
        sample.Point = panel.ControlPoint;
        sample.Normal = panel.Normal;
        sample.Tangential = velocity - panel.Normal * Dot(velocity, panel.Normal);
        sample.Speed = sample.Tangential.Norm();
        const double ratio = (vRef > Math::Tiny) ? velocity.Norm() / vRef : 0.0;
        sample.Cp = 1.0 - ratio * ratio;
        sample.PanelIndex = static_cast<int>(k);
    }

    // The metric needs neighbours, so it is a second pass over a full grid.
    double speedSum = 0.0;
    for (int i = 0; i < stations; ++i) {
        for (int j = 0; j < sectors; ++j) {
            Vec3 ai, aj;
            Detail::SurfaceBasis(grid, i, j, ai, aj);
            SurfaceSample& sample =
                grid.Samples[static_cast<std::size_t>(i) * static_cast<std::size_t>(sectors) +
                             static_cast<std::size_t>(j)];
            sample.RootG = Cross(ai, aj).Norm();
            sample.Degenerate = !Detail::Contravariant(ai, aj, sample.Tangential, sample.Ui, sample.Vj);
            speedSum += sample.Speed;
        }
    }
    grid.MeanSpeed = speedSum / static_cast<double>(grid.Samples.size());
    return grid;
}

// ------------------------------------------------- interpolation in index space
/** A surface state at a non-integer (station, sector) coordinate. */
struct SurfaceState {
    double Station = 0.0, Sector = 0.0; ///< Index-space position.
    Vec3 Point{0, 0, 0};
    Vec3 Normal{0, 0, 0};
    Vec3 Tangential{0, 0, 0};
    double Speed = 0.0;
    double Cp = 0.0;
    double Ui = 0.0, Vj = 0.0;
    double RootG = 0.0;
};

/**
 * Bilinear interpolation over the patch. Sectors wrap; stations clamp, so a
 * streamline that runs off the nose or tail is held at the boundary and the
 * caller's own range test ends the trace.
 */
[[nodiscard]] inline SurfaceState SampleSurface(const SurfaceGrid& grid, double station, double sector) {
    SurfaceState state;
    state.Station = station;
    state.Sector = sector;
    if (!grid.Valid()) return state;

    const double clampedStation = std::clamp(station, 0.0, static_cast<double>(grid.Stations - 1));
    const int i0 = std::min(static_cast<int>(std::floor(clampedStation)), grid.Stations - 2);
    const double fi = clampedStation - i0;
    const int j0 = static_cast<int>(std::floor(sector));
    const double fj = sector - j0;

    const double w00 = (1.0 - fi) * (1.0 - fj), w10 = fi * (1.0 - fj);
    const double w01 = (1.0 - fi) * fj, w11 = fi * fj;

    const SurfaceSample& s00 = grid.At(i0, j0);
    const SurfaceSample& s10 = grid.At(i0 + 1, j0);
    const SurfaceSample& s01 = grid.At(i0, j0 + 1);
    const SurfaceSample& s11 = grid.At(i0 + 1, j0 + 1);

    const auto blend = [&](auto member) {
        return s00.*member * w00 + s10.*member * w10 + s01.*member * w01 + s11.*member * w11;
    };
    state.Point = s00.Point * w00 + s10.Point * w10 + s01.Point * w01 + s11.Point * w11;
    state.Normal = (s00.Normal * w00 + s10.Normal * w10 + s01.Normal * w01 + s11.Normal * w11).Normalized();
    state.Tangential =
        s00.Tangential * w00 + s10.Tangential * w10 + s01.Tangential * w01 + s11.Tangential * w11;
    state.Speed = state.Tangential.Norm();
    state.Cp = blend(&SurfaceSample::Cp);
    state.Ui = blend(&SurfaceSample::Ui);
    state.Vj = blend(&SurfaceSample::Vj);
    state.RootG = blend(&SurfaceSample::RootG);
    return state;
}

// ------------------------------------------------------------- streamlines ---
/** One point along a traced surface streamline. */
struct StreamlinePoint {
    Vec3 Point{0, 0, 0};
    double ArcLength = 0.0; ///< s from the seed [m], always increasing regardless of trace direction.
    double Speed = 0.0;     ///< Edge speed U_e at this station [m/s].
    double Cp = 0.0;
    double Divergence = 0.0; ///< Metric factor h(s), see SurfaceStreamline::Points.
    double Station = 0.0, Sector = 0.0;
};

/** Why a trace stopped -- the caller usually needs to know. */
enum class StreamlineExit {
    StepLimit,     ///< Ran out of steps; the line is incomplete.
    LeftPatch,     ///< Reached the nose or tail boundary of the station range.
    Stalled,       ///< Surface speed collapsed: converged onto a critical point.
    InvalidGrid,   ///< Nothing to trace.
};

/**
 * A traced surface streamline: the path, and along it the two quantities an
 * integral boundary layer marches on.
 *
 * The first is the edge speed U_e(s). The second is Divergence, the
 * streamline-spreading metric h(s) of the axisymmetric analogue: the
 * momentum-integral equation on a three-dimensional surface differs from
 * its two-dimensional form only by
 *
 *     dtheta/ds + (H + 2) theta/U_e dU_e/ds + theta/h dh/ds = cf/2,
 *
 * so a BL march needs h as much as it needs U_e. It is integrated along the
 * line from the surface divergence of the unit flow direction,
 *
 *     d(ln h)/ds = div_s(e),   e = V_t / |V_t|,
 *
 * normalized to h = 1 at the seed. Where streamlines fan apart h grows and
 * the layer is thinned; where they converge h collapses, which is the
 * signature of a separation line and the reason h is worth carrying rather
 * than assuming.
 */
struct SurfaceStreamline {
    std::vector<StreamlinePoint> Points;
    StreamlineExit Exit = StreamlineExit::InvalidGrid;
    double Length = 0.0; ///< Total arc length [m].
    bool Upstream = false; ///< Traced against the flow (see StreamlineOptions).
};

/** How a streamline trace is run. */
struct StreamlineOptions {
    /**
     * Trace against the flow rather than with it. Downstream answers "where
     * does this go"; UPSTREAM answers "where did the flow at this point
     * come from, and how far has it run" -- which is the question a
     * boundary layer at a given surface station actually asks, since its
     * state is an integral over everything upstream of it.
     */
    bool Upstream = false;
    double Step = DefaultStreamlineStep; ///< Index units per RK4 step.
    int MaxSteps = DefaultStreamlineSteps;
    /** Stop when the local speed falls below this fraction of the patch mean. */
    double StallFraction = StreamlineStallFraction;
};

namespace Detail {

// --- differentiating the surface flow ----------------------------------------
// Anything that needs a DERIVATIVE of the surface flow -- classifying a
// critical point, or measuring how fast streamlines fan apart -- takes it
// from the tangential velocity vector in a local orthonormal frame, never
// from the contravariant components Ui/Vj.
//
// The reason is that Ui and Vj are not smooth fields even where the flow
// is. On a body of revolution the azimuthal basis vector shrinks like the
// local radius, so Vj carries a 1/r factor that varies by tens of percent
// across a single cell near a nose. Interpolating it is fine (it is exact
// at the nodes and that is all a streamline integrator needs), but
// DIFFERENCING it mixes the chart's own distortion into the answer: the
// sphere's stagnation node, which is exactly isotropic, comes out with a
// 2:1 eigenvalue split.
//
// The tangential velocity itself has no such problem -- it is a smooth
// vector field sampled at smooth points -- so it is what gets differenced,
// and the frame it is differenced in is orthonormal, which makes the
// resulting 2x2 Jacobian a genuinely physical strain rate in 1/s.

/** A local orthonormal tangent frame, with the index-space steps that walk along it. */
struct LocalFrame {
    bool Valid = false;
    Vec3 E1{0, 0, 0}, E2{0, 0, 0}; ///< Orthonormal, tangent to the surface.
    double Step = 0.0;             ///< Physical half-step used for the differences [m].
    double I1 = 0.0, J1 = 0.0;     ///< Index displacement of one Step along E1.
    double I2 = 0.0, J2 = 0.0;     ///< Index displacement of one Step along E2.
};

[[nodiscard]] inline LocalFrame MakeLocalFrame(const SurfaceGrid& grid, double station, double sector) {
    LocalFrame frame;
    if (!grid.Valid()) return frame;

    Vec3 ai, aj;
    SurfaceBasis(grid, static_cast<int>(std::lround(station)), static_cast<int>(std::lround(sector)), ai,
                 aj);
    const SurfaceState state = SampleSurface(grid, station, sector);
    const Vec3 normal = state.Normal;

    Vec3 e1 = ai - normal * Dot(ai, normal);
    if (!(e1.Norm() > Math::Tiny)) return frame;
    e1 = e1.Normalized();
    const Vec3 e2 = Cross(normal, e1).Normalized();
    if (!(e2.Norm() > Math::Half)) return frame;

    // Half a meridional cell. The meridional spacing is the well-behaved one
    // -- the azimuthal spacing collapses toward a nose -- so it sets the
    // scale, and the azimuthal step follows from the metric.
    const double meridional = ai.Norm();
    if (!(meridional > Math::Tiny)) return frame;
    frame.Step = Math::Half * meridional;

    if (!Contravariant(ai, aj, e1 * frame.Step, frame.I1, frame.J1)) return frame;
    if (!Contravariant(ai, aj, e2 * frame.Step, frame.I2, frame.J2)) return frame;
    frame.E1 = e1;
    frame.E2 = e2;
    frame.Valid = true;
    return frame;
}

/**
 * The 2x2 Jacobian of the tangential flow in the frame's own axes, by
 * central differences over one physical step. `unit` differentiates the
 * flow DIRECTION instead of the velocity, whose trace is the surface
 * divergence that drives streamline spreading.
 *
 *     J = [ d(V.E1)/ds1  d(V.E1)/ds2 ]      [1/s], or [1/m] when unit
 *         [ d(V.E2)/ds1  d(V.E2)/ds2 ]
 */
inline void PhysicalJacobian(const SurfaceGrid& grid, const LocalFrame& frame, double station,
                             double sector, bool unit, double J[4]) {
    J[0] = J[1] = J[2] = J[3] = 0.0;
    if (!frame.Valid) return;

    const auto sampleVector = [&](double s, double t) -> Vec3 {
        const SurfaceState state = SampleSurface(grid, s, t);
        if (!unit) return state.Tangential;
        return (state.Speed > SurfaceSpeedFloor) ? state.Tangential * (1.0 / state.Speed) : Vec3(0, 0, 0);
    };

    const Vec3 along1 = sampleVector(station + frame.I1, sector + frame.J1) -
                        sampleVector(station - frame.I1, sector - frame.J1);
    const Vec3 along2 = sampleVector(station + frame.I2, sector + frame.J2) -
                        sampleVector(station - frame.I2, sector - frame.J2);

    const double inv = 1.0 / (2.0 * frame.Step);
    J[0] = Dot(along1, frame.E1) * inv;
    J[1] = Dot(along2, frame.E1) * inv;
    J[2] = Dot(along1, frame.E2) * inv;
    J[3] = Dot(along2, frame.E2) * inv;
}

/** Surface divergence of the unit flow direction: d(ln h)/ds along a streamline. */
[[nodiscard]] inline double UnitFlowDivergence(const SurfaceGrid& grid, double station, double sector) {
    const LocalFrame frame = MakeLocalFrame(grid, station, sector);
    if (!frame.Valid) return 0.0;
    double J[4];
    PhysicalJacobian(grid, frame, station, sector, true, J);
    return J[0] + J[3];
}

} // namespace Detail

/**
 * Trace a surface streamline from an index-space seed. The integration runs
 * in index space with RK4 on the contravariant velocity, so a step is a
 * fixed fraction of a CELL rather than of a metre -- the natural step for a
 * mesh whose cells vary in size along the body.
 */
[[nodiscard]] inline SurfaceStreamline TraceSurfaceStreamline(const SurfaceGrid& grid, double seedStation,
                                                              double seedSector,
                                                              const StreamlineOptions& options = {}) {
    SurfaceStreamline line;
    line.Upstream = options.Upstream;
    if (!grid.Valid()) return line;

    const double sign = options.Upstream ? -1.0 : 1.0;
    const double stallSpeed = std::max(grid.MeanSpeed * options.StallFraction, SurfaceSpeedFloor);
    const double topStation = static_cast<double>(grid.Stations - 1);

    // Index-space velocity, already carrying the trace direction.
    const auto derivative = [&](double station, double sector, double& dI, double& dJ) {
        const SurfaceState state = SampleSurface(grid, station, sector);
        dI = sign * state.Ui;
        dJ = sign * state.Vj;
    };

    double station = seedStation, sector = seedSector;
    double arc = 0.0, logH = 0.0;
    line.Exit = StreamlineExit::StepLimit;

    for (int step = 0; step <= options.MaxSteps; ++step) {
        const SurfaceState state = SampleSurface(grid, station, sector);

        StreamlinePoint point;
        point.Point = state.Point;
        point.ArcLength = arc;
        point.Speed = state.Speed;
        point.Cp = state.Cp;
        point.Divergence = std::exp(logH);
        point.Station = station;
        point.Sector = sector;
        line.Points.push_back(point);

        if (state.Speed < stallSpeed) { line.Exit = StreamlineExit::Stalled; break; }
        if (station <= 0.0 || station >= topStation) {
            // The seed itself may sit on the boundary; only a LATER arrival ends the trace.
            if (step > 0) { line.Exit = StreamlineExit::LeftPatch; break; }
        }
        if (step == options.MaxSteps) break;

        // RK4 in index space. The parameter is not time: it is scaled so a
        // step covers `Step` index units at the LOCAL flow speed, which
        // keeps the step a constant fraction of a cell wherever the mesh is
        // fine or coarse.
        //
        // That scaling alone is not enough near a critical point, and this
        // is the one place the integrator can go badly wrong. There the
        // field is xi_dot = lambda xi, so dividing by a velocity that is
        // itself proportional to the distance from the zero makes dt grow
        // without bound -- seeded a thousandth of a cell out, the nominal
        // step is hundreds of e-folding times, and RK4 answers with a number
        // like 1e40. So the step is additionally capped by DISPLACEMENT,
        // halving until the move is within a small multiple of one nominal
        // step. Away from critical points the cap never binds and the
        // integration is plain RK4.
        double k1i, k1j, k2i, k2j, k3i, k3j, k4i, k4j;
        derivative(station, sector, k1i, k1j);
        const double scale = std::hypot(k1i, k1j);
        if (!(scale > SurfaceSpeedFloor)) { line.Exit = StreamlineExit::Stalled; break; }

        double dt = options.Step / scale; // index units per step / (index units per second)
        double stepI = 0.0, stepJ = 0.0;
        for (int trial = 0; trial < StreamlineStepHalvings; ++trial) {
            derivative(station + Math::Half * dt * k1i, sector + Math::Half * dt * k1j, k2i, k2j);
            derivative(station + Math::Half * dt * k2i, sector + Math::Half * dt * k2j, k3i, k3j);
            derivative(station + dt * k3i, sector + dt * k3j, k4i, k4j);
            stepI = dt / 6.0 * (k1i + 2.0 * k2i + 2.0 * k3i + k4i);
            stepJ = dt / 6.0 * (k1j + 2.0 * k2j + 2.0 * k3j + k4j);
            if (std::hypot(stepI, stepJ) <= StreamlineStepOvershoot * options.Step) break;
            dt *= Math::Half;
        }

        const double nextStation = station + stepI;
        const double nextSector = sector + stepJ;
        const SurfaceState next = SampleSurface(grid, nextStation, nextSector);

        // Physical arc length of the step, and the spreading metric over it.
        const double ds = (next.Point - state.Point).Norm();
        arc += ds;
        logH += Math::Half *
                (Detail::UnitFlowDivergence(grid, station, sector) +
                 Detail::UnitFlowDivergence(grid, nextStation, nextSector)) *
                ds * sign;

        station = nextStation;
        sector = nextSector;
    }

    line.Length = arc;
    return line;
}

// -------------------------------------------------------- critical points ----
/** What a zero of the surface flow is. */
enum class CriticalPointType {
    AttachmentNode,  ///< Streamlines emanate: the flow arrives here and spreads.
    SeparationNode,  ///< Streamlines converge: the flow leaves the surface here.
    Saddle,          ///< Two in, two out; its separatrices are the dividing lines.
    AttachmentFocus, ///< Spiral source.
    SeparationFocus, ///< Spiral sink -- the classic vortical separation.
    Degenerate,      ///< Determinant too near zero to classify honestly.
};

/** A zero of the surface velocity field, located and classified. */
struct CriticalPoint {
    Vec3 Point{0, 0, 0};
    double Station = 0.0, Sector = 0.0; ///< Index-space position.
    CriticalPointType Type = CriticalPointType::Degenerate;

    double Trace = 0.0;       ///< tr(J), [1/s] -- sign separates attachment from separation.
    double Determinant = 0.0; ///< det(J); negative is a saddle.
    double Discriminant = 0.0; ///< tr^2 - 4 det; negative is a focus.

    /**
     * The real eigen-directions as UNIT 3-D vectors tangent to the surface,
     * ordered fast (larger |eigenvalue|) first. At a node these are the
     * principal spreading directions; at a saddle they are the separatrices,
     * which are the dividing streamlines. Empty for a focus, which has none.
     */
    std::vector<Vec3> PrincipalDirections;
    std::vector<double> Eigenvalues; ///< Aligned with PrincipalDirections [1/s].

    [[nodiscard]] bool IsAttachment() const {
        return Type == CriticalPointType::AttachmentNode || Type == CriticalPointType::AttachmentFocus;
    }
    [[nodiscard]] bool IsSeparation() const {
        return Type == CriticalPointType::SeparationNode || Type == CriticalPointType::SeparationFocus;
    }
};

namespace Detail {

/** Bilinear (u, v) and their index-space Jacobian inside cell (i0, j0). */
inline void CellField(const SurfaceGrid& grid, int i0, int j0, double fi, double fj, double& u, double& v,
                      double J[4]) {
    const SurfaceSample& s00 = grid.At(i0, j0);
    const SurfaceSample& s10 = grid.At(i0 + 1, j0);
    const SurfaceSample& s01 = grid.At(i0, j0 + 1);
    const SurfaceSample& s11 = grid.At(i0 + 1, j0 + 1);

    u = s00.Ui * (1 - fi) * (1 - fj) + s10.Ui * fi * (1 - fj) + s01.Ui * (1 - fi) * fj + s11.Ui * fi * fj;
    v = s00.Vj * (1 - fi) * (1 - fj) + s10.Vj * fi * (1 - fj) + s01.Vj * (1 - fi) * fj + s11.Vj * fi * fj;

    J[0] = (s10.Ui - s00.Ui) * (1 - fj) + (s11.Ui - s01.Ui) * fj; // du/di
    J[1] = (s01.Ui - s00.Ui) * (1 - fi) + (s11.Ui - s10.Ui) * fi; // du/dj
    J[2] = (s10.Vj - s00.Vj) * (1 - fj) + (s11.Vj - s01.Vj) * fj; // dv/di
    J[3] = (s01.Vj - s00.Vj) * (1 - fi) + (s11.Vj - s10.Vj) * fi; // dv/dj
}

} // namespace Detail

/**
 * Every zero of the surface velocity field inside the patch, classified.
 *
 * A cell is a candidate when both contravariant components change sign
 * across its corners; a 2-D Newton on the bilinear interpolant then locates
 * the zero, and it is kept only if it lands inside the cell. Coarse
 * panelling can merge two nearby critical points into one or miss a pair
 * entirely -- the usual resolution caveat, and the reason the returned
 * eigenvalues are worth looking at rather than just the count.
 */
[[nodiscard]] inline std::vector<CriticalPoint> FindCriticalPoints(const SurfaceGrid& grid) {
    std::vector<CriticalPoint> found;
    if (!grid.Valid()) return found;

    for (int i0 = 0; i0 + 1 < grid.Stations; ++i0) {
        for (int j0 = 0; j0 < grid.Sectors; ++j0) {
            const SurfaceSample* corner[4] = {&grid.At(i0, j0), &grid.At(i0 + 1, j0), &grid.At(i0, j0 + 1),
                                              &grid.At(i0 + 1, j0 + 1)};
            bool uPositive = false, uNegative = false, vPositive = false, vNegative = false;
            bool degenerate = false;
            for (const SurfaceSample* c : corner) {
                uPositive |= (c->Ui > 0.0); uNegative |= (c->Ui < 0.0);
                vPositive |= (c->Vj > 0.0); vNegative |= (c->Vj < 0.0);
                degenerate |= c->Degenerate;
            }
            // A degenerate corner's zeroed components are an absence of an
            // answer, not a stagnation point (see SurfaceSample::Degenerate).
            if (degenerate) continue;
            if (!(uPositive && uNegative && vPositive && vNegative)) continue;

            // Newton from the cell centre on the bilinear interpolant.
            double fi = Math::Half, fj = Math::Half;
            double u = 0.0, v = 0.0, J[4] = {0, 0, 0, 0};
            bool converged = false;
            for (int iter = 0; iter < CriticalPointIterations; ++iter) {
                Detail::CellField(grid, i0, j0, fi, fj, u, v, J);
                if (std::fabs(u) + std::fabs(v) < CriticalPointTolerance * std::max(grid.MeanSpeed, 1.0)) {
                    converged = true;
                    break;
                }
                const double det = J[0] * J[3] - J[1] * J[2];
                if (!(std::fabs(det) > SurfaceMetricFloor)) break;
                fi -= (J[3] * u - J[1] * v) / det;
                fj -= (-J[2] * u + J[0] * v) / det;
                if (fi < -0.5 || fi > 1.5 || fj < -0.5 || fj > 1.5) break; // wandered out; not this cell's
            }
            if (!converged || fi < 0.0 || fi > 1.0 || fj < 0.0 || fj > 1.0) continue;

            CriticalPoint point;
            point.Station = i0 + fi;
            point.Sector = j0 + fj;
            point.Point = SampleSurface(grid, point.Station, point.Sector).Point;

            // The index-space Jacobian above LOCATED the zero; it does not
            // get to classify it. Strain rates come from the tangential
            // velocity differenced in a local orthonormal frame, for the
            // reason set out at Detail::MakeLocalFrame.
            const Detail::LocalFrame frame = Detail::MakeLocalFrame(grid, point.Station, point.Sector);
            if (!frame.Valid) continue;
            double P[4];
            Detail::PhysicalJacobian(grid, frame, point.Station, point.Sector, false, P);

            point.Trace = P[0] + P[3];
            point.Determinant = P[0] * P[3] - P[1] * P[2];
            point.Discriminant = point.Trace * point.Trace - 4.0 * point.Determinant;

            // A focus needs REAL rotation, not just a negative discriminant
            // (see FocusSpiralThreshold).
            const double growth = std::fabs(Math::Half * point.Trace);
            const double rotation =
                (point.Discriminant < 0.0) ? Math::Half * std::sqrt(-point.Discriminant) : 0.0;
            const bool spirals = rotation > FocusSpiralThreshold * std::max(growth, Math::Tiny);

            const double scale = std::max(std::fabs(point.Trace), std::sqrt(std::fabs(point.Determinant)));
            if (!(scale > Math::Tiny)) {
                point.Type = CriticalPointType::Degenerate;
            } else if (point.Determinant < 0.0) {
                point.Type = CriticalPointType::Saddle;
            } else if (spirals) {
                point.Type = (point.Trace > 0.0) ? CriticalPointType::AttachmentFocus
                                                 : CriticalPointType::SeparationFocus;
            } else {
                point.Type = (point.Trace > 0.0) ? CriticalPointType::AttachmentNode
                                                 : CriticalPointType::SeparationNode;
            }

            // A node whose discriminant came out slightly negative is a STAR
            // node -- two equal eigenvalues, and genuinely no distinguished
            // direction, which is the correct answer for a sphere rather
            // than a failure to find one. Both frame axes are reported as
            // principal, so a caller still gets four seeds to trace from.
            if (point.Discriminant < 0.0 && !spirals &&
                point.Type != CriticalPointType::Degenerate) {
                point.PrincipalDirections = {frame.E1, frame.E2};
                point.Eigenvalues = {Math::Half * point.Trace, Math::Half * point.Trace};
            }

            // Real eigen-directions. The frame is orthonormal, so an
            // eigenvector's components ARE its 3-D direction -- no metric to
            // push it through.
            if (point.Discriminant >= 0.0 && point.Type != CriticalPointType::Degenerate) {
                const double root = std::sqrt(point.Discriminant);
                double lambda[2] = {Math::Half * (point.Trace + root), Math::Half * (point.Trace - root)};
                if (std::fabs(lambda[1]) > std::fabs(lambda[0])) std::swap(lambda[0], lambda[1]);
                for (double value : lambda) {
                    // (P - lambda I) e = 0 -- take whichever row is better conditioned.
                    double e1 = P[1], e2 = value - P[0];
                    if (std::fabs(P[1]) + std::fabs(value - P[0]) <
                        std::fabs(value - P[3]) + std::fabs(P[2])) {
                        e1 = value - P[3];
                        e2 = P[2];
                    }
                    const Vec3 direction = (frame.E1 * e1 + frame.E2 * e2).Normalized();
                    if (direction.Norm() > Math::Half) { // a real direction, not a collapsed one
                        point.PrincipalDirections.push_back(direction);
                        point.Eigenvalues.push_back(value);
                    }
                }
            }
            found.push_back(point);
        }
    }
    return found;
}

// ---------------------------------------------- attachment / separation lines
/**
 * The distinguished streamlines through a critical point: the ones leaving
 * (or arriving) along its principal directions.
 *
 * These are not ordinary streamlines. At a node every neighbouring
 * streamline eventually aligns with the SLOW eigenvector, and only the two
 * along the fast one are distinct -- those two are the attachment (or
 * separation) line. At a saddle the four separatrices divide the surface
 * into regions no streamline crosses, which is what makes a saddle's
 * separatrix the honest definition of a dividing line.
 *
 * Traced away from an attachment point and toward a separation point, so
 * every returned line runs in the direction the physical flow does.
 */
[[nodiscard]] inline std::vector<SurfaceStreamline> TraceFromCriticalPoint(
    const SurfaceGrid& grid, const CriticalPoint& point, const StreamlineOptions& options = {}) {
    std::vector<SurfaceStreamline> lines;
    if (!grid.Valid()) return lines;

    StreamlineOptions local = options;
    // Streamlines run OUT of an attachment point and IN to a separation
    // point, so the latter is reached by tracing upstream.
    local.Upstream = point.IsSeparation();

    // Index-space eigen-directions, recovered by projecting the reported
    // 3-D directions back onto the local basis.
    Vec3 ai, aj;
    Detail::SurfaceBasis(grid, static_cast<int>(std::lround(point.Station)),
                         static_cast<int>(std::lround(point.Sector)), ai, aj);

    for (const Vec3& direction : point.PrincipalDirections) {
        double di = 0.0, dj = 0.0;
        if (!Detail::Contravariant(ai, aj, direction, di, dj)) continue;
        const double norm = std::hypot(di, dj);
        if (!(norm > Math::Tiny)) continue;
        di /= norm;
        dj /= norm;
        for (double side : {1.0, -1.0}) {
            const double seedStation = point.Station + side * di * CriticalPointSeedOffset;
            const double seedSector = point.Sector + side * dj * CriticalPointSeedOffset;
            if (seedStation < 0.0 || seedStation > static_cast<double>(grid.Stations - 1)) continue;
            lines.push_back(TraceSurfaceStreamline(grid, seedStation, seedSector, local));
        }
    }
    return lines;
}

/**
 * The complete skin-flow picture of one surface at one flight condition.
 *
 * `AttachesUpstream` is the case that would otherwise read as "no
 * attachment point exists": at zero incidence the stagnation point of a
 * closed body of revolution sits exactly on the nose apex, which is not a
 * panel control point but the collapsed forward edge of the first ring --
 * upstream of every sample, so no interior zero is found. Rather than
 * inventing one, the flag says so, and NoseAttachment holds the apex. Any
 * incidence at all moves the stagnation point onto the surface proper,
 * where it is found normally.
 */
struct SurfaceFlowTopology {
    std::vector<CriticalPoint> CriticalPoints;
    bool AttachesUpstream = false;
    Vec3 NoseAttachment{0, 0, 0}; ///< Meaningful only when AttachesUpstream.
    std::string Surface;
};

/** Locate and classify the whole skin-flow topology of one surface. */
[[nodiscard]] inline SurfaceFlowTopology AnalyzeSurfaceFlow(const SurfaceGrid& grid) {
    SurfaceFlowTopology topology;
    topology.Surface = grid.Surface;
    if (!grid.Valid()) return topology;
    topology.CriticalPoints = FindCriticalPoints(grid);

    const bool hasAttachment =
        std::ranges::any_of(topology.CriticalPoints, [](const CriticalPoint& c) { return c.IsAttachment(); });
    if (hasAttachment) return topology;

    // No interior attachment. If the whole first ring flows aft, the flow
    // entered through the patch's forward boundary.
    bool allAft = true;
    Vec3 apex(0, 0, 0);
    for (int j = 0; j < grid.Sectors; ++j) {
        if (!(grid.At(0, j).Ui > 0.0)) { allAft = false; break; }
        apex = apex + grid.At(0, j).Point;
    }
    if (allAft) {
        topology.AttachesUpstream = true;
        // The forward boundary's centroid. On a closed nose the first ring's
        // control points are a small circle about the apex, so their mean is
        // the apex itself to the panelling's own accuracy.
        topology.NoseAttachment = apex * (1.0 / static_cast<double>(grid.Sectors));
    }
    return topology;
}

/** The most upstream attachment point, or null when there is none inside the patch. */
[[nodiscard]] inline const CriticalPoint* PrimaryAttachment(const SurfaceFlowTopology& topology) {
    const CriticalPoint* best = nullptr;
    for (const CriticalPoint& point : topology.CriticalPoints) {
        if (!point.IsAttachment()) continue;
        if (!best || point.Station < best->Station) best = &point;
    }
    return best;
}

/** The most upstream separation point, or null when there is none inside the patch. */
[[nodiscard]] inline const CriticalPoint* PrimarySeparation(const SurfaceFlowTopology& topology) {
    const CriticalPoint* best = nullptr;
    for (const CriticalPoint& point : topology.CriticalPoints) {
        if (!point.IsSeparation()) continue;
        if (!best || point.Station < best->Station) best = &point;
    }
    return best;
}

} // namespace Aeolion::Solver
