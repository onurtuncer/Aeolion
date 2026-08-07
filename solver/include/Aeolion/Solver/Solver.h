// Solver/Solver.h
//
// Minimal steady vortex lattice method (VLM) for a single lifting surface.
//
// Convention (body/wind-ish axes):
//   x : downstream / aft
//   y : spanwise, right positive
//   z : up
//
// Each spanwise panel is represented by ONE horseshoe vortex whose bound
// segment sits on the local quarter-chord line; the control point sits on
// the local three-quarter-chord line (classic Katz & Plotkin VLM1). This
// is a "flat" lattice: camber is not modeled, but sweep, taper, dihedral,
// and linear washout/twist are all captured through panel geometry and
// local normal-vector orientation.
//
// Forces are computed with the NEAR-FIELD method: at each bound-vortex
// midpoint we sum the velocity induced by every OTHER horseshoe (all of
// them, including their trailing legs) plus freestream, then apply
// Kutta-Joukowski to that local segment. This gives lift AND induced drag
// directly, without a separate Trefftz-plane integration.
//
// The data model lives outside this header: Vec3 in Aeolion/Math, Panel in
// Aeolion/Lattice (shared vocabulary, since PanelBuilder produces panels
// and the viewer draws them), and the result/parameter structs in this
// header's siblings under Aeolion/Solver. This header owns the algorithms
// and the LAPACK-backed dense solver. The dense solve uses LAPACK
// (dgetrf/dgetrs via the Fortran ABI); link against a LAPACK provider such
// as OpenBLAS.

#pragma once
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <map>
#include <functional>
#include <numbers>
#include <ranges>
#include <mdspan>

#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Lattice/Panel.h"
#include "Aeolion/Lattice/SourcePanel.h"
#include "Aeolion/Solver/SourceInfluence.h"
#include "Aeolion/Solver/WingParams.h"
#include "Aeolion/Solver/FreestreamConditions.h"
#include "Aeolion/Solver/ReferenceGeometry.h"
#include "Aeolion/Solver/StationResult.h"
#include "Aeolion/Solver/SolveResult.h"
#include "Aeolion/Solver/StabilityDerivatives.h"

// The dense solve calls LAPACK's double-precision general LU routines
// directly through their Fortran ABI (dgetrf_ / dgetrs_). We bind these
// rather than the LAPACKE C wrapper so the only build requirement is a
// LAPACK library (OpenBLAS, reference LAPACK, MKL, ...) -- no lapacke.h,
// which several LAPACK distributions (incl. the common Windows/vcpkg
// OpenBLAS build) don't ship. These prototypes must sit at global scope.
extern "C" {
    void dgetrf_(const int* m, const int* n, double* a, const int* lda,
                 int* ipiv, int* info);
    void dgetrs_(const char* trans, const int* n, const int* nrhs,
                 const double* a, const int* lda, const int* ipiv,
                 double* b, const int* ldb, int* info);
}

namespace Aeolion::Solver {

using Math::Half;
using Math::Two;
using Math::DegToRad;
using Math::RadToDeg;

// --- tuning / tolerance constants (no bare literals in the algorithms) ----
inline constexpr double VortexCoreCutoff    = 1e-8;   // cross-product^2 cutoff regularizing a filament's own singularity
inline constexpr double FilamentEndpointEps = 1e-10;  // "point is on the filament endpoint" guard
inline constexpr double InvFourPi           = 1.0 / (4.0 * std::numbers::pi); // Biot-Savart kernel coefficient
inline constexpr int    SemiSpanCount       = 2;      // a full span is two semi-spans
inline constexpr double DefaultTrailSpanFactor = 50.0; // auto trailing-leg length = this * span
inline constexpr double DenormFloor         = 1e-300; // keeps the matrix-scale normalization finite
inline constexpr double CoeffDenomEps       = 1e-12;  // guard before dividing by dynamic-pressure * area
inline constexpr double GeometryEps         = 1e-9;   // "practically zero" length/area guard
inline constexpr double UnitFallbackLength  = 1.0;    // reference length used when none is derivable
inline constexpr double DefaultAlphaStepDeg = 0.5;    // central-difference step, angle of attack
inline constexpr double DefaultBetaStepDeg  = 0.5;    // central-difference step, sideslip
inline constexpr double DefaultRateStepFraction = 0.002; // central-difference rate step, as a fraction of Vinf/length

// --------------------------------------------- finite vortex segment -----
/**
 * Biot-Savart induced velocity at point P from a straight vortex filament
 * running P1 -> P2 with circulation strength gamma. A small viscous-core
 * cutoff regularizes the singularity on the filament itself; a nonzero
 * coreRadius [m] additionally applies a finite Vatistas-style vortex
 * core, v ~ v_thin * h^2/(h^2 + rc^2) with h the perpendicular distance
 * -- required once curved wake legs can thread between panels, where
 * grazing an unregularized filament is numerically vicious.
 */
[[nodiscard]] inline Vec3 SegmentVelocity(const Vec3& P, const Vec3& P1, const Vec3& P2, double gamma,
                                          double coreRadius = 0.0) {
    Vec3 r1 = P - P1;
    Vec3 r2 = P - P2;
    Vec3 r0 = P2 - P1;

    Vec3 c = Cross(r1, r2);
    double cmag2 = Dot(c, c);

    double r1n = r1.Norm();
    double r2n = r2.Norm();

    if (cmag2 < VortexCoreCutoff || r1n < FilamentEndpointEps || r2n < FilamentEndpointEps) {
        return {0, 0, 0};   // point on (or essentially on) the filament -> no contribution
    }

    double K = gamma * InvFourPi / cmag2 * (Dot(r0, r1) / r1n - Dot(r0, r2) / r2n);
    if (coreRadius > 0.0) {
        const double r0n2 = Dot(r0, r0);
        if (r0n2 > VortexCoreCutoff) {
            const double h2 = cmag2 / r0n2; // squared perpendicular distance to the filament line
            K *= h2 / (h2 + coreRadius * coreRadius);
        }
    }
    return c * K;
}

/** Linearly-tapered local chord at absolute spanwise coordinate yAbs. */
[[nodiscard]] inline double ChordAt(const WingParams& w, double yAbs) {
    double eta = yAbs / (w.Span * Half);
    return w.RootChord + (w.TipChord - w.RootChord) * eta;
}

/**
 * Builds the panel lattice for a single trapezoidal, linearly-twisted,
 * swept, dihedral wing. Trailing legs are aligned with the global +x axis
 * (standard small-to-moderate AoA simplification).
 */
[[nodiscard]] inline std::vector<Panel> BuildWing(const WingParams& w) {
    if (w.NPanelsSemiSpan < 1) throw std::invalid_argument("NPanelsSemiSpan must be >= 1");

    int N = SemiSpanCount * w.NPanelsSemiSpan;
    std::vector<double> yStations(N + 1);
    double halfSpan = w.Span * Half;

    for (int i = 0; i <= N; ++i) {
        double u = -1.0 + Two * i / N; // -1..1
        if (w.CosineSpacing) {
            // cosine clustering toward tips: theta in [0,pi], y = -cos(theta)
            double theta = (u + 1.0) * Half * std::numbers::pi; // 0..pi
            yStations[i] = -halfSpan * std::cos(theta);
        } else {
            yStations[i] = halfSpan * u;
        }
    }

    double sweep = DegToRad(w.SweepQuarterChordDeg);
    double dihedral = DegToRad(w.DihedralDeg);

    auto qcX = [&](double y) { return std::fabs(y) * std::tan(sweep); };
    auto zOf = [&](double y) { return std::fabs(y) * std::tan(dihedral); };
    auto twistOf = [&](double y) {
        double eta = std::fabs(y) / halfSpan;
        return DegToRad(w.TwistTipDeg) * eta;
    };

    std::vector<Panel> panels;
    panels.reserve(N);

    for (int i = 0; i < N; ++i) {
        double yA = yStations[i];
        double yB = yStations[i + 1];
        if (yA > yB) std::swap(yA, yB); // keep A.y < B.y

        double cA = ChordAt(w, std::fabs(yA));
        double cB = ChordAt(w, std::fabs(yB));

        Vec3 A(qcX(yA), yA, zOf(yA));
        Vec3 B(qcX(yB), yB, zOf(yB));

        double yc = Half * (yA + yB);
        double cc = ChordAt(w, std::fabs(yc));
        double eps = twistOf(yc);

        // Twist is a physical rotation of the chordline about the local
        // quarter-chord point -- rotate the quarter->three-quarter-chord
        // offset by eps (about the spanwise axis) rather than just tilting
        // the boundary-condition normal. This keeps the panel geometry
        // itself consistent with what a real (or mesh-derived) twisted
        // wing looks like, not just a small-angle linearization.
        Vec3 qcToCp = RotateAboutY(Vec3(Half * cc, 0, 0), eps);
        Vec3 cp = Vec3(qcX(yc), yc, zOf(yc)) + qcToCp;

        Vec3 n(0, 0, 1);
        n = RotateAboutY(n, eps);       // apply local twist (about spanwise axis)
        n = RotateAboutX(n, dihedral * (yc >= 0 ? 1.0 : -1.0)); // apply dihedral tilt
        n = n.Normalized();

        Panel p;
        p.A = A; p.B = B;
        p.ControlPoint = cp;
        p.Normal = n;
        p.TrailDirA = Vec3(1, 0, 0);
        p.TrailDirB = Vec3(1, 0, 0);
        p.SpanwiseWidth = yB - yA;
        p.PlanformArea = Half * (cA + cB) * p.SpanwiseWidth;
        // This lattice is flat, so the only thing separating the true
        // surface from its projection is the dihedral tilt.
        p.Area = p.PlanformArea / std::cos(dihedral);
        panels.push_back(p);
    }
    return panels;
}

// A trailing leg: from its root (A or B), through the panel's optional
// curved near-wake path, then the straight far tail along the stated
// direction. Circulation continuity fixes the segment orientations: the
// A-leg flows from downstream infinity INTO A, the B-leg OUT of B --
// `inbound` selects which. The panel's WakeCoreRadius regularizes every
// leg segment (curved wakes thread between panels; see SegmentVelocity).
[[nodiscard]] inline Vec3 TrailingLegVelocity(const Vec3& P, const Vec3& root,
                                              const std::vector<Vec3>& path, const Vec3& farDir,
                                              double trail, double gamma, double coreRadius,
                                              bool inbound) {
    Vec3 v(0, 0, 0);
    Vec3 previous = root;
    for (const Vec3& point : path) {
        v = v + (inbound ? SegmentVelocity(P, point, previous, gamma, coreRadius)
                         : SegmentVelocity(P, previous, point, gamma, coreRadius));
        previous = point;
    }
    const Vec3 far = previous + farDir * trail;
    v = v + (inbound ? SegmentVelocity(P, far, previous, gamma, coreRadius)
                     : SegmentVelocity(P, previous, far, gamma, coreRadius));
    return v;
}

/**
 * Full horseshoe induced velocity (unit gamma) at point P for panel j,
 * using far-downstream points computed on the fly with `trail`. Legs
 * follow the panel's curved near-wake paths when it carries them.
 */
[[nodiscard]] inline Vec3 HorseshoeVelocity(const Vec3& P, const Panel& pj, double gamma, double trail) {
    Vec3 v = TrailingLegVelocity(P, pj.A, pj.TrailPathA, pj.TrailDirA, trail, gamma,
                                 pj.WakeCoreRadius, true); // leg: downstream-inf -> A
    v = v + SegmentVelocity(P, pj.A, pj.B, gamma);         // bound: A -> B
    v = v + TrailingLegVelocity(P, pj.B, pj.TrailPathB, pj.TrailDirB, trail, gamma,
                                pj.WakeCoreRadius, false); // leg: B -> downstream-inf
    return v;
}

/**
 * Same but skipping the bound segment (used for near-field self-induced
 * velocity at a panel's own bound-vortex midpoint, where the bound segment
 * itself would be singular).
 */
[[nodiscard]] inline Vec3 HorseshoeVelocityNoBound(const Vec3& P, const Panel& pj, double gamma, double trail) {
    Vec3 v = TrailingLegVelocity(P, pj.A, pj.TrailPathA, pj.TrailDirA, trail, gamma,
                                 pj.WakeCoreRadius, true);
    v = v + TrailingLegVelocity(P, pj.B, pj.TrailPathB, pj.TrailDirB, trail, gamma,
                                pj.WakeCoreRadius, false);
    return v;
}

// ------------------------------------------------------- Linear solve ----
// LU decomposition with partial pivoting, factorization separated from
// solve. This matters because the influence matrix depends ONLY on panel
// geometry -- NOT on alpha/beta/rates -- so anything that solves the same
// geometry against multiple right-hand sides (ComputeDerivatives() does
// ~11 per call) should factorize once and reuse it, rather than re-running
// full O(N^3) elimination from scratch every time.
//
// The factorization is done by LAPACK's dgetrf (partial-pivoting LU) and
// each solve by dgetrs: condition-aware pivoting, vectorized BLAS-level
// kernels, and a real INFO-based singularity report. Link against any
// LAPACK provider (OpenBLAS, reference LAPACK, ...), e.g.
//   g++ -std=c++23 -O2 -Iinclude -o vlm_demo src/main.cpp -llapack -lblas
//
// The influence matrix is a single contiguous, row-major NxN buffer viewed
// through a 2D std::mdspan (DenseMatrixView) -- one allocation with good
// cache locality, no vector-of-vectors indirection.
//
// LAPACK is column-major; a row-major NxN buffer IS the transpose A^T in
// column-major terms, so we factorize A^T and then ask dgetrs to solve
// op(A^T) x = b with op = 'T', i.e. A x = b -- no explicit transpose needed
// on either side. The default mdspan layout (layout_right) is exactly this
// row-major order.
//
// Nothing above this point in the file (PreparedSystem, SolveWithSystem,
// ComputeDerivatives) needs to know how the solve is implemented -- it only
// uses the LUFactorization / LuFactorize / LuSolve API below.

// Non-owning 2D views over a caller-owned contiguous buffer. Storage is a
// plain std::vector<double>; the mdspan just re-interprets it as a matrix.
using DenseMatrixView      = std::mdspan<double, std::dextents<std::size_t, 2>>;
using ConstDenseMatrixView = std::mdspan<const double, std::dextents<std::size_t, 2>>;

/** LAPACK dgetrf LU factorization of the influence matrix, ready for repeated dgetrs solves. */
struct LUFactorization {
    int N = 0;
    std::vector<double> Flat; ///< Row-major N*N (== column-major A^T), overwritten in place with L/U by dgetrf.
    std::vector<int> Ipiv;
    bool NearSingular = false;
    // Smallest LU pivot divided by the largest entry ANYWHERE in A.
    //
    // Only meaningful when A's entries are commensurable, which is true of a
    // pure lifting-surface system and false of a mixed one. A horseshoe's
    // normal influence per unit circulation has units of 1/length; a source
    // panel's per unit strength is dimensionless. On a wing-body system the
    // vortex columns therefore dwarf the source columns -- measured at 149:1
    // on the reference airframe -- and this ratio collapses by about that
    // factor while the solve stays exact (relative residual 2e-15).
    //
    // So: do not read this as a conditioning number for a coupled system.
    // Check the residual of A x - b instead, which is what TestAirframe
    // asserts. See also SingularAtIndex, which is dgetrf's own verdict and
    // carries no such caveat.
    double MinPivotRatio = 1.0;
    int SingularAtIndex = 0; ///< dgetrf's INFO: >0 means U(info,info) is exactly zero -- 0 = fully nonsingular.
};

/** Factorize A in place via LAPACK dgetrf (partial-pivoting LU). */
[[nodiscard]] inline LUFactorization LuFactorize(ConstDenseMatrixView A) {
    int n = static_cast<int>(A.extent(0));
    LUFactorization f;
    f.N = n;
    f.Flat.resize((size_t)n * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            f.Flat[(size_t)i * n + j] = A[i, j];

    double scale = std::ranges::fold_left(f.Flat, 0.0,
        [](double acc, double v) { return std::max(acc, std::fabs(v)); });
    scale = std::max(scale, DenormFloor);

    f.Ipiv.resize(n);
    int info = 0;
    dgetrf_(&n, &n, f.Flat.data(), &n, f.Ipiv.data(), &info);
    f.SingularAtIndex = info;
    f.NearSingular = (info != 0);

    for (int i = 0; i < n; ++i) {
        double d = std::fabs(f.Flat[(size_t)i * n + i]) / scale;
        f.MinPivotRatio = std::min(f.MinPivotRatio, d);
    }
    return f;
}

/** Solve A x = b against an already-factorized system via LAPACK dgetrs. */
[[nodiscard]] inline std::vector<double> LuSolve(const LUFactorization& f, const std::vector<double>& b) {
    std::vector<double> x = b; // dgetrs solves in place
    if (f.N == 0) return x;
    const char trans = 'T'; // Flat holds A^T column-major, so 'T' solves A x = b
    const int nrhs = 1;
    int info = 0;
    dgetrs_(&trans, &f.N, &nrhs, f.Flat.data(), &f.N, f.Ipiv.data(),
            x.data(), &f.N, &info);
    return x;
}

// ---------------------------------------------------- the panel system ----
// Everything a solve is posed on: lifting surfaces carrying horseshoe
// vortices, and body surfaces carrying sources. Either list may be empty.
//
// There is ONE set of equations, not two treatments joined together. Every
// unknown is a singularity strength, every equation is flow tangency at a
// control point, and the whole thing is
//
//     [ A_vv  A_vs ] [ Gamma ]   [ w_v - Vkin . n_v ]
//     [ A_sv  A_ss ] [ sigma ] = [ w_s - Vkin . n_s ]
//
// where A_xy is the normal velocity induced at an x control point by unit
// strength on a y element, and w is the normal velocity the boundary
// condition prescribes -- zero for any solid surface, nonzero only where a
// panel transpires (the open base carrying propeller efflux, see
// Lattice::SourcePanel).
//
// The blocking is emergent, not written out: unknowns are indexed vortices
// first then sources, so the four blocks are just regions of one matrix
// filled by one loop. A wing-only case is Sources.empty(), which makes the
// system exactly A_vv Gamma = -Vkin . n_v -- the same equations, with the
// source rows and columns absent rather than skipped by a branch. That is
// why there is no separate wing-only path to keep in step.
/** Everything a solve is posed on: lifting surfaces (vortices) and body surfaces (sources). */
struct PanelSystem {
    std::vector<Panel> Vortices;
    std::vector<SourcePanel> Sources;

    [[nodiscard]] int VortexCount() const { return static_cast<int>(Vortices.size()); }
    [[nodiscard]] int SourceCount() const { return static_cast<int>(Sources.size()); }
    [[nodiscard]] int UnknownCount() const { return VortexCount() + SourceCount(); }
    [[nodiscard]] bool IsVortex(int index) const { return index < VortexCount(); }

    [[nodiscard]] const Vec3& ControlPoint(int equation) const {
        return IsVortex(equation) ? Vortices[static_cast<std::size_t>(equation)].ControlPoint
                                  : Sources[static_cast<std::size_t>(equation - VortexCount())].ControlPoint;
    }
    [[nodiscard]] const Vec3& Normal(int equation) const {
        return IsVortex(equation) ? Vortices[static_cast<std::size_t>(equation)].Normal
                                  : Sources[static_cast<std::size_t>(equation - VortexCount())].Normal;
    }
    // Zero for a solid surface; a lifting panel is always solid.
    [[nodiscard]] double PrescribedNormalVelocity(int equation) const {
        return IsVortex(equation)
                   ? 0.0
                   : Sources[static_cast<std::size_t>(equation - VortexCount())].PrescribedNormalVelocity;
    }

    // Velocity at `point` from unit strength on element `index`.
    [[nodiscard]] Vec3 UnitVelocity(int index, const Vec3& point, double trailLength) const {
        if (IsVortex(index))
            return HorseshoeVelocity(point, Vortices[static_cast<std::size_t>(index)], 1.0, trailLength);
        return SourcePanelVelocity(point, Sources[static_cast<std::size_t>(index - VortexCount())]);
    }

    // One entry of the influence matrix.
    [[nodiscard]] double NormalInfluence(int equation, int element, double trailLength) const {
        // A source panel's own control point lies exactly on its sheet,
        // where the induced velocity is indeterminate. The value there is
        // analytic (see SourceInfluence.h), so it is substituted rather
        // than evaluated.
        if (equation == element && !IsVortex(element)) return SourceSelfInfluence;
        return Dot(UnitVelocity(element, ControlPoint(equation), trailLength), Normal(equation));
    }
};

// A geometry-only factorization: the influence matrix depends solely on
// panel positions/normals, never on alpha/beta/rates/externalField, so
// this is what you build ONCE and reuse across many flight conditions
// (see ComputeDerivatives() below for the motivating case).
/** A geometry-only factorization, built once and reused across many flight conditions. */
struct PreparedSystem {
    PanelSystem System;
    double TrailLength = 0.0;
    LUFactorization Factorization;
};

/** Assemble and factorize the influence matrix for system. */
[[nodiscard]] inline PreparedSystem Prepare(const PanelSystem& system, double trailLength) {
    const int N = system.UnknownCount();
    std::vector<double> aStorage(static_cast<std::size_t>(N) * N, 0.0);
    DenseMatrixView Amat(aStorage.data(), N, N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            Amat[i, j] = system.NormalInfluence(i, j, trailLength);

    PreparedSystem sys;
    sys.System = system;
    sys.TrailLength = trailLength;
    sys.Factorization = LuFactorize(Amat);
    return sys;
}

/**
 * Wing-only convenience. This is not a second code path: it states the
 * reduction literally, by posing the same system with no source panels.
 */
[[nodiscard]] inline PreparedSystem Prepare(const std::vector<Panel>& panels, double trailLength) {
    return Prepare(PanelSystem{panels, {}}, trailLength);
}

// --------------------------------------------------------- the flow field ----
// A SOLVED system's velocity field, evaluable anywhere.
//
// Every force the solve reports comes from asking one question -- "what is
// the velocity here?" -- and that question used to be answered by three
// separate inline loops buried inside SolveWithSystem. Nothing downstream
// could reach the field: a surface streamline, a stagnation point, a wake
// survey, a probe rake all had to restate those loops and then drift out of
// step with them. This states them once, and SolveWithSystem itself now
// goes through it, so there is no second definition left to drift.
//
// THREE evaluations rather than one, because two points in the field are
// genuinely special and neither is an approximation of the general case:
//
//   Velocity(P)                the plain field at an arbitrary point.
//
//   BoundMidpointVelocity(i)   at panel i's own bound-vortex midpoint. Its
//                              own bound segment is singular there and is
//                              EXCLUDED: Kutta-Joukowski needs the flow the
//                              segment sits in, not the flow it makes.
//
//   SourceSurfaceVelocity(k)   at source panel k's own control point, which
//                              lies exactly on its own sheet where the
//                              induced velocity is indeterminate. The
//                              analytic +1/2 sheet jump is substituted; the
//                              in-plane part cancels by symmetry.
//
// The field holds a NON-OWNING pointer to the system it was solved on --
// keep the PreparedSystem alive for as long as the field is used.
/** A solved system's velocity field, evaluable at any point. */
struct FlowField {
    const PanelSystem* System = nullptr; ///< Non-owning; must outlive this field.
    double TrailLength = 0.0;
    std::vector<double> gamma; ///< Circulations, aligned with System->Vortices.
    std::vector<double> sigma; ///< Source strengths, aligned with System->Sources.

    Vec3 Vinf{0, 0, 0};     ///< Translational freestream (see FreestreamVelocity).
    Vec3 Omega{0, 0, 0};    ///< Body rates (p, q, r) [rad/s].
    Vec3 RefPoint{0, 0, 0}; ///< Rotation center.
    /** Optional background PERTURBATION, added on top of freestream+rotation. */
    std::function<Vec3(const Vec3&)> External;

    /** Oncoming flow at a point: freestream, body rotation, external field. */
    [[nodiscard]] Vec3 KinematicVelocity(const Vec3& point) const {
        Vec3 v = Vinf - Cross(Omega, point - RefPoint);
        if (External) v = v + External(point);
        return v;
    }

    /** Velocity induced by every singularity in the system at an arbitrary point. */
    [[nodiscard]] Vec3 InducedVelocity(const Vec3& point) const {
        Vec3 v(0, 0, 0);
        if (!System) return v;
        const int nv = System->VortexCount();
        for (int j = 0; j < nv; ++j)
            v = v + HorseshoeVelocity(point, System->Vortices[static_cast<std::size_t>(j)],
                                      gamma[static_cast<std::size_t>(j)], TrailLength);
        const int ns = System->SourceCount();
        for (int k = 0; k < ns; ++k)
            v = v + SourcePanelVelocity(point, System->Sources[static_cast<std::size_t>(k)]) *
                        sigma[static_cast<std::size_t>(k)];
        return v;
    }

    /** Total velocity at an arbitrary point: kinematic plus induced. */
    [[nodiscard]] Vec3 Velocity(const Vec3& point) const {
        return KinematicVelocity(point) + InducedVelocity(point);
    }

    /** Total velocity at vortex panel i's bound-vortex midpoint, its own bound segment excluded. */
    [[nodiscard]] Vec3 BoundMidpointVelocity(int i) const {
        if (!System) return {0, 0, 0};
        const std::size_t self = static_cast<std::size_t>(i);
        const Vec3 mid = (System->Vortices[self].A + System->Vortices[self].B) * Half;
        Vec3 v = KinematicVelocity(mid);
        const int nv = System->VortexCount();
        for (int j = 0; j < nv; ++j) {
            const std::size_t jj = static_cast<std::size_t>(j);
            v = v + ((j == i) ? HorseshoeVelocityNoBound(mid, System->Vortices[jj], gamma[jj], TrailLength)
                              : HorseshoeVelocity(mid, System->Vortices[jj], gamma[jj], TrailLength));
        }
        const int ns = System->SourceCount();
        for (int k = 0; k < ns; ++k)
            v = v + SourcePanelVelocity(mid, System->Sources[static_cast<std::size_t>(k)]) *
                        sigma[static_cast<std::size_t>(k)];
        return v;
    }

    /** Total velocity on source panel k's own face, its self-influence substituted analytically. */
    [[nodiscard]] Vec3 SourceSurfaceVelocity(int k) const {
        if (!System) return {0, 0, 0};
        const std::size_t self = static_cast<std::size_t>(k);
        const SourcePanel& panel = System->Sources[self];
        Vec3 v = KinematicVelocity(panel.ControlPoint);
        const int nv = System->VortexCount();
        for (int j = 0; j < nv; ++j)
            v = v + HorseshoeVelocity(panel.ControlPoint, System->Vortices[static_cast<std::size_t>(j)],
                                      gamma[static_cast<std::size_t>(j)], TrailLength);
        const int ns = System->SourceCount();
        for (int m = 0; m < ns; ++m) {
            const std::size_t mm = static_cast<std::size_t>(m);
            v = v + ((m == k) ? panel.Normal * (sigma[mm] * SourceSelfInfluence)
                              : SourcePanelVelocity(panel.ControlPoint, System->Sources[mm]) * sigma[mm]);
        }
        return v;
    }
};

/**
 * The flow field of a system already solved for `gamma` (and `sigma`, empty
 * for a wing-only system). `sys` must outlive the returned field.
 */
[[nodiscard]] inline FlowField MakeFlowField(const PreparedSystem& sys, const FreestreamConditions& fc,
                                             std::vector<double> gamma, std::vector<double> sigma,
                                             std::function<Vec3(const Vec3&)> externalField = nullptr) {
    FlowField field;
    field.System = &sys.System;
    field.TrailLength = sys.TrailLength;
    field.gamma = std::move(gamma);
    field.sigma = std::move(sigma);
    field.Vinf = FreestreamVelocity(fc);
    field.Omega = Vec3(fc.p, fc.q, fc.r);
    field.RefPoint = fc.RefPoint;
    field.External = std::move(externalField);
    return field;
}

/**
 * Core solver, RHS/force part: takes an already-factorized system and just
 * needs the flight condition -- O(N^2) instead of O(N^3), since the
 * expensive part (matrix assembly + factorization) already happened in
 * Prepare(). This is the function to call in a loop over many conditions
 * against the same geometry.
 */
[[nodiscard]] inline SolveResult SolveWithSystem(const PreparedSystem& sys, const FreestreamConditions& fc,
                                    const ReferenceGeometry& ref,
                                    const std::function<Vec3(const Vec3&)>& externalField = nullptr) {
    const PanelSystem& system = sys.System;
    const std::vector<Panel>& panels = system.Vortices;
    const std::vector<SourcePanel>& bodies = system.Sources;
    int N = system.VortexCount();
    int NB = system.SourceCount();

    // The field starts strengthless: only its KINEMATIC half is needed to
    // pose the right-hand side, and the solved strengths are filled in
    // below. From that point on it is the ONE place the velocity anywhere
    // in this solve is defined (see FlowField).
    FlowField field = MakeFlowField(sys, fc, {}, {}, externalField);
    const Vec3 Vinf = field.Vinf;

    // One right-hand side for both row types: the normal velocity the
    // boundary condition prescribes, less what the oncoming flow already
    // supplies. A solid surface prescribes zero, so a wing-only system's
    // rows are exactly -Vkin.n as before.
    const int unknowns = system.UnknownCount();
    std::vector<double> rhs(static_cast<std::size_t>(unknowns), 0.0);
    for (int i = 0; i < unknowns; ++i)
        rhs[static_cast<std::size_t>(i)] =
            system.PrescribedNormalVelocity(i) -
            Dot(field.KinematicVelocity(system.ControlPoint(i)), system.Normal(i));

    const std::vector<double> strengths = LuSolve(sys.Factorization, rhs);
    // Vortices first, then sources -- the ordering the blocking assumes.
    std::vector<double> gamma(strengths.begin(), strengths.begin() + N);
    std::vector<double> sigma(strengths.begin() + N, strengths.end());
    field.gamma = gamma;
    field.sigma = sigma;

    // --- near-field forces & moments ---
    double rho = fc.rho;
    Vec3 dragDir = Vinf.Normalized();     // wind-axes directions use the TRANSLATIONAL
    Vec3 liftDir = Cross(dragDir, Vec3(0, 1, 0)).Normalized(); // freestream only, not rotation
    Vec3 sideDir = Cross(liftDir, dragDir).Normalized();

    Vec3 totalForce(0, 0, 0);
    Vec3 totalMoment(0, 0, 0);
    SolveResult res;
    res.gamma = gamma;
    res.sigma = sigma;
    std::vector<double> panelLift(N, 0.0);

    double S = 0.0;
    // Reference area is PLANFORM area -- the coefficient convention (see
    // Panel::Area vs Panel::PlanformArea).
    for (auto& p : panels) S += p.PlanformArea;
    if (ref.Area > 0.0) S = ref.Area; // caller override (e.g. wing-only area for a multi-surface aircraft)
    res.ReferenceArea = S;
    res.ReferenceChord = (ref.Chord > 0.0) ? ref.Chord : (S > 0.0 ? S / std::max(GeometryEps, Two) : UnitFallbackLength);
    res.ReferenceSpan = (ref.Span > 0.0) ? ref.Span : UnitFallbackLength;

    for (int i = 0; i < N; ++i) {
        Vec3 mid = (panels[i].A + panels[i].B) * Half;
        // The body's sources are part of the flow the bound vortex sits in,
        // so they enter Kutta-Joukowski here too -- this is the wing-body
        // interference (upwash/blockage) acting on the lifting surface.
        Vec3 Vlocal = field.BoundMidpointVelocity(i);
        Vec3 dl = panels[i].B - panels[i].A;
        Vec3 F = Cross(Vlocal, dl) * (rho * gamma[i]); // Kutta-Joukowski, rho * gamma * (V x dl)
        totalForce = totalForce + F;
        totalMoment = totalMoment + Cross(mid - fc.RefPoint, F);

        panelLift[i] = Dot(F, liftDir);

        const std::string& sname = panels[i].Surface;
        res.LiftBySurface[sname] += Dot(F, liftDir);
        res.DragBySurface[sname] += Dot(F, dragDir);
        res.AreaBySurface[sname] += panels[i].PlanformArea;
    }

    // --- body surface pressure ---
    // A source panel carries no bound circulation, so Kutta-Joukowski says
    // nothing about it; its force is the pressure integral over its face.
    // On a closed body the net force very nearly cancels -- that is
    // d'Alembert, and a large residual here means the panelling is too
    // coarse rather than that the body is producing lift. The MOMENT does
    // not cancel: that is the Munk moment, the dominant effect a fuselage
    // has on an airframe's stability, and it is why the body is modelled
    // at all.
    const double qInf = Half * rho * fc.Vinf * fc.Vinf;
    for (int k = 0; k < NB; ++k) {
        const SourcePanel& body = bodies[static_cast<std::size_t>(k)];
        // A permeable face still shaped the flow through its boundary
        // condition, but there is no wall there for pressure to push on.
        // See Lattice::SourcePanel::Permeable.
        if (body.Permeable) continue;
        // The panel's own contribution at its centroid is purely normal and
        // analytic; FlowField::SourceSurfaceVelocity substitutes it.
        const Vec3 Vlocal = field.SourceSurfaceVelocity(k);

        const double speedRatio = (fc.Vinf > CoeffDenomEps) ? Vlocal.Norm() / fc.Vinf : 0.0;
        const double cp = 1.0 - speedRatio * speedRatio;
        // Pressure acts inward along the surface normal.
        const Vec3 F = body.Normal * (-cp * qInf * body.Area);
        totalForce = totalForce + F;
        totalMoment = totalMoment + Cross(body.ControlPoint - fc.RefPoint, F);

        const std::string& sname = body.Surface;
        res.LiftBySurface[sname] += Dot(F, liftDir);
        res.DragBySurface[sname] += Dot(F, dragDir);
    }

    // --- collapse the lattice into spanwise stations ---
    // A chordwise stack of panels sharing a StripIndex is ONE spanwise
    // station: its circulation is the sum of the stack's horseshoe
    // strengths, its lift the sum of theirs, and the chord that
    // nondimensionalizes cl_local is the whole section chord (the stack's
    // areas summed over the strip width) -- not one panel's slice of it.
    // A negative StripIndex keeps a panel its own station, so a
    // single-row lattice reports exactly as it did before.
    std::map<int, std::vector<int>> strips;
    for (int i = 0; i < N; ++i) {
        int key = (panels[i].StripIndex >= 0) ? panels[i].StripIndex : -(i + 1);
        strips[key].push_back(i);
    }

    double qLocal = Half * rho * fc.Vinf * fc.Vinf;
    res.Stations.reserve(strips.size());
    for (const auto& [key, members] : strips) {
        const Panel& first = panels[members.front()];
        StationResult sr;
        sr.y = Half * (first.A.y + first.B.y);
        sr.Surface = first.Surface;
        sr.PanelIndex = members.front();

        double lift = 0.0;
        double area = 0.0;
        for (int i : members) {
            sr.gamma += gamma[i];
            lift += panelLift[i];
            area += panels[i].PlanformArea;
        }
        double width = first.SpanwiseWidth;
        sr.LiftPerSpan = (width > GeometryEps) ? lift / width : 0.0;
        double sectionChord = (width > GeometryEps) ? area / width : 0.0;
        sr.cl_local = (qLocal > CoeffDenomEps && sectionChord > CoeffDenomEps)
                          ? sr.LiftPerSpan / (qLocal * sectionChord)
                          : 0.0;
        res.Stations.push_back(sr);
    }

    res.L = Dot(totalForce, liftDir);
    res.Di = Dot(totalForce, dragDir);
    res.Y = Dot(totalForce, sideDir);
    res.Mx = totalMoment.x; res.My = totalMoment.y; res.Mz = totalMoment.z;

    double q = Half * rho * fc.Vinf * fc.Vinf;
    res.CL = (q * S > CoeffDenomEps) ? res.L / (q * S) : 0.0;
    res.CDi = (q * S > CoeffDenomEps) ? res.Di / (q * S) : 0.0;
    res.CY = (q * S > CoeffDenomEps) ? res.Y / (q * S) : 0.0;
    res.Croll = (q * S * res.ReferenceSpan > CoeffDenomEps) ? res.Mx / (q * S * res.ReferenceSpan) : 0.0;
    res.Cm = (q * S * res.ReferenceChord > CoeffDenomEps) ? res.My / (q * S * res.ReferenceChord) : 0.0;
    res.Cn = (q * S * res.ReferenceSpan > CoeffDenomEps) ? res.Mz / (q * S * res.ReferenceSpan) : 0.0;

    std::ranges::sort(res.Stations, {}, &StationResult::y);

    return res;
}

/**
 * Core solver: works on an arbitrary panel list (one or many lifting
 * surfaces already merged together -- e.g. wing + horizontal tail +
 * vertical tail panels all in one vector).
 *
 * externalField: optional background velocity PERTURBATION as a function
 * of position, added on top of the freestream+rotation kinematic velocity
 * at every panel. This is how you inject something like a propeller's
 * slipstream (see the external BEMT project's SlipstreamField), ground
 * effect, or a gust field into an otherwise-ordinary lifting-surface
 * solve -- the field only needs to return the extra velocity it
 * contributes, not the total.
 *
 * This is a thin convenience wrapper around Prepare()+SolveWithSystem();
 * if you're going to solve the SAME geometry against multiple flight
 * conditions (a sweep, stability derivatives, an optimization loop), call
 * Prepare() once yourself and use SolveWithSystem() directly instead --
 * see ComputeDerivatives() for exactly that pattern.
 */
[[nodiscard]] inline SolveResult Solve(const std::vector<Panel>& panels, const FreestreamConditions& fc,
                          const ReferenceGeometry& ref, double trailLength,
                          const std::function<Vec3(const Vec3&)>& externalField = nullptr) {
    PreparedSystem sys = Prepare(panels, trailLength);
    return SolveWithSystem(sys, fc, ref, externalField);
}

/**
 * Convenience overload: build a single parametric wing and solve it (this
 * is the original single-surface entry point). Reference chord defaults to
 * the mean geometric chord (S/span); moment reference point comes from
 * fc.RefPoint (default: origin -- set fc.RefPoint to your CG).
 */
[[nodiscard]] inline SolveResult Solve(const WingParams& wp, const FreestreamConditions& fc) {
    std::vector<Panel> panels = BuildWing(wp);
    for (auto& p : panels) p.Surface = "wing";
    double trail = (wp.TrailLength > 0.0) ? wp.TrailLength : 50.0 * wp.Span;
    double S = 0.0;
    for (auto& p : panels) S += p.PlanformArea;
    ReferenceGeometry ref;
    ref.Area = S;
    ref.Span = wp.Span;
    ref.Chord = S / wp.Span;
    return Solve(panels, fc, ref, trail);
}

/**
 * Central-difference stability & control derivatives about base, by
 * factorizing the geometry once and re-solving it against 11 perturbed
 * flight conditions via SolveWithSystem().
 */
[[nodiscard]] inline StabilityDerivatives ComputeDerivatives(const std::vector<Panel>& panels,
                                                const FreestreamConditions& base,
                                                const ReferenceGeometry& ref,
                                                double trailLength,
                                                double dAlphaDeg = DefaultAlphaStepDeg,
                                                double dBetaDeg = DefaultBetaStepDeg,
                                                double rateStepFrac = DefaultRateStepFraction) {
    StabilityDerivatives d;
    // Factorize the influence matrix ONCE -- it depends only on panel
    // geometry, never on alpha/beta/rates -- and reuse it for all 11
    // solves below via SolveWithSystem() instead of re-factorizing (an
    // O(N^3) operation) from scratch every time.
    PreparedSystem sys = Prepare(panels, trailLength);
    SolveResult s0 = SolveWithSystem(sys, base, ref);
    d.CL0 = s0.CL; d.CDi0 = s0.CDi; d.CY0 = s0.CY; d.Cm0 = s0.Cm; d.Croll0 = s0.Croll; d.Cn0 = s0.Cn;

    double dA = DegToRad(dAlphaDeg);
    {
        FreestreamConditions fp = base; fp.alphaDeg = base.alphaDeg + dAlphaDeg;
        FreestreamConditions fm = base; fm.alphaDeg = base.alphaDeg - dAlphaDeg;
        SolveResult rp = SolveWithSystem(sys, fp, ref);
        SolveResult rm = SolveWithSystem(sys, fm, ref);
        d.CL_alpha = (rp.CL - rm.CL) / (Two * dA);
        d.CDi_alpha = (rp.CDi - rm.CDi) / (Two * dA);
        d.Cm_alpha = (rp.Cm - rm.Cm) / (Two * dA);
    }
    double dB = DegToRad(dBetaDeg);
    {
        FreestreamConditions fp = base; fp.betaDeg = base.betaDeg + dBetaDeg;
        FreestreamConditions fm = base; fm.betaDeg = base.betaDeg - dBetaDeg;
        SolveResult rp = SolveWithSystem(sys, fp, ref);
        SolveResult rm = SolveWithSystem(sys, fm, ref);
        d.CY_beta = (rp.CY - rm.CY) / (Two * dB);
        d.Croll_beta = (rp.Croll - rm.Croll) / (Two * dB);
        d.Cn_beta = (rp.Cn - rm.Cn) / (Two * dB);
    }

    double cbar = (ref.Chord > 0.0) ? ref.Chord : UnitFallbackLength;
    double b = (ref.Span > 0.0) ? ref.Span : UnitFallbackLength;

    double dq = rateStepFrac * base.Vinf / cbar;
    {
        FreestreamConditions fp = base; fp.q = base.q + dq;
        FreestreamConditions fm = base; fm.q = base.q - dq;
        SolveResult rp = SolveWithSystem(sys, fp, ref);
        SolveResult rm = SolveWithSystem(sys, fm, ref);
        d.CL_q = (rp.CL - rm.CL) / (Two * dq);
        d.CDi_q = (rp.CDi - rm.CDi) / (Two * dq);
        d.Cm_q = (rp.Cm - rm.Cm) / (Two * dq);
        d.CL_q_nd = d.CL_q * (cbar / (Two * base.Vinf));
        d.Cm_q_nd = d.Cm_q * (cbar / (Two * base.Vinf));
    }

    double dp = rateStepFrac * base.Vinf / b;
    {
        FreestreamConditions fp = base; fp.p = base.p + dp;
        FreestreamConditions fm = base; fm.p = base.p - dp;
        SolveResult rp = SolveWithSystem(sys, fp, ref);
        SolveResult rm = SolveWithSystem(sys, fm, ref);
        d.CY_p = (rp.CY - rm.CY) / (Two * dp);
        d.Croll_p = (rp.Croll - rm.Croll) / (Two * dp);
        d.Cn_p = (rp.Cn - rm.Cn) / (Two * dp);
        d.Croll_p_nd = d.Croll_p * (b / (Two * base.Vinf));
        d.Cn_p_nd = d.Cn_p * (b / (Two * base.Vinf));
    }

    double dr = rateStepFrac * base.Vinf / b;
    {
        FreestreamConditions fp = base; fp.r = base.r + dr;
        FreestreamConditions fm = base; fm.r = base.r - dr;
        SolveResult rp = SolveWithSystem(sys, fp, ref);
        SolveResult rm = SolveWithSystem(sys, fm, ref);
        d.CY_r = (rp.CY - rm.CY) / (Two * dr);
        d.Croll_r = (rp.Croll - rm.Croll) / (Two * dr);
        d.Cn_r = (rp.Cn - rm.Cn) / (Two * dr);
        d.Croll_r_nd = d.Croll_r * (b / (Two * base.Vinf));
        d.Cn_r_nd = d.Cn_r * (b / (Two * base.Vinf));
    }

    return d;
}

} // namespace Aeolion::Solver
