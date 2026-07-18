// Vlm.h
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
// The dense linear solve uses LAPACK (dgetrf/dgetrs via the Fortran ABI);
// link against a LAPACK provider such as OpenBLAS or reference LAPACK.

#pragma once
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <map>
#include <functional>

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

namespace Aeolion { namespace Vlm {

constexpr double Pi = 3.14159265358979323846;

// ---------------------------------------------------------------- Vec3 ---
struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator-() const { return {-x, -y, -z}; }

    double Norm() const { return std::sqrt(x * x + y * y + z * z); }

    Vec3 Normalized() const {
        double n = Norm();
        if (n < 1e-14) return {0, 0, 0};
        return {x / n, y / n, z / n};
    }
};

inline double Dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

inline Vec3 Cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline Vec3 RotateAboutX(const Vec3& v, double angleRad) {
    double c = std::cos(angleRad), s = std::sin(angleRad);
    return {v.x, c * v.y - s * v.z, s * v.y + c * v.z};
}

// Rotate about a local spanwise axis (approximated as global y) — used to
// apply geometric twist to the panel normal vector.
inline Vec3 RotateAboutY(const Vec3& v, double angleRad) {
    double c = std::cos(angleRad), s = std::sin(angleRad);
    return {c * v.x + s * v.z, v.y, -s * v.x + c * v.z};
}

// --------------------------------------------- finite vortex segment -----
// Biot-Savart induced velocity at point P from a straight vortex filament
// running P1 -> P2 with circulation strength gamma. A small viscous-core
// cutoff (rc) regularizes the singularity on the filament itself.
inline Vec3 SegmentVelocity(const Vec3& P, const Vec3& P1, const Vec3& P2, double gamma) {
    const double rc = 1e-8;   // core radius, relative-scale-independent cutoff on cross product
    Vec3 r1 = P - P1;
    Vec3 r2 = P - P2;
    Vec3 r0 = P2 - P1;

    Vec3 c = Cross(r1, r2);
    double cmag2 = Dot(c, c);

    double r1n = r1.Norm();
    double r2n = r2.Norm();

    if (cmag2 < rc || r1n < 1e-10 || r2n < 1e-10) {
        return {0, 0, 0};   // point on (or essentially on) the filament -> no contribution
    }

    double K = gamma / (4.0 * Pi * cmag2) * (Dot(r0, r1) / r1n - Dot(r0, r2) / r2n);
    return c * K;
}

// ------------------------------------------------------------- Panel -----
struct Panel {
    Vec3 A, B;              // bound vortex endpoints (quarter-chord), A.y < B.y
    Vec3 ControlPoint;      // three-quarter-chord, mid-span
    Vec3 Normal;            // unit outward normal at control point (includes twist+dihedral)
    Vec3 TrailDirA, TrailDirB; // unit direction of trailing legs from A and from B (downstream)
    double Area = 0.0;      // panel planform area (for Cl-per-station bookkeeping)
    double SpanwiseWidth = 0.0;
    std::string Surface;    // which lifting surface this panel belongs to (e.g. "wing", "htail")
};

// ------------------------------------------------------- Wing geometry ---
struct WingParams {
    double Span = 2.0;            // full span, tip-to-tip [m]
    double RootChord = 0.3;       // [m]
    double TipChord = 0.3;        // [m]
    double SweepQuarterChordDeg = 0.0; // sweep of the quarter-chord line
    double DihedralDeg = 0.0;
    double TwistTipDeg = 0.0;     // linear washout: root=0, tip=TwistTipDeg (negative = washout)
    int NPanelsSemiSpan = 15;     // panels per semi-span (total panels = 2x this)
    bool CosineSpacing = true;    // cluster panels near tips
    double TrailLength = 0.0;     // 0 => auto (50x span)
};

inline double ChordAt(const WingParams& w, double yAbs) {
    double eta = yAbs / (w.Span / 2.0);
    return w.RootChord + (w.TipChord - w.RootChord) * eta;
}

// Builds the panel lattice for a single trapezoidal, linearly-twisted,
// swept, dihedral wing. Trailing legs are aligned with the global +x axis
// (standard small-to-moderate AoA simplification).
inline std::vector<Panel> BuildWing(const WingParams& w) {
    if (w.NPanelsSemiSpan < 1) throw std::invalid_argument("NPanelsSemiSpan must be >= 1");

    int N = 2 * w.NPanelsSemiSpan;
    std::vector<double> yStations(N + 1);
    double halfSpan = w.Span / 2.0;

    for (int i = 0; i <= N; ++i) {
        double u = -1.0 + 2.0 * i / N; // -1..1
        if (w.CosineSpacing) {
            // cosine clustering toward tips: theta in [0,pi], y = -cos(theta)
            double theta = (u + 1.0) * 0.5 * Pi; // 0..pi
            yStations[i] = -halfSpan * std::cos(theta);
        } else {
            yStations[i] = halfSpan * u;
        }
    }

    double sweep = w.SweepQuarterChordDeg * Pi / 180.0;
    double dihedral = w.DihedralDeg * Pi / 180.0;
    double trail = (w.TrailLength > 0.0) ? w.TrailLength : 50.0 * w.Span;

    auto qcX = [&](double y) { return std::fabs(y) * std::tan(sweep); };
    auto zOf = [&](double y) { return std::fabs(y) * std::tan(dihedral); };
    auto twistOf = [&](double y) {
        double eta = std::fabs(y) / halfSpan;
        return (w.TwistTipDeg * Pi / 180.0) * eta;
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

        double yc = 0.5 * (yA + yB);
        double cc = ChordAt(w, std::fabs(yc));
        double eps = twistOf(yc);

        // Twist is a physical rotation of the chordline about the local
        // quarter-chord point -- rotate the quarter->three-quarter-chord
        // offset by eps (about the spanwise axis) rather than just tilting
        // the boundary-condition normal. This keeps the panel geometry
        // itself consistent with what a real (or mesh-derived) twisted
        // wing looks like, not just a small-angle linearization.
        Vec3 qcToCp = RotateAboutY(Vec3(0.5 * cc, 0, 0), eps);
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
        p.Area = 0.5 * (cA + cB) * p.SpanwiseWidth;
        panels.push_back(p);

        (void)trail;
    }
    return panels;
}

// Full horseshoe induced velocity (unit gamma) at point P for panel j,
// using far-downstream points computed on the fly with `trail`.
inline Vec3 HorseshoeVelocity(const Vec3& P, const Panel& pj, double gamma, double trail) {
    Vec3 A_inf = pj.A + pj.TrailDirA * trail;
    Vec3 B_inf = pj.B + pj.TrailDirB * trail;
    Vec3 v(0, 0, 0);
    v = v + SegmentVelocity(P, A_inf, pj.A, gamma); // leg: downstream-inf -> A
    v = v + SegmentVelocity(P, pj.A, pj.B, gamma);  // bound: A -> B
    v = v + SegmentVelocity(P, pj.B, B_inf, gamma); // leg: B -> downstream-inf
    return v;
}

// Same but skipping the bound segment (used for near-field self-induced
// velocity at a panel's own bound-vortex midpoint, where the bound segment
// itself would be singular).
inline Vec3 HorseshoeVelocityNoBound(const Vec3& P, const Panel& pj, double gamma, double trail) {
    Vec3 A_inf = pj.A + pj.TrailDirA * trail;
    Vec3 B_inf = pj.B + pj.TrailDirB * trail;
    Vec3 v(0, 0, 0);
    v = v + SegmentVelocity(P, A_inf, pj.A, gamma);
    v = v + SegmentVelocity(P, pj.B, B_inf, gamma);
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
//   g++ -std=c++17 -O2 -o geometry_vlm GeometryVlm.cpp -llapack -lblas
//
// LAPACK is column-major; our matrices are row-major std::vector<vector>.
// Storing a row-major NxN buffer hands LAPACK the transpose A^T, so we
// factorize A^T and then ask dgetrs to solve op(A^T) x = b with op = 'T',
// i.e. A x = b -- no explicit transpose needed on either side.
//
// Nothing above this point in the file (PreparedSystem, SolveWithSystem,
// ComputeDerivatives) needs to know how the solve is implemented -- it only
// uses the LUFactorization / LuFactorize / LuSolve API below.

struct LUFactorization {
    int N = 0;
    std::vector<double> Flat; // row-major N*N (== column-major A^T), overwritten in place with L/U by dgetrf
    std::vector<int> Ipiv;
    bool NearSingular = false;
    double MinPivotRatio = 1.0; // diagonal-vs-scale heuristic, for a quick eyeballed diagnostic (see also SingularAtIndex)
    int SingularAtIndex = 0; // dgetrf's INFO: >0 means U(info,info) is exactly zero -- 0 = fully nonsingular
};

inline LUFactorization LuFactorize(const std::vector<std::vector<double>>& A) {
    int n = static_cast<int>(A.size());
    LUFactorization f;
    f.N = n;
    f.Flat.resize((size_t)n * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            f.Flat[(size_t)i * n + j] = A[i][j];

    double scale = 0.0;
    for (double v : f.Flat) scale = std::max(scale, std::fabs(v));
    scale = std::max(scale, 1e-300);

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

inline std::vector<double> LuSolve(const LUFactorization& f, const std::vector<double>& b) {
    std::vector<double> x = b; // dgetrs solves in place
    if (f.N == 0) return x;
    const char trans = 'T'; // Flat holds A^T column-major, so 'T' solves A x = b
    const int nrhs = 1;
    int info = 0;
    dgetrs_(&trans, &f.N, &nrhs, f.Flat.data(), &f.N, f.Ipiv.data(),
            x.data(), &f.N, &info);
    return x;
}

// Convenience one-shot wrapper (factorize + solve a single RHS) -- kept for
// call sites that only ever need one solve against a given matrix.
inline std::vector<double> SolveLinear(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
    return LuSolve(LuFactorize(A), b);
}

// --------------------------------------------------------- Solver I/O ----
struct FreestreamConditions {
    double Vinf = 20.0;      // [m/s]
    double alphaDeg = 5.0;   // angle of attack [deg]
    double betaDeg = 0.0;    // sideslip [deg] -- positive beta: relative wind has a +y (right) component
    double rho = 1.225;      // [kg/m^3]
    double p = 0.0;          // body roll  rate about x (aft axis)  [rad/s]
    double q = 0.0;          // body pitch rate about y (right axis) [rad/s]
    double r = 0.0;          // body yaw   rate about z (up axis)    [rad/s]
    Vec3 RefPoint = Vec3(0, 0, 0); // moment reference point AND rotation center (e.g. CG)
};

// Normalization constants for force/moment coefficients. A generic panel
// list carries no planform metadata of its own, so these must be supplied
// by the caller (the WingParams convenience overload fills them in
// automatically for a single parametric wing).
struct ReferenceGeometry {
    double Area = 0.0;   // S, wing reference area [m^2]
    double Chord = 0.0;  // reference chord (mean aerodynamic chord) -- normalizes Cm [m]
    double Span = 0.0;   // reference span -- normalizes Croll, Cn [m]
};

struct StationResult {
    double y = 0.0;
    double gamma = 0.0;
    double LiftPerSpan = 0.0; // dL/dy  [N/m]
    double cl_local = 0.0;    // local (sectional) lift coefficient
    std::string Surface;
    int PanelIndex = -1;      // index into the original (unsorted) panel/gamma array
};

struct SolveResult {
    std::vector<double> gamma;
    std::vector<StationResult> Stations;
    double CL = 0.0;
    double CDi = 0.0;   // INDUCED drag only -- VLM is a potential-flow method and cannot predict
                        // viscous/profile drag. Total CD = CDi + your own CD0 estimate.
    double CY = 0.0;
    double Croll = 0.0; // rolling moment coefficient about RefPoint (positive: right wingtip up)
    double Cm = 0.0;    // pitching moment coefficient about RefPoint (positive: nose up)
    double Cn = 0.0;    // yawing moment coefficient about RefPoint (positive: nose right, about +z)
    double L = 0.0, Di = 0.0, Y = 0.0;
    double Mx = 0.0, My = 0.0, Mz = 0.0; // moments about RefPoint [N*m], same axis convention as Croll/Cm/Cn
    double ReferenceArea = 0.0, ReferenceChord = 0.0, ReferenceSpan = 0.0;
    std::map<std::string, double> LiftBySurface;  // [N], per source surface tag
    std::map<std::string, double> DragBySurface;  // [N]
    std::map<std::string, double> AreaBySurface;  // [m^2]
};

// A geometry-only factorization: the influence matrix depends solely on
// panel positions/normals, never on alpha/beta/rates/externalField, so
// this is what you build ONCE and reuse across many flight conditions
// (see ComputeDerivatives() below for the motivating case).
struct PreparedSystem {
    std::vector<Panel> Panels;
    double TrailLength = 0.0;
    LUFactorization Factorization;
};

inline PreparedSystem Prepare(const std::vector<Panel>& panels, double trailLength) {
    int N = static_cast<int>(panels.size());
    std::vector<std::vector<double>> Amat(N, std::vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) {
        const Vec3& cp = panels[i].ControlPoint;
        for (int j = 0; j < N; ++j) {
            Vec3 v = HorseshoeVelocity(cp, panels[j], 1.0, trailLength);
            Amat[i][j] = Dot(v, panels[i].Normal);
        }
    }
    PreparedSystem sys;
    sys.Panels = panels;
    sys.TrailLength = trailLength;
    sys.Factorization = LuFactorize(Amat);
    return sys;
}

// Core solver, RHS/force part: takes an already-factorized system and just
// needs the flight condition -- O(N^2) instead of O(N^3), since the
// expensive part (matrix assembly + factorization) already happened in
// Prepare(). This is the function to call in a loop over many conditions
// against the same geometry.
inline SolveResult SolveWithSystem(const PreparedSystem& sys, const FreestreamConditions& fc,
                                    const ReferenceGeometry& ref,
                                    const std::function<Vec3(const Vec3&)>& externalField = nullptr) {
    const std::vector<Panel>& panels = sys.Panels;
    int N = static_cast<int>(panels.size());
    double trail = sys.TrailLength;

    double alpha = fc.alphaDeg * Pi / 180.0;
    double beta = fc.betaDeg * Pi / 180.0;
    Vec3 Vinf(fc.Vinf * std::cos(alpha) * std::cos(beta),
              fc.Vinf * std::sin(beta),
              fc.Vinf * std::sin(alpha) * std::cos(beta));

    Vec3 omega(fc.p, fc.q, fc.r);
    auto kinematicVelocity = [&](const Vec3& p) -> Vec3 {
        Vec3 v = Vinf - Cross(omega, p - fc.RefPoint);
        if (externalField) v = v + externalField(p);
        return v;
    };

    std::vector<double> rhs(N, 0.0);
    for (int i = 0; i < N; ++i) rhs[i] = -Dot(kinematicVelocity(panels[i].ControlPoint), panels[i].Normal);

    std::vector<double> gamma = LuSolve(sys.Factorization, rhs);

    // --- near-field forces & moments ---
    double rho = fc.rho;
    Vec3 dragDir = Vinf.Normalized();     // wind-axes directions use the TRANSLATIONAL
    Vec3 liftDir = Cross(dragDir, Vec3(0, 1, 0)).Normalized(); // freestream only, not rotation
    Vec3 sideDir = Cross(liftDir, dragDir).Normalized();

    Vec3 totalForce(0, 0, 0);
    Vec3 totalMoment(0, 0, 0);
    SolveResult res;
    res.gamma = gamma;
    res.Stations.resize(N);

    double S = 0.0;
    for (auto& p : panels) S += p.Area;
    if (ref.Area > 0.0) S = ref.Area; // caller override (e.g. wing-only area for a multi-surface aircraft)
    res.ReferenceArea = S;
    res.ReferenceChord = (ref.Chord > 0.0) ? ref.Chord : (S > 0.0 ? S / std::max(1e-9, 2.0) : 1.0);
    res.ReferenceSpan = (ref.Span > 0.0) ? ref.Span : 1.0;

    for (int i = 0; i < N; ++i) {
        Vec3 mid = (panels[i].A + panels[i].B) * 0.5;
        Vec3 vind(0, 0, 0);
        for (int j = 0; j < N; ++j) {
            if (j == i) {
                vind = vind + HorseshoeVelocityNoBound(mid, panels[j], gamma[j], trail);
            } else {
                vind = vind + HorseshoeVelocity(mid, panels[j], gamma[j], trail);
            }
        }
        Vec3 Vlocal = kinematicVelocity(mid) + vind;
        Vec3 dl = panels[i].B - panels[i].A;
        Vec3 F = Cross(Vlocal, dl) * (rho * gamma[i]); // Kutta-Joukowski, rho * gamma * (V x dl)
        totalForce = totalForce + F;
        totalMoment = totalMoment + Cross(mid - fc.RefPoint, F);

        StationResult sr;
        sr.y = mid.y;
        sr.gamma = gamma[i];
        sr.LiftPerSpan = Dot(F, liftDir) / panels[i].SpanwiseWidth;
        double qLocal = 0.5 * rho * fc.Vinf * fc.Vinf;
        double localChord = panels[i].Area / panels[i].SpanwiseWidth;
        sr.cl_local = (qLocal > 1e-12 && localChord > 1e-12) ? sr.LiftPerSpan / (qLocal * localChord) : 0.0;
        sr.Surface = panels[i].Surface;
        sr.PanelIndex = i;
        res.Stations[i] = sr;

        const std::string& sname = panels[i].Surface;
        res.LiftBySurface[sname] += Dot(F, liftDir);
        res.DragBySurface[sname] += Dot(F, dragDir);
        res.AreaBySurface[sname] += panels[i].Area;
    }

    res.L = Dot(totalForce, liftDir);
    res.Di = Dot(totalForce, dragDir);
    res.Y = Dot(totalForce, sideDir);
    res.Mx = totalMoment.x; res.My = totalMoment.y; res.Mz = totalMoment.z;

    double q = 0.5 * rho * fc.Vinf * fc.Vinf;
    res.CL = (q * S > 1e-12) ? res.L / (q * S) : 0.0;
    res.CDi = (q * S > 1e-12) ? res.Di / (q * S) : 0.0;
    res.CY = (q * S > 1e-12) ? res.Y / (q * S) : 0.0;
    res.Croll = (q * S * res.ReferenceSpan > 1e-12) ? res.Mx / (q * S * res.ReferenceSpan) : 0.0;
    res.Cm = (q * S * res.ReferenceChord > 1e-12) ? res.My / (q * S * res.ReferenceChord) : 0.0;
    res.Cn = (q * S * res.ReferenceSpan > 1e-12) ? res.Mz / (q * S * res.ReferenceSpan) : 0.0;

    std::sort(res.Stations.begin(), res.Stations.end(),
              [](const StationResult& a, const StationResult& b) { return a.y < b.y; });

    return res;
}

// Core solver: works on an arbitrary panel list (one or many lifting
// surfaces already merged together -- e.g. wing + horizontal tail +
// vertical tail panels all in one vector).
//
// externalField: optional background velocity PERTURBATION as a function
// of position, added on top of the freestream+rotation kinematic velocity
// at every panel. This is how you inject something like a propeller's
// slipstream (see Bemt.h), ground effect, or a gust field into an
// otherwise-ordinary lifting-surface solve -- the field only needs to
// return the extra velocity it contributes, not the total.
//
// This is a thin convenience wrapper around Prepare()+SolveWithSystem();
// if you're going to solve the SAME geometry against multiple flight
// conditions (a sweep, stability derivatives, an optimization loop), call
// Prepare() once yourself and use SolveWithSystem() directly instead --
// see ComputeDerivatives() for exactly that pattern.
inline SolveResult Solve(const std::vector<Panel>& panels, const FreestreamConditions& fc,
                          const ReferenceGeometry& ref, double trailLength,
                          const std::function<Vec3(const Vec3&)>& externalField = nullptr) {
    PreparedSystem sys = Prepare(panels, trailLength);
    return SolveWithSystem(sys, fc, ref, externalField);
}

// Convenience overload: build a single parametric wing and solve it (this
// is the original single-surface entry point). Reference chord defaults to
// the mean geometric chord (S/span); moment reference point comes from
// fc.RefPoint (default: origin -- set fc.RefPoint to your CG).
inline SolveResult Solve(const WingParams& wp, const FreestreamConditions& fc) {
    std::vector<Panel> panels = BuildWing(wp);
    for (auto& p : panels) p.Surface = "wing";
    double trail = (wp.TrailLength > 0.0) ? wp.TrailLength : 50.0 * wp.Span;
    double S = 0.0;
    for (auto& p : panels) S += p.Area;
    ReferenceGeometry ref;
    ref.Area = S;
    ref.Span = wp.Span;
    ref.Chord = S / wp.Span;
    return Solve(panels, fc, ref, trail);
}

// ------------------------------------------------- Stability derivatives -
// Central-difference stability & control derivatives about a baseline
// flight condition. Angular derivatives (CL_alpha, CY_beta, ...) are per
// RADIAN. Rate derivatives are provided two ways:
//   - "dimensional": d(coefficient) / d(rate in rad/s)
//   - "_nd"         : the conventional nondimensional stability-derivative
//                      form, e.g. Cm_q_nd = dCm / d(q*cbar/(2*Vinf)), which
//                      is what you'll want for a 6-DOF sim or a DAVE-ML
//                      style derivative table.
// NOTE: CDi derivatives reflect INDUCED drag only (see SolveResult::CDi).
struct StabilityDerivatives {
    double CL0 = 0, CDi0 = 0, CY0 = 0, Cm0 = 0, Croll0 = 0, Cn0 = 0;

    double CL_alpha = 0, CDi_alpha = 0, Cm_alpha = 0;      // per rad
    double CY_beta = 0, Croll_beta = 0, Cn_beta = 0;       // per rad

    double CL_q = 0, CDi_q = 0, Cm_q = 0;                  // per rad/s
    double CY_p = 0, Croll_p = 0, Cn_p = 0;                // per rad/s
    double CY_r = 0, Croll_r = 0, Cn_r = 0;                // per rad/s

    double CL_q_nd = 0, Cm_q_nd = 0;                       // x cbar/(2V)
    double Croll_p_nd = 0, Cn_p_nd = 0;                    // x b/(2V)
    double Croll_r_nd = 0, Cn_r_nd = 0;                    // x b/(2V)
};

inline StabilityDerivatives ComputeDerivatives(const std::vector<Panel>& panels,
                                                const FreestreamConditions& base,
                                                const ReferenceGeometry& ref,
                                                double trailLength,
                                                double dAlphaDeg = 0.5,
                                                double dBetaDeg = 0.5,
                                                double rateStepFrac = 0.002) {
    StabilityDerivatives d;
    // Factorize the influence matrix ONCE -- it depends only on panel
    // geometry, never on alpha/beta/rates -- and reuse it for all 11
    // solves below via SolveWithSystem() instead of re-factorizing (an
    // O(N^3) operation) from scratch every time.
    PreparedSystem sys = Prepare(panels, trailLength);
    SolveResult s0 = SolveWithSystem(sys, base, ref);
    d.CL0 = s0.CL; d.CDi0 = s0.CDi; d.CY0 = s0.CY; d.Cm0 = s0.Cm; d.Croll0 = s0.Croll; d.Cn0 = s0.Cn;

    double dA = dAlphaDeg * Pi / 180.0;
    {
        FreestreamConditions fp = base; fp.alphaDeg = base.alphaDeg + dAlphaDeg;
        FreestreamConditions fm = base; fm.alphaDeg = base.alphaDeg - dAlphaDeg;
        SolveResult rp = SolveWithSystem(sys, fp, ref);
        SolveResult rm = SolveWithSystem(sys, fm, ref);
        d.CL_alpha = (rp.CL - rm.CL) / (2 * dA);
        d.CDi_alpha = (rp.CDi - rm.CDi) / (2 * dA);
        d.Cm_alpha = (rp.Cm - rm.Cm) / (2 * dA);
    }
    double dB = dBetaDeg * Pi / 180.0;
    {
        FreestreamConditions fp = base; fp.betaDeg = base.betaDeg + dBetaDeg;
        FreestreamConditions fm = base; fm.betaDeg = base.betaDeg - dBetaDeg;
        SolveResult rp = SolveWithSystem(sys, fp, ref);
        SolveResult rm = SolveWithSystem(sys, fm, ref);
        d.CY_beta = (rp.CY - rm.CY) / (2 * dB);
        d.Croll_beta = (rp.Croll - rm.Croll) / (2 * dB);
        d.Cn_beta = (rp.Cn - rm.Cn) / (2 * dB);
    }

    double cbar = (ref.Chord > 0.0) ? ref.Chord : 1.0;
    double b = (ref.Span > 0.0) ? ref.Span : 1.0;

    double dq = rateStepFrac * base.Vinf / cbar;
    {
        FreestreamConditions fp = base; fp.q = base.q + dq;
        FreestreamConditions fm = base; fm.q = base.q - dq;
        SolveResult rp = SolveWithSystem(sys, fp, ref);
        SolveResult rm = SolveWithSystem(sys, fm, ref);
        d.CL_q = (rp.CL - rm.CL) / (2 * dq);
        d.CDi_q = (rp.CDi - rm.CDi) / (2 * dq);
        d.Cm_q = (rp.Cm - rm.Cm) / (2 * dq);
        d.CL_q_nd = d.CL_q * (cbar / (2 * base.Vinf));
        d.Cm_q_nd = d.Cm_q * (cbar / (2 * base.Vinf));
    }

    double dp = rateStepFrac * base.Vinf / b;
    {
        FreestreamConditions fp = base; fp.p = base.p + dp;
        FreestreamConditions fm = base; fm.p = base.p - dp;
        SolveResult rp = SolveWithSystem(sys, fp, ref);
        SolveResult rm = SolveWithSystem(sys, fm, ref);
        d.CY_p = (rp.CY - rm.CY) / (2 * dp);
        d.Croll_p = (rp.Croll - rm.Croll) / (2 * dp);
        d.Cn_p = (rp.Cn - rm.Cn) / (2 * dp);
        d.Croll_p_nd = d.Croll_p * (b / (2 * base.Vinf));
        d.Cn_p_nd = d.Cn_p * (b / (2 * base.Vinf));
    }

    double dr = rateStepFrac * base.Vinf / b;
    {
        FreestreamConditions fp = base; fp.r = base.r + dr;
        FreestreamConditions fm = base; fm.r = base.r - dr;
        SolveResult rp = SolveWithSystem(sys, fp, ref);
        SolveResult rm = SolveWithSystem(sys, fm, ref);
        d.CY_r = (rp.CY - rm.CY) / (2 * dr);
        d.Croll_r = (rp.Croll - rm.Croll) / (2 * dr);
        d.Cn_r = (rp.Cn - rm.Cn) / (2 * dr);
        d.Croll_r_nd = d.Croll_r * (b / (2 * base.Vinf));
        d.Cn_r_nd = d.Cn_r * (b / (2 * base.Vinf));
    }

    return d;
}

}} // namespace Aeolion::Vlm
