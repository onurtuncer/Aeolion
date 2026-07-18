// AirfoilPanel.h
//
// 2D panel method for a single airfoil section: constant-strength source
// panels (capture thickness/camber displacement) plus ONE constant vortex
// strength shared by the whole surface (captures circulation/lift), with
// the Kutta condition setting that vortex strength. This is the classic
// "Hess-Smith" style formulation (see e.g. Moran, "An Introduction to
// Theoretical and Computational Aerodynamics", ch. 3) -- simple, robust
// closed-form influence coefficients, good surface velocity distribution
// for feeding a boundary-layer method downstream.
//
// This is a POTENTIAL-FLOW method: it has no viscosity of its own. Its job
// here is only to produce Ue(s)/Vinf around the section, which
// BoundaryLayer.h then uses to solve the actual viscous boundary layer.
//
// No external dependencies beyond the standard library.

#pragma once
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace Aeolion { namespace Airfoil {

constexpr double Pi = 3.14159265358979323846;

struct Vec2 { double x = 0, y = 0; };
inline Vec2 operator-(const Vec2& a, const Vec2& b) { return {a.x - b.x, a.y - b.y}; }
inline Vec2 operator+(const Vec2& a, const Vec2& b) { return {a.x + b.x, a.y + b.y}; }
inline Vec2 operator*(const Vec2& a, double s) { return {a.x * s, a.y * s}; }
inline double Dot2(const Vec2& a, const Vec2& b) { return a.x * b.x + a.y * b.y; }
inline double Norm2(const Vec2& a) { return std::sqrt(a.x * a.x + a.y * a.y); }

// ------------------------------------------------------- NACA 4-digit ----
// Standard closed-form NACA 4-digit thickness/camber formulas. Returns a
// closed loop of N+1 points (N panels), starting and ending at the
// trailing edge, going around the LOWER surface first then the UPPER
// surface (this ordering matters for the panel method's normal-direction
// convention below). Cosine-spaced for panel clustering near LE/TE.
inline std::vector<Vec2> Naca4(const std::string& code, int nPanelsPerSide = 80, double chord = 1.0) {
    if (code.size() != 4) throw std::invalid_argument("Naca4: code must be 4 digits, e.g. \"2412\"");
    int m = code[0] - '0';        // max camber, in % chord
    int p = code[1] - '0';        // location of max camber, in tenths of chord
    int tt = (code[2] - '0') * 10 + (code[3] - '0'); // thickness, in % chord
    double mC = m / 100.0, pC = p / 10.0, t = tt / 100.0;

    auto thickness = [&](double x) {
        return 5.0 * t * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1015 * x * x * x * x);
    };
    auto camber = [&](double x) -> std::pair<double, double> { // returns (yc, dyc/dx)
        if (pC <= 1e-9 || mC <= 1e-9) return {0.0, 0.0};
        if (x < pC) {
            double yc = mC / (pC * pC) * (2 * pC * x - x * x);
            double dyc = 2 * mC / (pC * pC) * (pC - x);
            return {yc, dyc};
        } else {
            double yc = mC / ((1 - pC) * (1 - pC)) * ((1 - 2 * pC) + 2 * pC * x - x * x);
            double dyc = 2 * mC / ((1 - pC) * (1 - pC)) * (pC - x);
            return {yc, dyc};
        }
    };

    int N = nPanelsPerSide;
    std::vector<Vec2> lower(N + 1), upper(N + 1);
    for (int i = 0; i <= N; ++i) {
        double beta = Pi * i / N;
        double x = 0.5 * (1 - std::cos(beta)); // cosine spacing, 0..1
        double yt = thickness(x);
        auto [yc, dyc] = camber(x);
        double theta = std::atan(dyc);
        upper[i] = {(x - yt * std::sin(theta)) * chord, (yc + yt * std::cos(theta)) * chord};
        lower[i] = {(x + yt * std::sin(theta)) * chord, (yc - yt * std::cos(theta)) * chord};
    }

    std::vector<Vec2> loop;
    loop.reserve(2 * N + 1);
    for (int i = N; i >= 0; --i) loop.push_back(lower[i]);       // TE -> LE along lower
    for (int i = 1; i <= N; ++i) loop.push_back(upper[i]);       // LE -> TE along upper
    return loop;
}

// ------------------------------------------------------------- Panel -----
struct Panel {
    Vec2 P1, P2;      // endpoints (loop order)
    Vec2 Mid;         // control point
    Vec2 Tangent;     // unit, P1->P2
    Vec2 Normal;      // unit, outward (see orientation note below)
    double Length = 0.0;
};

// Builds panels from a closed loop. Orientation convention: for a loop
// traversed lower-surface (TE->LE) then upper-surface (LE->TE) -- i.e.
// CLOCKWISE when the airfoil is drawn with LE left, TE right, upper
// surface up -- the outward normal is (-tangent.y, tangent.x).
inline std::vector<Panel> BuildPanels(const std::vector<Vec2>& loop) {
    if (loop.size() < 4) throw std::invalid_argument("BuildPanels: need at least 3 panels");
    std::vector<Panel> panels;
    panels.reserve(loop.size() - 1);
    for (size_t i = 0; i + 1 < loop.size(); ++i) {
        Panel pnl;
        pnl.P1 = loop[i]; pnl.P2 = loop[i + 1];
        Vec2 d = pnl.P2 - pnl.P1;
        pnl.Length = Norm2(d);
        if (pnl.Length < 1e-12) continue; // skip degenerate (duplicate point)
        pnl.Tangent = d * (1.0 / pnl.Length);
        pnl.Normal = Vec2{-pnl.Tangent.y, pnl.Tangent.x};
        pnl.Mid = (pnl.P1 + pnl.P2) * 0.5;
        panels.push_back(pnl);
    }
    return panels;
}

// Induced velocity (in GLOBAL coords) at point P from panel j's unit-
// strength constant source and unit-strength constant vortex, evaluated by
// transforming into the panel's local frame (panel along local x, 0..L).
// Classic closed-form (see e.g. Katz & Plotkin SOR2DC/VOR2DC, or Moran ch.3).
struct PanelInfluence { Vec2 FromSource, FromVortex; };

inline PanelInfluence ComputePanelInfluence(const Vec2& P, const Panel& pnl) {
    // local coords: x' along tangent, y' along normal-ish... use standard
    // (tangent, normal) basis with normal pointing "up" from the panel.
    Vec2 rel = P - pnl.P1;
    double xl = Dot2(rel, pnl.Tangent);
    double yl = Dot2(rel, pnl.Normal); // local "up" = +normal (normal is outward-pointing)
    double L = pnl.Length;

    // Nudge yl away from exactly zero when evaluating a panel's own control
    // point (or anything numerically on its line): atan2's branch cut at
    // y=0 means floating-point noise of either sign there flips (theta2 -
    // theta1) between +pi and -pi, corrupting the otherwise-exact self-
    // influence result (which must always be +-0.5, a well known closed
    // form). Forcing a consistent tiny positive offset picks the correct
    // branch every time.
    if (std::fabs(yl) < 1e-9 * std::max(L, 1.0)) yl = 1e-9 * std::max(L, 1.0);

    double r1 = std::sqrt(xl * xl + yl * yl);
    double r2 = std::sqrt((xl - L) * (xl - L) + yl * yl);
    r1 = std::max(r1, 1e-12); r2 = std::max(r2, 1e-12);
    double th1 = std::atan2(yl, xl);
    double th2 = std::atan2(yl, xl - L);

    double u_src = (1.0 / (2.0 * Pi)) * std::log(r1 / r2);
    double w_src = (1.0 / (2.0 * Pi)) * (th2 - th1);
    double u_vtx = (1.0 / (2.0 * Pi)) * (th2 - th1);
    double w_vtx = -(1.0 / (2.0 * Pi)) * std::log(r1 / r2);

    // transform (u,w) in (tangent, normal) local basis back to global.
    PanelInfluence out;
    out.FromSource = pnl.Tangent * u_src + pnl.Normal * w_src;
    out.FromVortex = pnl.Tangent * u_vtx + pnl.Normal * w_vtx;
    return out;
}

struct SolveResult {
    std::vector<double> sigma;     // source strength per panel
    double gamma = 0.0;            // single shared vortex strength
    std::vector<double> Ue;        // tangential surface speed / Vinf, per panel (signed: + in +tangent direction)
    std::vector<double> Cp;        // pressure coefficient per panel
    double Cl = 0.0;
};

// Solve for the given freestream angle of attack (radians). Vinf magnitude
// doesn't matter for the potential-flow solve (linear), so results are
// normalized as Ue/Vinf; Cl is dimensionless.
inline SolveResult Solve(const std::vector<Panel>& panels, double alphaRad) {
    int N = static_cast<int>(panels.size());
    Vec2 Vinf{std::cos(alphaRad), std::sin(alphaRad)};

    // N+1 unknowns: sigma_0..sigma_{N-1}, gamma. N tangency equations +
    // 1 Kutta condition (equal speed magnitude at the two panels touching
    // the trailing edge, panel 0 and panel N-1 given our loop ordering).
    std::vector<std::vector<double>> A(N + 1, std::vector<double>(N + 1, 0.0));
    std::vector<double> b(N + 1, 0.0);

    std::vector<std::vector<PanelInfluence>> infl(N, std::vector<PanelInfluence>(N));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            infl[i][j] = ComputePanelInfluence(panels[i].Mid, panels[j]);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = Dot2(infl[i][j].FromSource, panels[i].Normal);
        }
        double gammaCoeff = 0.0;
        for (int j = 0; j < N; ++j) gammaCoeff += Dot2(infl[i][j].FromVortex, panels[i].Normal);
        A[i][N] = gammaCoeff;
        b[i] = -Dot2(Vinf, panels[i].Normal);
    }

    // Kutta condition: tangential speed at panel 0 equals (minus) tangential
    // speed at panel N-1 in a consistent sense -- equivalently, sum of the
    // two panels' tangential-velocity equations set equal. Using surface
    // speed magnitude equality: Vt(panel0) + Vt(panelN-1) = 0 in this
    // panel-normal convention (adjacent TE panels have opposing tangent
    // directions in a closed loop), which is the standard closed-panel
    // Kutta implementation.
    {
        int i0 = 0, iN = N - 1;
        for (int j = 0; j < N; ++j) {
            double c0 = Dot2(infl[i0][j].FromSource, panels[i0].Tangent);
            double cN = Dot2(infl[iN][j].FromSource, panels[iN].Tangent);
            A[N][j] = c0 + cN;
        }
        double g0 = 0.0, gN = 0.0;
        for (int j = 0; j < N; ++j) {
            g0 += Dot2(infl[i0][j].FromVortex, panels[i0].Tangent);
            gN += Dot2(infl[iN][j].FromVortex, panels[iN].Tangent);
        }
        A[N][N] = g0 + gN;
        b[N] = -(Dot2(Vinf, panels[i0].Tangent) + Dot2(Vinf, panels[iN].Tangent));
    }

    // Gaussian elimination with partial pivoting.
    int M = N + 1;
    for (int col = 0; col < M; ++col) {
        int piv = col; double best = std::fabs(A[col][col]);
        for (int r = col + 1; r < M; ++r) if (std::fabs(A[r][col]) > best) { best = std::fabs(A[r][col]); piv = r; }
        if (piv != col) { std::swap(A[piv], A[col]); std::swap(b[piv], b[col]); }
        if (std::fabs(A[col][col]) < 1e-13) continue;
        for (int r = col + 1; r < M; ++r) {
            double f = A[r][col] / A[col][col];
            if (f == 0.0) continue;
            for (int c = col; c < M; ++c) A[r][c] -= f * A[col][c];
            b[r] -= f * b[col];
        }
    }
    std::vector<double> x(M, 0.0);
    for (int r = M - 1; r >= 0; --r) {
        double s = b[r];
        for (int c = r + 1; c < M; ++c) s -= A[r][c] * x[c];
        x[r] = (std::fabs(A[r][r]) > 1e-13) ? s / A[r][r] : 0.0;
    }

    SolveResult res;
    res.sigma.assign(x.begin(), x.begin() + N);
    res.gamma = x[N];
    res.Ue.resize(N);
    res.Cp.resize(N);

    for (int i = 0; i < N; ++i) {
        Vec2 v = Vinf;
        for (int j = 0; j < N; ++j) {
            v = v + infl[i][j].FromSource * res.sigma[j];
            v = v + infl[i][j].FromVortex * res.gamma;
        }
        double vt = Dot2(v, panels[i].Tangent);
        res.Ue[i] = vt;
        res.Cp[i] = 1.0 - vt * vt;
    }

    double totalLength = 0.0;
    for (auto& p : panels) totalLength += p.Length; // chord-ish normalization handled by caller via actual chord
    (void)totalLength;

    // Circulation-based Cl (Kutta-Joukowski): Gamma = gamma * (sum of panel
    // lengths) since vortex strength is constant per unit length over the
    // whole surface; Cl = 2*Gamma / (Vinf_mag * chord) -- caller passes
    // chord=1 convention (Naca4() chord scaling handled upstream), so here
    // we just return Cl assuming unit freestream and let the caller divide
    // by chord if they built the airfoil at chord != 1.
    double totalPanelLength = 0.0;
    for (auto& p : panels) totalPanelLength += p.Length;
    double Gamma = res.gamma * totalPanelLength;
    res.Cl = 2.0 * Gamma; // (Vinf magnitude = 1 by construction; divide by chord externally if chord != 1)

    return res;
}

}} // namespace Aeolion::Airfoil
