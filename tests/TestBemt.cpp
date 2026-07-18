// TestBemt.cpp -- validates Bemt.h against HARD physical bounds, not
// just "looks plausible" checks:
//   - Figure of Merit (ideal hover power / actual power) must be <= 1.0.
//     This is a real thermodynamic constraint (actual power can never be
//     less than the ideal-uniform-loading minimum for the same thrust),
//     so it's a strong bug detector -- this caught two real sign/exponent
//     bugs during development.
//   - propulsive efficiency (thrust*V/power) must be <= 1.0 for all
//     forward-flight points where the prop is actually absorbing power.
//   - every station must converge.
#include "Aeolion/Bemt.h"
#include <numbers>
#include <iostream>
#include <cmath>

using namespace Aeolion;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

static Bemt::PropGeometry MakeTestProp() {
    Bemt::PropGeometry geom;
    geom.NBlades = 2; geom.Radius = 0.15; geom.HubRadius = 0.02;
    int N = 12;
    for (int i = 0; i < N; ++i) {
        // interior stations only -- never exactly at hub/tip radius, where
        // the Prandtl loss factor is deliberately singular (see Bemt.h)
        double r = geom.HubRadius + (geom.Radius - geom.HubRadius) * (i + 0.5) / (N + 1);
        double eta = r / geom.Radius;
        geom.Stations.push_back({r, 0.025 * (1 - 0.3 * eta), 30.0 - 20.0 * eta});
    }
    return geom;
}

int main() {
    auto geom = MakeTestProp();
    Bemt::Polar polar;
    double rho = 1.225, rpm = 6000.0;

    auto hover = Bemt::Solve(geom, polar, rpm, 0.0, rho);
    std::cout << "hover: converged=" << hover.Converged << " thrust=" << hover.Thrust
              << " torque=" << hover.Torque << " power=" << hover.Power << "\n";
    CHECK(hover.Converged, "hover case should converge");
    CHECK(hover.Thrust > 0, "hover thrust should be positive for this blade pitch schedule");
    CHECK(hover.Power > 0, "hover power should be positive (prop absorbing shaft power)");

    double A = std::numbers::pi * geom.Radius * geom.Radius;
    double idealPower = hover.Thrust * std::sqrt(hover.Thrust / (2 * rho * A));
    double FOM = idealPower / hover.Power;
    std::cout << "Figure of Merit = " << FOM << "\n";
    CHECK(FOM > 0.0 && FOM < 1.0, "Figure of Merit MUST be <= 1.0 -- this is a hard physical bound, not a tuning target");
    CHECK(FOM > 0.3, "Figure of Merit implausibly low for a reasonable blade design -- suspect a sign/scale bug");

    for (double V : {5.0, 10.0, 15.0}) {
        auto r = Bemt::Solve(geom, polar, rpm, V, rho);
        CHECK(r.Converged, "forward-flight case should converge, V=" << V);
        if (r.Power > 1e-6 && r.Thrust > 0) {
            double eff = r.Thrust * V / r.Power;
            std::cout << "V=" << V << "  T=" << r.Thrust << "  P=" << r.Power << "  eta=" << eff << "\n";
            CHECK(eff > 0.0 && eff < 1.0, "propulsive efficiency MUST be <= 1.0 -- hard physical bound, V=" << V);
        }
    }

    if (failures == 0) { std::cout << "PASS: TestBemt\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestBemt\n";
    return 1;
}
