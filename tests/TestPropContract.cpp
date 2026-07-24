// TestPropContract.cpp -- the handoff-to-propeller bridge, and the hub
// station that broke the solve when it was first wired up.
//
// The bridge is unit conversion, which is exactly the sort of code that
// looks obviously right and silently isn't: the contract states radii as
// fractions of the disk radius and BEMT wants metres. A factor-of-R error
// here would shift every radius and still produce a plausible-looking
// thrust.
#include "Aeolion/BEMT/BEMT.h"
#include "Aeolion/Geometry/HandoffContract.h"

#include <cmath>
#include <iostream>
#include <string>

using Aeolion::Geometry::ContractError;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

Aeolion::Geometry::HandoffContract Fixture() {
    return Aeolion::Geometry::LoadHandoff(std::string(AEOLION_TEST_DATA_DIR) +
                                          "/AeolionGeometryHandoff-1.5.0.json");
}

void TestBridgeConvertsUnits() {
    const auto contract = Fixture();
    const auto& spec = contract.Propulsion;
    const auto prop = Aeolion::Geometry::ToPropGeometry(spec);

    CHECK(prop.NBlades == 3, "blade count should come from the contract");
    CHECK(prop.Stations.size() == spec.BladeStations.size(), "every blade station should survive the bridge");
    CHECK(std::fabs(prop.Radius - spec.DiskRadius) < 1e-15, "tip radius is the disk radius");

    // The contract defines the blade only outboard of its first station, so
    // that station IS the hub cutout.
    CHECK(std::fabs(prop.HubRadius - spec.BladeStations.front().RadiusFraction * spec.DiskRadius) < 1e-15,
          "hub radius is the innermost station, in metres");

    // Fractions became metres, and chord/twist passed through untouched.
    for (std::size_t i = 0; i < prop.Stations.size(); ++i) {
        const double expected = spec.BladeStations[i].RadiusFraction * spec.DiskRadius;
        CHECK(std::fabs(prop.Stations[i].r - expected) < 1e-15,
              "station " + std::to_string(i) + " radius should be r/R times R");
        CHECK(prop.Stations[i].Chord == spec.BladeStations[i].Chord, "chord is already metric");
        CHECK(prop.Stations[i].TwistDeg == spec.BladeStations[i].TwistDeg, "twist passes through");
        if (i > 0)
            CHECK(prop.Stations[i].r > prop.Stations[i - 1].r, "radii must increase outboard");
    }
    CHECK(prop.Stations.back().r <= prop.Radius + 1e-15, "no station may sit outboard of the tip");
}

void TestBridgeRefusesWhatItCannotKnow() {
    auto contract = Fixture();

    // Rotation sense is an installation choice, not blade geometry, and a
    // wrong one flips the slipstream swirl a control vane responds to.
    bool rejected = false;
    try {
        (void)Aeolion::Geometry::ToPropGeometry(contract.Propulsion, 0);
    } catch (const ContractError&) { rejected = true; }
    CHECK(rejected, "a rotation sign of 0 should be refused");
    CHECK(Aeolion::Geometry::ToPropGeometry(contract.Propulsion, -1).RotationSign == -1,
          "left-handed rotation should be accepted");

    // Blade count arrived at schema 1.4.0. An older contract states none,
    // and guessing would scale thrust roughly linearly while looking fine.
    contract.Propulsion.BladeCount = 0;
    rejected = false;
    try {
        (void)Aeolion::Geometry::ToPropGeometry(contract.Propulsion);
    } catch (const ContractError&) { rejected = true; }
    CHECK(rejected, "a contract with no blade count should be refused, not defaulted");
    CHECK(Aeolion::Geometry::ToPropGeometry(contract.Propulsion, 1, 3).NBlades == 3,
          "the caller should be able to supply the missing blade count");
}

// --- the hub station -------------------------------------------------------
// A contract states its blade from the hub cutout outward, so the first
// station lands exactly ON the hub, where Prandtl's loss factor vanishes
// and the momentum balance divides by it. Before this was handled, that one
// station never converged -- at every airspeed, with residuals of hundreds
// of m/s -- while the other twelve were fine.
void TestEveryStationConverges() {
    const auto contract = Fixture();
    const auto prop = Aeolion::Geometry::ToPropGeometry(contract.Propulsion);
    const Aeolion::BEMT::Polar polar;

    for (double speed : {0.0, 5.0, 15.0, 30.0}) {
        const auto result =
            Aeolion::BEMT::Solve(prop, polar, contract.Propulsion.ReferenceRpm, speed);
        CHECK(result.Converged, "every station should converge at V=" + std::to_string(speed));

        for (const auto& station : result.Stations)
            CHECK(station.Converged, "station at r/R " +
                                         std::to_string(station.r / prop.Radius) +
                                         " failed to converge at V=" + std::to_string(speed));

        // The hub and tip stations are unloaded boundaries, which is what
        // the tip/hub relief says, not an artifact of skipping them.
        CHECK(result.Stations.front().dT_dr == 0.0, "the hub station should carry no load");
        CHECK(result.Stations.back().dT_dr == 0.0, "the tip station should carry no load");
    }
}

// Physics that must hold regardless of the blade: thrust falls with
// airspeed at fixed rpm, and both efficiency measures respect their bounds.
void TestPerformanceIsPhysical() {
    const auto contract = Fixture();
    const auto prop = Aeolion::Geometry::ToPropGeometry(contract.Propulsion);
    const Aeolion::BEMT::Polar polar;
    const double rpm = contract.Propulsion.ReferenceRpm;

    const auto hover = Aeolion::BEMT::Solve(prop, polar, rpm, 0.0);
    const double figureOfMerit = Aeolion::BEMT::FigureOfMerit(hover);
    std::cout << "prop: hover T=" << hover.Thrust << " N, P=" << hover.Power << " W, FOM=" << figureOfMerit
              << "\n";
    CHECK(hover.Thrust > 0.0, "the propeller should make thrust in hover");
    CHECK(figureOfMerit > 0.0 && figureOfMerit < 1.0,
          "figure of merit MUST be below 1 -- a thermodynamic bound, not a target");
    CHECK(Aeolion::BEMT::PropulsiveEfficiency(hover, 0.0) == 0.0,
          "propulsive efficiency is zero in hover by definition");

    double previousThrust = hover.Thrust;
    double bestEfficiency = 0.0;
    for (double speed : {5.0, 10.0, 15.0, 20.0, 25.0}) {
        const auto result = Aeolion::BEMT::Solve(prop, polar, rpm, speed);
        const double efficiency = Aeolion::BEMT::PropulsiveEfficiency(result, speed);
        CHECK(result.Thrust < previousThrust,
              "thrust must fall as the propeller unloads with airspeed, at V=" + std::to_string(speed));
        CHECK(efficiency > 0.0 && efficiency < 1.0,
              "propulsive efficiency MUST be below 1 at V=" + std::to_string(speed));
        previousThrust = result.Thrust;
        bestEfficiency = std::max(bestEfficiency, efficiency);
    }
    // A propeller with a sane twist schedule should reach a useful peak
    // efficiency somewhere in its speed range.
    CHECK(bestEfficiency > 0.5, "peak propulsive efficiency implausibly low, got " +
                                    std::to_string(bestEfficiency));
}

} // namespace

int main() {
    TestBridgeConvertsUnits();
    TestBridgeRefusesWhatItCannotKnow();
    TestEveryStationConverges();
    TestPerformanceIsPhysical();

    if (failures == 0) { std::cout << "PASS: TestPropContract\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestPropContract\n";
    return 1;
}
