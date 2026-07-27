// PropellerCLI.cpp -- run the handoff's propeller through BEMT.
//
//   aeolion_prop <aeolion_geometry.json> [--rpm N] [--speed V] [--rho R]
//                                        [--blades N] [--reverse]
//
// Prints hover and forward-flight performance and the radial distribution
// the blade actually works at, so a design can be sanity-checked without
// writing code. With no --speed it sweeps a range, which is the more useful
// default: a single operating point says little about a propeller.

#include "Aeolion/BEMT/BEMT.h"
#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/Logger/Log.h"

#include <charconv>
#include <exception>
#include <format>
#include <string>
#include <string_view>
#include <vector>

namespace {

// Speeds swept when the caller names none. Spread rather than fine: the
// point is to show the trend from hover to cruise, not to resolve it.
constexpr double SweepSpeeds[] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0};

struct Options {
    std::string Path;
    double Rpm = 0.0;      // 0 => the contract's reference rpm
    double Speed = -1.0;   // <0 => sweep
    double Density = Aeolion::BEMT::SeaLevelDensity;
    int Blades = 0;        // 0 => whatever the contract states
    int RotationSign = 1;
};

[[nodiscard]] bool ParseDouble(std::string_view text, double& out) {
    // from_chars for double is not universally available on every stdlib
    // build this targets, so go through the string overload.
    try {
        std::size_t consumed = 0;
        out = std::stod(std::string(text), &consumed);
        return consumed == text.size();
    } catch (const std::exception&) {
        return false;
    }
}

[[nodiscard]] bool ParseInt(std::string_view text, int& out) {
    const auto result = std::from_chars(text.data(), text.data() + text.size(), out);
    return result.ec == std::errc{} && result.ptr == text.data() + text.size();
}

[[nodiscard]] bool ParseArguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        const auto next = [&](double& target) {
            return i + 1 < argc && ParseDouble(argv[++i], target);
        };
        if (arg == "--reverse") { options.RotationSign = -1; continue; }
        if (arg == "--rpm") { if (!next(options.Rpm)) return false; continue; }
        if (arg == "--speed") { if (!next(options.Speed)) return false; continue; }
        if (arg == "--rho") { if (!next(options.Density)) return false; continue; }
        if (arg == "--blades") {
            if (i + 1 >= argc || !ParseInt(argv[++i], options.Blades)) return false;
            continue;
        }
        if (arg.starts_with("--")) return false;
        if (!options.Path.empty()) return false;
        options.Path = arg;
    }
    return !options.Path.empty();
}

void ReportOperatingPoint(const Aeolion::BEMT::PropGeometry& prop, const Aeolion::BEMT::Polar& polar,
                          double rpm, double speed, double density) {
    const auto result = Aeolion::BEMT::Solve(prop, polar, rpm, speed, density);
    const double efficiency = Aeolion::BEMT::PropulsiveEfficiency(result, speed);
    const double figureOfMerit = Aeolion::BEMT::FigureOfMerit(result, density);

    // Figure of merit is a hover metric and propulsive efficiency a
    // forward-flight one; neither means anything in the other regime, so
    // only the applicable one is shown.
    std::string row = std::format("{:>7.1f} {:>10.3f} {:>10.4f} {:>10.1f} {:>9.1f}", speed, result.Thrust,
                                  result.Torque, result.Power, Aeolion::BEMT::DiskLoading(result));
    if (speed > 0.0) row += std::format(" {:>10.3f}", efficiency);
    else             row += std::format(" {:>10}", "-");
    if (speed > 0.0) row += std::format(" {:>8}", "-");
    else             row += std::format(" {:>8.3f}", figureOfMerit);
    if (result.Converged) {
        AE_INFO("{}", row);
    } else {
        int failed = 0;
        double worstRadius = 0.0, worstResidual = 0.0;
        for (const auto& station : result.Stations)
            if (!station.Converged) {
                ++failed;
                if (station.Residual > worstResidual) {
                    worstResidual = station.Residual;
                    worstRadius = station.r;
                }
            }
        AE_WARN("{}   {} of {} stations unconverged (worst r/R {:.2f}, residual {:.2e} m/s)", row, failed,
                result.Stations.size(), worstRadius / prop.Radius, worstResidual);
    }
}

} // namespace

int main(int argc, char** argv) {
    Aeolion::Logger::Log::Init();

    Options options;
    if (!ParseArguments(argc, argv, options)) {
        AE_ERROR("usage: aeolion_prop <aeolion_geometry.json> [--rpm N] [--speed V] [--rho R]\n"
                 "                    [--blades N] [--reverse]\n"
                 "\n"
                 "  --rpm     shaft speed; defaults to the contract's reference_rpm\n"
                 "  --speed   a single airspeed [m/s]; omitted, a range is swept\n"
                 "  --rho     air density [kg/m^3]; defaults to sea-level ISA\n"
                 "  --blades  blade count, for contracts older than schema 1.4.0\n"
                 "  --reverse left-handed rotation (flips slipstream swirl)");
        return 2;
    }

    try {
        const auto contract = Aeolion::Geometry::LoadHandoff(options.Path);
        if (contract.Propulsion.BladeStations.empty()) {
            AE_ERROR("{}: this contract carries no propulsion_bemt block", options.Path);
            return 1;
        }

        const auto prop = Aeolion::Geometry::ToPropGeometry(contract.Propulsion, options.RotationSign,
                                                           options.Blades);
        const double rpm = (options.Rpm > 0.0) ? options.Rpm : contract.Propulsion.ReferenceRpm;

        AE_INFO("design_id      {}", contract.DesignId);
        AE_INFO("blades         {}{}", prop.NBlades,
                contract.Propulsion.BladeCount > 0 ? " (from contract)" : " (supplied)");
        AE_INFO("radius         {:.4f} m   hub {:.4f} m   stations {}", prop.Radius, prop.HubRadius,
                prop.Stations.size());
        AE_INFO("rpm            {:.1f}{}", rpm,
                options.Rpm > 0.0 ? "" : "   (contract reference)");
        AE_INFO("density        {:.4f} kg/m^3", options.Density);

        // The analytic polar is a placeholder, and saying so matters: blade
        // section data is what separates a believable propeller prediction
        // from a plausible-looking one.
        const Aeolion::BEMT::Polar polar;
        AE_WARN("using the analytic default polar (thin-airfoil lift slope, parabolic drag);\n"
                "thrust is far less sensitive to this than torque and efficiency are");

        AE_INFO("{:>7} {:>10} {:>10} {:>10} {:>9} {:>10} {:>8}", "V[m/s]", "T[N]", "Q[N.m]", "P[W]",
                "T/A[Pa]", "eta", "FOM");
        if (options.Speed >= 0.0) {
            ReportOperatingPoint(prop, polar, rpm, options.Speed, options.Density);
        } else {
            for (const double speed : SweepSpeeds)
                ReportOperatingPoint(prop, polar, rpm, speed, options.Density);
        }
    } catch (const std::exception& error) {
        AE_ERROR("propeller error: {}", error.what());
        return 1;
    }
    return 0;
}
