// Versioned JSON geometry contract shared with vbat-uav-notebooks.
#pragma once

#include "Aeolion/Types/WingParams.h"
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <nlohmann/json.hpp>

namespace Aeolion::Geometry {

inline constexpr int CurrentSchemaVersion = 1;

struct Contract {
    VLM::WingParams Wing;
    std::string Airfoil;
    std::string SurfaceName = "wing";
};

[[nodiscard]] inline Contract Parse(const nlohmann::json& root) {
    if (root.at("schema_version").get<int>() != CurrentSchemaVersion)
        throw std::invalid_argument("unsupported Aeolion geometry schema_version");
    if (root.at("units").at("length") != "m" || root.at("units").at("angle") != "deg")
        throw std::invalid_argument("Aeolion geometry v1 requires metres and degrees");

    const auto& surfaces = root.at("surfaces");
    auto it = std::find_if(surfaces.begin(), surfaces.end(), [](const auto& surface) {
        return surface.value("name", "") == "wing" &&
               surface.value("kind", "") == "lifting_surface";
    });
    if (it == surfaces.end()) throw std::invalid_argument("geometry has no lifting surface named wing");

    const auto& planform = it->at("planform");
    const auto& mesh = it->at("discretization");
    Contract result;
    result.SurfaceName = it->at("name").template get<std::string>();
    result.Airfoil = it->value("airfoil", "");
    result.Wing.Span = planform.at("span_m").template get<double>();
    result.Wing.RootChord = planform.at("root_chord_m").template get<double>();
    result.Wing.TipChord = planform.at("tip_chord_m").template get<double>();
    result.Wing.SweepQuarterChordDeg = planform.at("quarter_chord_sweep_deg").template get<double>();
    result.Wing.DihedralDeg = planform.at("dihedral_deg").template get<double>();
    result.Wing.TwistTipDeg = planform.at("tip_twist_deg").template get<double>();
    result.Wing.NPanelsSemiSpan = mesh.at("panels_per_semispan").template get<int>();
    result.Wing.CosineSpacing = mesh.at("span_spacing").template get<std::string>() == "cosine";
    if (result.Wing.Span <= 0 || result.Wing.RootChord <= 0 || result.Wing.TipChord <= 0)
        throw std::invalid_argument("wing dimensions must be positive");
    return result;
}

[[nodiscard]] inline Contract Load(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open geometry contract: " + path);
    nlohmann::json root;
    input >> root;
    return Parse(root);
}

} // namespace Aeolion::Geometry
