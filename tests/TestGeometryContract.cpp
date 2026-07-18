#include "Aeolion/GeometryContract.h"
#include <cmath>
#include <iostream>

int main() {
    auto json = nlohmann::json::parse(R"({
      "schema_version": 1,
      "units": {"length": "m", "angle": "deg"},
      "surfaces": [{
        "name": "wing", "kind": "lifting_surface", "airfoil": "NACA 4412",
        "planform": {"span_m": 1.0, "root_chord_m": 0.2, "tip_chord_m": 0.2,
          "quarter_chord_sweep_deg": 0.0, "dihedral_deg": 0.0, "tip_twist_deg": 0.0},
        "discretization": {"panels_per_semispan": 12, "span_spacing": "cosine"}
      }]
    })");
    const auto contract = Aeolion::Geometry::Parse(json);
    if (std::fabs(contract.Wing.Span - 1.0) > 1e-12 ||
        contract.Wing.NPanelsSemiSpan != 12 || contract.Airfoil != "NACA 4412") {
        std::cerr << "geometry contract mapping failed\n";
        return 1;
    }
    return 0;
}
