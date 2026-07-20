#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/VLM/VLM.h"
#include <print>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::println(stderr, "usage: aeolion_geometry <aeolion_geometry.json>");
        return 2;
    }
    try {
        const auto contract = Aeolion::Geometry::LoadHandoff(argv[1]);
        const auto wing = Aeolion::Geometry::ToWingParams(contract);

        // Only wing-bound controls belong to this lattice; duct-jet vanes are
        // in the contract for the jet model, not for the VLM (see
        // ControlSurface.h).
        int wingControls = 0;
        for (const auto& surface : contract.ControlSurfaces)
            if (surface.Binding == Aeolion::Geometry::ControlSurfaceBinding::Wing) ++wingControls;

        Aeolion::VLM::FreestreamConditions conditions;
        const auto result = Aeolion::VLM::Solve(wing, conditions);
        std::print("design_id={}\n"
                   "schema_version={}\n"
                   "stations={}\n"
                   "wing_controls={}\n"
                   "panels={}\n"
                   "reference_area_m2={:.10g}\n"
                   "CL={:.10g}\n"
                   "CDi={:.10g}\n",
                   contract.DesignId, contract.SchemaVersion, contract.Stations.size(), wingControls,
                   result.gamma.size(), result.ReferenceArea, result.CL, result.CDi);
    } catch (const std::exception& error) {
        std::println(stderr, "geometry error: {}", error.what());
        return 1;
    }
}
