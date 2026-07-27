#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/Logger/Log.h"
#include "Aeolion/Solver/Solver.h"

int main(int argc, char** argv) {
    Aeolion::Logger::Log::Init();

    if (argc != 2) {
        AE_ERROR("usage: aeolion_geometry <aeolion_geometry.json>");
        return 2;
    }
    try {
        const auto contract = Aeolion::Geometry::LoadHandoff(argv[1]);
        const auto wing = Aeolion::Geometry::ToWingParams(contract);

        // Only wing-bound controls belong to this lattice; duct-jet vanes are
        // in the contract for the jet model, not for the solver (see
        // ControlSurface.h).
        int wingControls = 0;
        for (const auto& surface : contract.ControlSurfaces)
            if (surface.Binding == Aeolion::Geometry::ControlSurfaceBinding::Wing) ++wingControls;

        Aeolion::Solver::FreestreamConditions conditions;
        const auto result = Aeolion::Solver::Solve(wing, conditions);
        AE_INFO("design_id={}", contract.DesignId);
        AE_INFO("schema_version={}", contract.SchemaVersion);
        AE_INFO("stations={}", contract.Stations.size());
        AE_INFO("wing_controls={}", wingControls);
        AE_INFO("panels={}", result.gamma.size());
        AE_INFO("reference_area_m2={:.10g}", result.ReferenceArea);
        AE_INFO("CL={:.10g}", result.CL);
        AE_INFO("CDi={:.10g}", result.CDi);
    } catch (const std::exception& error) {
        AE_ERROR("geometry error: {}", error.what());
        return 1;
    }
}
