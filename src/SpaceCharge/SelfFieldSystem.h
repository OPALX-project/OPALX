/**
 * @file SelfFieldSystem.h
 * @brief Owns one configured self-field algorithm and validates solve requests.
 */

#ifndef OPALX_SPACE_CHARGE_SELF_FIELD_SYSTEM_H
#define OPALX_SPACE_CHARGE_SELF_FIELD_SYSTEM_H

#include "SpaceCharge/SelfFieldAlgorithm.h"
#include "SpaceCharge/SelfFieldConfig.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace opalx::spacecharge {

    /** @brief Stable tracker-facing entry point for all self-field algorithms. */
    class SelfFieldSystem {
    public:
        SelfFieldSystem(
                SelfFieldConfig config, std::unique_ptr<SelfFieldAlgorithm> algorithm,
                std::vector<ParticleFieldBinding3d> bindings, std::size_t primaryIndex = 0);

        SelfFieldSystem(const SelfFieldSystem&)            = delete;
        SelfFieldSystem& operator=(const SelfFieldSystem&) = delete;
        SelfFieldSystem(SelfFieldSystem&&)                 = delete;
        SelfFieldSystem& operator=(SelfFieldSystem&&)      = delete;

        void solve(SolveContext& context);

        [[nodiscard]] int reportedBinCount() const;
        [[nodiscard]] const SelfFieldDiagnostics& diagnostics() const { return diagnostics_m; }

    private:
        void validateConfiguration() const;
        void validateBindings(const ParticleSetView& particles) const;
        void validateRequest(const SolveContext& context) const;

        SelfFieldConfig config_m;
        std::unique_ptr<SelfFieldAlgorithm> algorithm_m;
        std::vector<ParticleFieldBinding3d> bindings_m;
        std::size_t primaryIndex_m = 0;
        SolverCapabilities capabilities_m;
        SelfFieldDiagnostics diagnostics_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SELF_FIELD_SYSTEM_H
