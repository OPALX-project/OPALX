/**
 * @file SpaceChargeSolver.h
 * @brief Owns one configured space-charge algorithm and validates solve requests.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVER_H
#define OPALX_SPACE_CHARGE_SOLVER_H

#include "SpaceCharge/SpaceChargeAlgorithm.h"
#include "SpaceCharge/SpaceChargeConfig.h"

#include <cstddef>
#include <memory>
#include <vector>

class BunchStateHandler;

namespace opalx::spacecharge {

    /** @brief Stable tracker-facing entry point for all space-charge algorithms. */
    class SpaceChargeSolver {
    public:
        SpaceChargeSolver(
                SpaceChargeConfig config, std::unique_ptr<SpaceChargeAlgorithm> algorithm,
                std::vector<ParticleFieldBinding3D> bindings,
                std::shared_ptr<const BunchStateHandler> bunchState, std::size_t primaryIndex = 0);

        SpaceChargeSolver(const SpaceChargeSolver&)            = delete;
        SpaceChargeSolver& operator=(const SpaceChargeSolver&) = delete;
        SpaceChargeSolver(SpaceChargeSolver&&)                 = delete;
        SpaceChargeSolver& operator=(SpaceChargeSolver&&)      = delete;

        /**
         * @brief Validate and dispatch one tracker-owned request.
         *
         * Exceptions are terminal for the current run. Once dispatch begins, transient particle,
         * mesh, field, frame, and backend state is unspecified if an operation throws.
         */
        void solve(SpaceChargeSolveContext& context);

        [[nodiscard]] int reportedBinCount() const;
        [[nodiscard]] const SpaceChargeDiagnostics& diagnostics() const { return diagnostics_m; }

    private:
        void validateConfiguration() const;
        void validateBindings(const ParticleFieldSet& particles) const;
        void validateRequest(const SpaceChargeSolveContext& context) const;

        SpaceChargeConfig config_m;
        std::unique_ptr<SpaceChargeAlgorithm> algorithm_m;
        std::vector<ParticleFieldBinding3D> bindings_m;
        std::size_t primaryIndex_m = 0;
        std::shared_ptr<const BunchStateHandler> bunchState_m;
        SpaceChargeCapabilities capabilities_m;
        SpaceChargeDiagnostics diagnostics_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SOLVER_H
