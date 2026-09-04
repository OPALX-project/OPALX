/**
 * @file SpaceChargeSolver.h
 * @brief Owns one configured algorithm and aggregates completed work.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVER_H
#define OPALX_SPACE_CHARGE_SOLVER_H

#include "SpaceCharge/SpaceChargeAlgorithm.h"

#include <cstddef>
#include <memory>

namespace opalx::spacecharge {

    /** @brief Stable tracker-facing entry point for all space-charge algorithms. */
    class SpaceChargeSolver {
    public:
        SpaceChargeSolver(
                std::unique_ptr<SpaceChargeAlgorithm> algorithm,
                std::size_t particleContainerCount);

        SpaceChargeSolver(const SpaceChargeSolver&)            = delete;
        SpaceChargeSolver& operator=(const SpaceChargeSolver&) = delete;
        SpaceChargeSolver(SpaceChargeSolver&&)                 = delete;
        SpaceChargeSolver& operator=(SpaceChargeSolver&&)      = delete;

        /**
         * @brief Validate per-container activity and dispatch one tracker step.
         *
         * Exceptions are terminal for the current run. Once dispatch begins, transient particle,
         * mesh, field, frame, and backend state is unspecified if an operation throws.
         */
        void solve(const SpaceChargeSolveContext& context);

        [[nodiscard]] int reportedBinCount() const { return reportedBinCount_m; }
        [[nodiscard]] std::size_t backendSolveCount() const { return backendSolveCount_m; }
        [[nodiscard]] std::size_t redistributionCount() const { return redistributionCount_m; }

    private:
        std::unique_ptr<SpaceChargeAlgorithm> algorithm_m;
        std::size_t particleContainerCount_m = 0;
        std::size_t backendSolveCount_m      = 0;
        std::size_t redistributionCount_m    = 0;
        int reportedBinCount_m               = 1;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SOLVER_H
