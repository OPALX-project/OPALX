/**
 * @file SpaceChargeSolver.h
 * @brief Owns one configured algorithm and aggregates completed work.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVER_H
#define OPALX_SPACE_CHARGE_SOLVER_H

#include "PartBunch/PartBunchFwd.h"
#include "SpaceCharge/SpaceChargeAlgorithm.h"
#include "Utilities/OpalException.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

namespace opalx::spacecharge {

    /** @brief Stable tracker-facing entry point for all space-charge algorithms. */
    class SpaceChargeSolver {
    public:
        SpaceChargeSolver(
                std::unique_ptr<SpaceChargeAlgorithm> algorithm, std::size_t particleContainerCount)
            : algorithm_m(std::move(algorithm)),
              particleContainerCount_m(particleContainerCount),
              containerActivity_m(particleContainerCount, 0) {
            if (algorithm_m == nullptr) {
                throw OpalException(
                        "SpaceChargeSolver::SpaceChargeSolver",
                        "The space-charge algorithm is null.");
            }
            if (particleContainerCount_m == 0) {
                throw OpalException(
                        "SpaceChargeSolver::SpaceChargeSolver",
                        "Space charge requires at least one particle container.");
            }
        }

        SpaceChargeSolver(const SpaceChargeSolver&)            = delete;
        SpaceChargeSolver& operator=(const SpaceChargeSolver&) = delete;
        SpaceChargeSolver(SpaceChargeSolver&&)                 = delete;
        SpaceChargeSolver& operator=(SpaceChargeSolver&&)      = delete;

        /** @brief Validate per-container activity and dispatch one tracker step. */
        void solve(const SpaceChargeSolveContext& context) {
            if (context.trackingActive().size() != particleContainerCount_m) {
                throw OpalException(
                        "SpaceChargeSolver::solve",
                        "The container activity set does not match solver construction.");
            }

            const SpaceChargeSolveResult result = algorithm_m->solve(context);
            backendSolveCount_m += result.backendSolves;
            redistributionCount_m += result.redistributions;
            reportedBinCount_m = result.reportedBins;
        }

        /**
         * @brief The space-charge solve of one tracker step for @p bunch, the bunch whose
         *        containers the algorithm borrows. Collective.
         *
         * Skips the step while the bunch holds at most @p minBinEmitted particles. Otherwise it
         * builds the solve frame the algorithm requests from the primary container's moments
         * (makeSolveFrame()), the container activity and the step state, and calls solve().
         * The particles are in their reference frame on entry and again on return.
         * @param emissionProgress Called once per solve, after the solve frame is built; it may
         *        reduce over the ranks.
         */
        void solveStep(
                const PartBunch_t& bunch, int minBinEmitted,
                const std::function<EmissionProgress()>& emissionProgress);

        [[nodiscard]] SolveFrame solveFrame() const { return algorithm_m->solveFrame(); }

        [[nodiscard]] int reportedBinCount() const { return reportedBinCount_m; }
        [[nodiscard]] std::size_t backendSolveCount() const { return backendSolveCount_m; }
        [[nodiscard]] std::size_t redistributionCount() const { return redistributionCount_m; }

    private:
        /// The solve frame @p request from the primary container's moments. Collective unless
        /// @p request is the reference frame.
        [[nodiscard]] ParticleFrame stepSolveFrame(
                const PartBunch_t& bunch, SolveFrame request) const;

        std::unique_ptr<SpaceChargeAlgorithm> algorithm_m;
        std::size_t particleContainerCount_m = 0;
        /// Per container: 1 if it takes part in the solve of the current step.
        std::vector<std::uint8_t> containerActivity_m;
        std::size_t backendSolveCount_m   = 0;
        std::size_t redistributionCount_m = 0;
        int reportedBinCount_m            = 1;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SOLVER_H
