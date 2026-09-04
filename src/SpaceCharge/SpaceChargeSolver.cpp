/**
 * @file SpaceChargeSolver.cpp
 * @brief Implements common space-charge dispatch and result aggregation.
 */

#include "SpaceCharge/SpaceChargeSolver.h"

#include "Utilities/OpalException.h"

#include <utility>

namespace opalx::spacecharge {

    SpaceChargeSolver::SpaceChargeSolver(
            std::unique_ptr<SpaceChargeAlgorithm> algorithm, std::size_t particleContainerCount)
        : algorithm_m(std::move(algorithm)), particleContainerCount_m(particleContainerCount) {
        if (algorithm_m == nullptr) {
            throw OpalException(
                    "SpaceChargeSolver::SpaceChargeSolver", "The space-charge algorithm is null.");
        }
        if (particleContainerCount_m == 0) {
            throw OpalException(
                    "SpaceChargeSolver::SpaceChargeSolver",
                    "Space charge requires at least one particle container.");
        }
    }

    void SpaceChargeSolver::solve(const SpaceChargeSolveContext& context) {
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

}  // namespace opalx::spacecharge
