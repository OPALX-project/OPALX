/**
 * @file SpaceChargeSolver.h
 * @brief Owns one configured algorithm and aggregates completed work.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVER_H
#define OPALX_SPACE_CHARGE_SOLVER_H

#include "SpaceCharge/SpaceChargeAlgorithm.h"
#include "Utilities/OpalException.h"

#include <cstddef>
#include <memory>

namespace opalx::spacecharge {

    /** @brief Stable tracker-facing entry point for all space-charge algorithms. */
    class SpaceChargeSolver {
    public:
        SpaceChargeSolver(
                std::unique_ptr<SpaceChargeAlgorithm> algorithm, std::size_t particleContainerCount,
                bool fft2d5 = false)
            : algorithm_m(std::move(algorithm)),
              particleContainerCount_m(particleContainerCount),
              fft2d5_m(fft2d5) {
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

        [[nodiscard]] int reportedBinCount() const { return reportedBinCount_m; }
        [[nodiscard]] std::size_t backendSolveCount() const { return backendSolveCount_m; }
        [[nodiscard]] std::size_t redistributionCount() const { return redistributionCount_m; }
        /// True only when the factory selected FFT2D5; enables its beam-parameter fast path.
        [[nodiscard]] bool usesFFT2D5() const { return fft2d5_m; }

    private:
        std::unique_ptr<SpaceChargeAlgorithm> algorithm_m;
        std::size_t particleContainerCount_m = 0;
        bool fft2d5_m                        = false;
        std::size_t backendSolveCount_m      = 0;
        std::size_t redistributionCount_m    = 0;
        int reportedBinCount_m               = 1;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SOLVER_H
