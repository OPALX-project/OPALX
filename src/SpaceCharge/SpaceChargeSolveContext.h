/**
 * @file SpaceChargeSolveContext.h
 * @brief Defines all borrowed and per-call state needed for one space-charge solve.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H
#define OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H

#include <cstddef>
#include <cstdint>
#include <span>
#include "Algorithms/CoordinateSystemTrafo.h"

namespace opalx::spacecharge {

    /** @brief Coordinate transforms owned by the current solve context. */
    struct CoordinateFrameTransforms {
        CoordinateSystemTrafo trackerToSolve;
        CoordinateSystemTrafo solveToTracker;
    };

    /** @brief Tracker state captured for one space-charge solve. */
    struct SpaceChargeStepState {
        std::size_t step       = 0;
        double time            = 0.0;
        double timeStep        = 0.0;
        bool emissionActive    = false;
        double emittedFraction = 1.0;
        int mpiSize            = 1;
        CoordinateFrameTransforms frames;
    };

    /**
     * @brief Borrowed container activity and immutable state for one space-charge call.
     *
     * Concrete algorithms borrow stable particle containers at construction. The context therefore
     * carries only per-step tracker state and one activity byte per container.
     */
    class SpaceChargeSolveContext {
    public:
        SpaceChargeSolveContext(
                std::span<const std::uint8_t> trackingActive, SpaceChargeStepState stepState);

        [[nodiscard]] std::span<const std::uint8_t> trackingActive() const {
            return trackingActive_m;
        }
        [[nodiscard]] const SpaceChargeStepState& stepState() const { return stepState_m; }

    private:
        std::span<const std::uint8_t> trackingActive_m;
        SpaceChargeStepState stepState_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H
