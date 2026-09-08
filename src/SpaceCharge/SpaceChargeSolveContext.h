/**
 * @file SpaceChargeSolveContext.h
 * @brief Defines all borrowed and per-call state needed for one space-charge solve.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H
#define OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H

#include "Algorithms/CoordinateSystemTrafo.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace opalx::spacecharge {

    /** @brief Inverse spatial transforms between tracker and solve axes; no Lorentz boost. */
    struct CoordinateFrameTransforms {
        CoordinateSystemTrafo trackerToSolve;
        CoordinateSystemTrafo solveToTracker;
    };

    /** @brief Tracker state captured for one space-charge solve. */
    struct SpaceChargeStepState {
        std::size_t step       = 0;    ///< Global tracker step.
        double time            = 0.0;  ///< Tracker time in seconds.
        double timeStep        = 0.0;  ///< Tracker time step in seconds.
        bool emissionActive    = false;
        double emittedFraction = 1.0;  ///< Least-complete active source fraction in [0, 1].
        int mpiSize            = 1;    ///< Active communicator size.
        CoordinateFrameTransforms frames;
    };

    /**
     * @brief Per-container activity and tracker state for one space-charge call.
     *
     * @c trackingActive[i] corresponds to the ith particle container supplied at construction;
     * zero means inactive. The activity span is borrowed for the call.
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
