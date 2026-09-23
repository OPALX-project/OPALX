/**
 * @file SpaceChargeSolveContext.cpp
 * @brief Implements construction of one-call space-charge contexts.
 */

#include "SpaceCharge/SpaceChargeSolveContext.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace opalx::spacecharge {

    SpaceChargeSolveContext::SpaceChargeSolveContext(
            std::span<const std::uint8_t> trackingActive, SpaceChargeStepState stepState)
        : trackingActive_m(trackingActive), stepState_m(std::move(stepState)) {
        if (!std::isfinite(stepState_m.time) || !std::isfinite(stepState_m.timeStep)) {
            throw std::invalid_argument(
                    "SpaceChargeSolveContext requires finite time and time step");
        }
        if (stepState_m.mpiSize < 1) {
            throw std::invalid_argument(
                    "SpaceChargeSolveContext communicator size must be positive");
        }
        if (!(stepState_m.emittedFraction >= 0.0 && stepState_m.emittedFraction <= 1.0)) {
            throw std::invalid_argument(
                    "SpaceChargeSolveContext emitted fraction must be in [0, 1]");
        }
    }

}  // namespace opalx::spacecharge
