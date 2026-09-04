/**
 * @file SpaceChargeSolveContext.cpp
 * @brief Implements construction of one-call space-charge contexts.
 */

#include "SpaceCharge/SpaceChargeSolveContext.h"

#include <stdexcept>
#include <utility>

namespace opalx::spacecharge {

    SpaceChargeSolveContext::SpaceChargeSolveContext(
            ParticleFieldSet particles, SpaceChargeStepState stepState, SpaceChargeRequest request)
        : particles_m(std::move(particles)),
          stepState_m(std::move(stepState)),
          request_m(std::move(request)) {
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
