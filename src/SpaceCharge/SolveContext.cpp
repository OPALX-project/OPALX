/**
 * @file SolveContext.cpp
 * @brief Implements construction of one-call self-field contexts.
 */

#include "SpaceCharge/SolveContext.h"

#include <stdexcept>
#include <utility>

namespace opalx::spacecharge {

    SolveContext::SolveContext(
            ParticleSetView particles, StepState stepState, RequestedPhysics requestedPhysics)
        : particles_m(std::move(particles)),
          stepState_m(std::move(stepState)),
          requestedPhysics_m(std::move(requestedPhysics)) {
        if (stepState_m.communicator.rank < 0
            || stepState_m.communicator.rank >= stepState_m.communicator.size) {
            throw std::invalid_argument("SolveContext communicator rank is outside its size");
        }
        if (stepState_m.communicator.size < 1) {
            throw std::invalid_argument("SolveContext communicator size must be positive");
        }
    }

}  // namespace opalx::spacecharge
