/**
 * @file ParticleFieldSet.cpp
 * @brief Implements particle-set selection and structural validation.
 */

#include "SpaceCharge/ParticleFieldSet.h"

#include <stdexcept>

namespace opalx::spacecharge {

    ParticleFieldSet::ParticleFieldSet(
            std::span<ParticleFieldBinding3D> bindings, std::size_t primaryIndex)
        : bindings_m(bindings), primaryIndex_m(primaryIndex) {
        if (bindings_m.empty()) {
            throw std::invalid_argument("ParticleFieldSet requires at least one binding");
        }
        if (primaryIndex_m >= bindings_m.size()) {
            throw std::out_of_range("ParticleFieldSet primary index is outside the binding set");
        }
        for (const ParticleFieldBinding3D& binding : bindings_m) {
            if (!binding.hasCompleteIdentity()) {
                throw std::invalid_argument("ParticleFieldSet requires complete particle bindings");
            }
        }
    }

    void ParticleFieldSet::selectForSolve(ParticleSelectionMode mode) {
        for (std::size_t index = 0; index < bindings_m.size(); ++index) {
            bindings_m[index].selectedForSolve = mode == ParticleSelectionMode::PrimaryOnly
                                                         ? index == primaryIndex_m
                                                         : bindings_m[index].trackingActive;
        }
    }

}  // namespace opalx::spacecharge
