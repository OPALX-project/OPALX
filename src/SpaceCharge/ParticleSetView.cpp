/**
 * @file ParticleSetView.cpp
 * @brief Implements particle-set selection and structural validation.
 */

#include "SpaceCharge/ParticleSetView.h"

#include <stdexcept>

namespace opalx::spacecharge {

    ParticleSetView::ParticleSetView(
            std::span<ParticleFieldBinding3d> bindings, std::size_t primaryIndex)
        : bindings_m(bindings), primaryIndex_m(primaryIndex) {
        if (bindings_m.empty()) {
            throw std::invalid_argument("ParticleSetView requires at least one binding");
        }
        if (primaryIndex_m >= bindings_m.size()) {
            throw std::out_of_range("ParticleSetView primary index is outside the binding set");
        }
        for (const ParticleFieldBinding3d& binding : bindings_m) {
            if (!binding.hasCompleteIdentity()) {
                throw std::invalid_argument("ParticleSetView requires complete particle bindings");
            }
        }
    }

    void ParticleSetView::applySelection(ParticleSelectionPolicy policy) {
        for (std::size_t index = 0; index < bindings_m.size(); ++index) {
            bindings_m[index].selectedForSolve = policy == ParticleSelectionPolicy::PrimaryOnly
                                                         ? index == primaryIndex_m
                                                         : bindings_m[index].trackingActive;
        }
    }

}  // namespace opalx::spacecharge
