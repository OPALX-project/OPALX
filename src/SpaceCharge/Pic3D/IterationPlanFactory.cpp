/**
 * @file IterationPlanFactory.cpp
 * @brief Implements immutable-config selection of Cartesian PIC iteration plans.
 */

#include "SpaceCharge/Pic3D/IterationPlanFactory.h"

#include "SpaceCharge/Pic3D/BinningPlan.h"
#include "SpaceCharge/Pic3D/WholeBunchPlan.h"

#include <memory>
#include <optional>
#include <utility>

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    std::unique_ptr<typename IterationPlanFactory<T, Dim>::Plan>
    IterationPlanFactory<T, Dim>::create(
            std::optional<BinningConfig> binning, ParticleContainer& particles) {
        if (binning.has_value()) {
            return std::make_unique<BinningPlan<T, Dim>>(particles, std::move(*binning));
        }
        return std::make_unique<WholeBunchPlan<T, Dim>>(particles);
    }

    template class IterationPlanFactory<double, 3>;

}  // namespace opalx::spacecharge
