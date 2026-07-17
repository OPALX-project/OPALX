/**
 * @file IterationPlanFactory.h
 * @brief Declares immutable-config selection of Cartesian PIC iteration plans.
 */

#ifndef OPALX_SPACE_CHARGE_PIC_ITERATION_PLAN_FACTORY_H
#define OPALX_SPACE_CHARGE_PIC_ITERATION_PLAN_FACTORY_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/Pic/IterationPlan.h"
#include "SpaceCharge/SelfFieldConfig.h"

#include <memory>
#include <optional>

namespace opalx::spacecharge {

    /** @brief Creates whole-bunch or binning traversal from immutable configuration. */
    template <typename T, unsigned Dim>
    class IterationPlanFactory final {
        static_assert(Dim == 3, "IterationPlanFactory currently supports Dim == 3 only.");

    public:
        using ParticleContainer = ::ParticleContainer<T, Dim>;
        using Plan              = IterationPlan<T, Dim>;

        [[nodiscard]] static std::unique_ptr<Plan> create(
                std::optional<BinningConfig> binning, ParticleContainer& particles);
    };

    extern template class IterationPlanFactory<double, 3>;

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC_ITERATION_PLAN_FACTORY_H
