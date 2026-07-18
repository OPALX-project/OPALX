/**
 * @file WholeBunchPlan.h
 * @brief Declares direct whole-bunch PIC iteration.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_WHOLE_BUNCH_PLAN_H
#define OPALX_SPACE_CHARGE_PIC3D_WHOLE_BUNCH_PLAN_H

#include "SpaceCharge/Pic3D/IterationPlan.h"

#include <cstdint>

namespace opalx::spacecharge {

    /** @brief Emits exactly one direct solve unit for the borrowed particle container. */
    template <typename T, unsigned Dim>
    class WholeBunchPlan final : public IterationPlan<T, Dim> {
        static_assert(Dim == 3, "WholeBunchPlan currently supports Dim == 3 only.");

    public:
        using Base              = IterationPlan<T, Dim>;
        using ParticleContainer = typename Base::ParticleContainer;
        using generation_type   = typename Base::generation_type;
        using Unit              = typename Base::Unit;
        using Selection         = typename Unit::Selection;
        using size_type         = typename Unit::size_type;

        explicit WholeBunchPlan(ParticleContainer& particles) : particles_m(particles) {}

        WholeBunchPlan(const WholeBunchPlan&)            = delete;
        WholeBunchPlan& operator=(const WholeBunchPlan&) = delete;
        WholeBunchPlan(WholeBunchPlan&&)                 = delete;
        WholeBunchPlan& operator=(WholeBunchPlan&&)      = delete;

        [[nodiscard]] IterationKind kind() const override { return IterationKind::WholeBunch; }
        [[nodiscard]] const std::string& diagnosticName() const override {
            static const std::string emptyName;
            return emptyName;
        }
        [[nodiscard]] std::size_t maximumBinCount() const override { return 0; }

        PreparedIteration prepare(
                generation_type particleGeneration,
                BinConfigurationObserver* observer = nullptr) override;

        [[nodiscard]] bool hasNext(
                const PreparedIteration& prepared,
                generation_type currentParticleGeneration) override;

        [[nodiscard]] std::optional<Unit> next(
                const PreparedIteration& prepared,
                generation_type currentParticleGeneration) override;

    private:
        void validatePrepared(
                const PreparedIteration& prepared, generation_type currentParticleGeneration) const;
        void advanceEpoch();

        ParticleContainer& particles_m;
        std::uint64_t preparationEpoch_m = 0;
        bool unitPending_m               = false;
    };

    extern template class WholeBunchPlan<double, 3>;

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_WHOLE_BUNCH_PLAN_H
