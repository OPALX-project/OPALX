/**
 * @file ParticleFieldSet.h
 * @brief Defines the native particle and field bindings used by one space-charge call.
 */

#ifndef OPALX_SPACE_CHARGE_PARTICLE_FIELD_SET_H
#define OPALX_SPACE_CHARGE_PARTICLE_FIELD_SET_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/SpaceChargeCapabilities.h"

#include <cstddef>
#include <span>

namespace opalx::spacecharge {

    /**
     * @brief Borrow the native particle storage and field attributes required by a solver.
     *
     * The binding owns no storage. Its pointers identify the exact container and attributes used
     * when the solver was constructed, while the two flags are mutable per-call selection state.
     */
    template <typename T, unsigned Dim>
    struct ParticleFieldBinding {
        using Container       = ::ParticleContainer<T, Dim>;
        using VectorAttribute = typename Container::particle_position_type;

        Container* container      = nullptr;
        VectorAttribute* position = nullptr;
        VectorAttribute* momentum = nullptr;
        VectorAttribute* electric = nullptr;
        VectorAttribute* magnetic = nullptr;
        bool trackingActive       = true;
        bool selectedForSolve     = true;

        [[nodiscard]] bool hasCompleteIdentity() const noexcept {
            return container != nullptr && position != nullptr && momentum != nullptr
                   && electric != nullptr && magnetic != nullptr;
        }

        [[nodiscard]] bool sameIdentity(const ParticleFieldBinding& other) const noexcept {
            return container == other.container && position == other.position
                   && momentum == other.momentum && electric == other.electric
                   && magnetic == other.magnetic;
        }
    };

    template <typename T, unsigned Dim>
    [[nodiscard]] ParticleFieldBinding<T, Dim> makeParticleFieldBinding(
            ::ParticleContainer<T, Dim>& container) {
        return {&container, &container.R, &container.P, &container.E, &container.B, true, true};
    }

    using ParticleFieldBinding3D = ParticleFieldBinding<double, 3>;

    /** @brief Non-owning collection of native particle bindings used by one space-charge solve. */
    class ParticleFieldSet {
    public:
        ParticleFieldSet(std::span<ParticleFieldBinding3D> bindings, std::size_t primaryIndex);

        [[nodiscard]] std::span<ParticleFieldBinding3D> bindings() const { return bindings_m; }
        [[nodiscard]] std::size_t primaryIndex() const { return primaryIndex_m; }

        /** @brief Apply an algorithm's selection mode without changing binding membership. */
        void selectForSolve(ParticleSelectionMode mode);

    private:
        std::span<ParticleFieldBinding3D> bindings_m;
        std::size_t primaryIndex_m = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PARTICLE_FIELD_SET_H
