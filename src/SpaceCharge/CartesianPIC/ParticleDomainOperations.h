/**
 * @file ParticleDomainOperations.h
 * @brief Declares the private bridge between Cartesian PIC domains and particle storage.
 */

#ifndef OPALX_SPACE_CHARGE_PARTICLE_DOMAIN_OPERATIONS_H
#define OPALX_SPACE_CHARGE_PARTICLE_DOMAIN_OPERATIONS_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/ParticleFieldSet.h"

#include <cstddef>
#include <span>
#include <vector>

namespace opalx::spacecharge {

    struct CartesianBounds {
        ippl::Vector<double, 3> lower{0.0};
        ippl::Vector<double, 3> upper{0.0};
    };

    /**
     * @brief Gives Cartesian PIC orchestration access to validated native particle storage.
     *
     * Kokkos views are reacquired after every layout update or migration. Binding identity is
     * checked once by SpaceChargeSolver, so these operations do not duplicate membership or
     * lifetime bookkeeping.
     */
    class ParticleDomainOperations final {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;
        using ParticlePosition  = typename ParticleContainer::particle_position_type;
        using Layout            = ippl::FieldLayout<3>;
        using Mesh              = ippl::UniformCartesian<double, 3>;
        using size_type         = typename ParticleContainer::size_type;

        explicit ParticleDomainOperations(
                std::span<const ParticleFieldBinding3D> bindings, std::size_t primaryIndex = 0);

        ParticleDomainOperations(const ParticleDomainOperations&)            = delete;
        ParticleDomainOperations& operator=(const ParticleDomainOperations&) = delete;
        ParticleDomainOperations(ParticleDomainOperations&&)                 = delete;
        ParticleDomainOperations& operator=(ParticleDomainOperations&&)      = delete;

        /** @brief Compute the global envelope of all globally nonempty containers. */
        [[nodiscard]] CartesianBounds computeBounds();
        [[nodiscard]] ParticlePosition& primaryPositions();
        [[nodiscard]] const ParticlePosition& primaryPositions() const;
        [[nodiscard]] size_type primaryLocalCount() const;
        [[nodiscard]] size_type primaryTotalCount() const;

        /**
         * @brief Whether primary-only ORB would be unsafe for the current bunch.
         *
         * A non-primary container blocks redistribution when it is tracking-active or nonempty,
         * independently of whether the current space-charge algorithm selected it for solving.
         */
        [[nodiscard]] bool isRedistributionBlocked(const ParticleFieldSet& particles) const;
        [[nodiscard]] bool loadIsImbalanced(double threshold, std::span<int> rankFlags) const;

        /** @brief Rebind every particle layout and complete migration before views are reused. */
        void updateLayoutsAndMigrate(Layout& fieldLayout, Mesh& mesh);

        /** @brief Rebind and migrate only the primary container used by Cartesian PIC. */
        void updatePrimaryLayoutAndMigrate(Layout& fieldLayout, Mesh& mesh);

        void updateMoments();
        void updatePrimaryMoments();

    private:
        [[nodiscard]] ParticleContainer& primary();
        [[nodiscard]] const ParticleContainer& primary() const;

        std::vector<ParticleContainer*> containers_m;
        std::size_t primaryIndex_m = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PARTICLE_DOMAIN_OPERATIONS_H
