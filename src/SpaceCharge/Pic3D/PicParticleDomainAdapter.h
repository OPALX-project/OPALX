/**
 * @file PicParticleDomainAdapter.h
 * @brief Declares the private bridge between Cartesian PIC domains and particle storage.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_PARTICLE_DOMAIN_ADAPTER_H
#define OPALX_SPACE_CHARGE_PIC3D_PARTICLE_DOMAIN_ADAPTER_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/ParticleSetView.h"

#include <cstddef>
#include <span>
#include <vector>

namespace opalx::spacecharge {

    struct PicDomainBounds {
        ippl::Vector<double, 3> lower{0.0};
        ippl::Vector<double, 3> upper{0.0};
    };

    /**
     * @brief Gives Cartesian PIC orchestration access to validated native particle storage.
     *
     * Kokkos views are reacquired after every layout update or migration. Binding identity is
     * checked once by SelfFieldSystem, so this private adapter does not duplicate membership or
     * lifetime bookkeeping.
     */
    class PicParticleDomainAdapter final {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;
        using ParticlePosition  = typename ParticleContainer::particle_position_type;
        using Layout            = ippl::FieldLayout<3>;
        using Mesh              = ippl::UniformCartesian<double, 3>;
        using size_type         = typename ParticleContainer::size_type;

        explicit PicParticleDomainAdapter(
                std::span<const ParticleFieldBinding3d> bindings, std::size_t primaryIndex = 0);

        PicParticleDomainAdapter(const PicParticleDomainAdapter&)            = delete;
        PicParticleDomainAdapter& operator=(const PicParticleDomainAdapter&) = delete;
        PicParticleDomainAdapter(PicParticleDomainAdapter&&)                 = delete;
        PicParticleDomainAdapter& operator=(PicParticleDomainAdapter&&)      = delete;

        [[nodiscard]] PicDomainBounds computeBounds();
        [[nodiscard]] ParticlePosition& primaryPositions();
        [[nodiscard]] const ParticlePosition& primaryPositions() const;
        [[nodiscard]] size_type primaryLocalCount() const;
        [[nodiscard]] size_type primaryTotalCount() const;

        [[nodiscard]] bool redistributionBlocked(const ParticleSetView& particles) const;
        [[nodiscard]] bool loadIsImbalanced(double threshold, std::span<int> rankFlags) const;

        /** @brief Rebind every particle layout and complete migration before views are reused. */
        void updateLayoutsAndMigrate(Layout& fieldLayout, Mesh& mesh);

        void updateMoments();

    private:
        [[nodiscard]] ParticleContainer& primary();
        [[nodiscard]] const ParticleContainer& primary() const;

        std::vector<ParticleContainer*> containers_m;
        std::size_t primaryIndex_m = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_PARTICLE_DOMAIN_ADAPTER_H
