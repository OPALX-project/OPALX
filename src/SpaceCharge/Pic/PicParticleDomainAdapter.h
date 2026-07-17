/**
 * @file PicParticleDomainAdapter.h
 * @brief Declares the private bridge between Cartesian PIC domains and particle storage.
 */

#ifndef OPALX_SPACE_CHARGE_PIC_PARTICLE_DOMAIN_ADAPTER_H
#define OPALX_SPACE_CHARGE_PIC_PARTICLE_DOMAIN_ADAPTER_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/ParticleSetView.h"

#include <cstddef>
#include <memory>
#include <span>
#include <vector>

namespace opalx::spacecharge {

    /** @brief Global particle envelope used to update a Cartesian PIC mesh. */
    struct PicDomainBounds {
        ippl::Vector<double, 3> lower{0.0};
        ippl::Vector<double, 3> upper{0.0};
    };

    /**
     * @brief Private host-side bridge to native 3D particle containers.
     *
     * The adapter borrows every container and never retains @c PartBunch or shared ownership.
     * Native particle-attribute objects have stable addresses, but their Kokkos views do not:
     * migration can replace or reorder storage. Consequently this class never caches a Kokkos
     * view. Callers of primaryPositions() must reacquire its view after generation() changes.
     *
     * The adapter is intentionally specific to Cartesian 3D PIC. Common self-field interfaces
     * remain independent of IPPL particle layouts and redistribution details.
     */
    class PicParticleDomainAdapter final {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;
        using ParticlePosition  = typename ParticleContainer::particle_position_type;
        using Layout            = ippl::FieldLayout<3>;
        using Mesh              = ippl::UniformCartesian<double, 3>;
        using size_type         = typename ParticleContainer::size_type;
        using generation_type   = ParticleSetView::generation_type;

        /**
         * @brief Borrow native containers from a non-owning pointer collection.
         * @param containers All particle containers in ParticleSetView order.
         * @param primaryIndex Explicit primary-container index.
         * @param generation Initial particle-storage generation.
         */
        explicit PicParticleDomainAdapter(
                std::span<ParticleContainer* const> containers, std::size_t primaryIndex = 0,
                generation_type generation = 0);

        /**
         * @brief Borrow pointees from the manager's shared-pointer collection.
         *
         * The shared pointers are inspected only during construction; the adapter stores raw
         * borrowed pointers and does not extend any container lifetime.
         */
        explicit PicParticleDomainAdapter(
                std::span<const std::shared_ptr<ParticleContainer>> containers,
                std::size_t primaryIndex = 0, generation_type generation = 0);

        PicParticleDomainAdapter(const PicParticleDomainAdapter&)            = delete;
        PicParticleDomainAdapter& operator=(const PicParticleDomainAdapter&) = delete;
        PicParticleDomainAdapter(PicParticleDomainAdapter&&)                 = delete;
        PicParticleDomainAdapter& operator=(PicParticleDomainAdapter&&)      = delete;

        /**
         * @brief Compute the global envelope of every globally nonempty container.
         *
         * If all containers are empty, this preserves the current behavior by asking the primary
         * container for its empty-distribution bounds.
         */
        [[nodiscard]] PicDomainBounds computeBounds();

        /**
         * @brief Return the stable primary position attribute, not its Kokkos view.
         *
         * Obtain a fresh view from the returned object after every generation change.
         */
        [[nodiscard]] ParticlePosition& primaryPositions();
        [[nodiscard]] const ParticlePosition& primaryPositions() const;
        [[nodiscard]] size_type primaryLocalCount() const;
        [[nodiscard]] size_type primaryTotalCount() const;

        /**
         * @brief Whether primary-only redistribution is unsafe for this particle set.
         *
         * A non-primary container blocks redistribution when it is still active in tracking or
         * contains particles. Tracking activity is deliberately independent of whether that
         * container is selected for the current self-field solve.
         */
        [[nodiscard]] bool redistributionBlocked(const ParticleSetView& particles) const;

        /**
         * @brief Collectively test the current primary load against @p threshold.
         * @param threshold Maximum per-rank deviation divided by the global primary count.
         * @param rankFlags Reusable host scratch with one entry per communicator rank.
         * @return True when at least one rank exceeds @p threshold.
         */
        [[nodiscard]] bool loadIsImbalanced(double threshold, std::span<int> rankFlags) const;

        /**
         * @brief Apply a rebuilt Cartesian layout to every container and migrate particles.
         *
         * Each native layout is updated, ParticleBase::update() performs migration, and moments
         * are marked dirty. The storage generation advances exactly once after all containers are
         * updated, and the borrowed ParticleSetView is advanced to the same generation. No view
         * obtained before this call may be reused afterward.
         */
        void updateLayoutsAndMigrate(Layout& fieldLayout, Mesh& mesh, ParticleSetView& particles);

        /** @brief Refresh distribution moments for every borrowed container. */
        void updateMoments();

        /** @brief Current monotonically increasing particle-storage generation. */
        [[nodiscard]] generation_type generation() const { return generation_m; }

    private:
        [[nodiscard]] ParticleContainer& primary();
        [[nodiscard]] const ParticleContainer& primary() const;
        void validateParticleSet(const ParticleSetView& particles) const;
        void advanceGeneration(ParticleSetView& particles);

        std::vector<ParticleContainer*> containers_m;
        std::size_t primaryIndex_m   = 0;
        generation_type generation_m = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC_PARTICLE_DOMAIN_ADAPTER_H
