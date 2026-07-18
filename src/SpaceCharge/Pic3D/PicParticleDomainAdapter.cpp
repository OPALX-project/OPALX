/**
 * @file PicParticleDomainAdapter.cpp
 * @brief Implements the private Cartesian PIC particle-domain bridge.
 */

#include "SpaceCharge/Pic3D/PicParticleDomainAdapter.h"

#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <typeindex>

namespace opalx::spacecharge {

    namespace {

        constexpr const char* adapterName = "PicParticleDomainAdapter";

    }  // namespace

    PicParticleDomainAdapter::PicParticleDomainAdapter(
            std::span<ParticleContainer* const> containers, std::size_t primaryIndex,
            generation_type generation)
        : containers_m(containers.begin(), containers.end()),
          primaryIndex_m(primaryIndex),
          generation_m(generation) {
        if (containers_m.empty()) {
            throw OpalException(adapterName, "At least one particle container is required.");
        }
        if (primaryIndex_m >= containers_m.size()) {
            throw OpalException(adapterName, "The primary particle-container index is invalid.");
        }
        if (std::any_of(
                    containers_m.begin(), containers_m.end(),
                    [](const ParticleContainer* container) {
                        return container == nullptr;
                    })) {
            throw OpalException(adapterName, "Null particle containers cannot be borrowed.");
        }
    }

    PicParticleDomainAdapter::PicParticleDomainAdapter(
            std::span<const std::shared_ptr<ParticleContainer>> containers,
            std::size_t primaryIndex, generation_type generation)
        : primaryIndex_m(primaryIndex), generation_m(generation) {
        containers_m.reserve(containers.size());
        for (const std::shared_ptr<ParticleContainer>& container : containers) {
            containers_m.push_back(container.get());
        }

        if (containers_m.empty()) {
            throw OpalException(adapterName, "At least one particle container is required.");
        }
        if (primaryIndex_m >= containers_m.size()) {
            throw OpalException(adapterName, "The primary particle-container index is invalid.");
        }
        if (std::any_of(
                    containers_m.begin(), containers_m.end(),
                    [](const ParticleContainer* container) {
                        return container == nullptr;
                    })) {
            throw OpalException(adapterName, "Null particle containers cannot be borrowed.");
        }
    }

    PicDomainBounds PicParticleDomainAdapter::computeBounds() {
        PicDomainBounds bounds;
        bool foundNonEmpty = false;

        for (ParticleContainer* container : containers_m) {
            if (container->getTotalNum() == 0) {
                continue;
            }

            // ParticleContainer reacquires R.getView() internally for this reduction.
            container->computeMinMaxR();
            const ippl::Vector<double, 3> lower = container->getMinR();
            const ippl::Vector<double, 3> upper = container->getMaxR();

            if (!foundNonEmpty) {
                bounds.lower  = lower;
                bounds.upper  = upper;
                foundNonEmpty = true;
                continue;
            }

            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                bounds.lower[dimension] = std::min(bounds.lower[dimension], lower[dimension]);
                bounds.upper[dimension] = std::max(bounds.upper[dimension], upper[dimension]);
            }
        }

        if (!foundNonEmpty) {
            // Preserve the established empty-bunch fallback. DistributionMoments maps its
            // all-empty reduction sentinels to a zero-sized box at the origin.
            primary().computeMinMaxR();
            bounds.lower = primary().getMinR();
            bounds.upper = primary().getMaxR();
        }

        return bounds;
    }

    PicParticleDomainAdapter::ParticlePosition& PicParticleDomainAdapter::primaryPositions() {
        return primary().R;
    }

    const PicParticleDomainAdapter::ParticlePosition& PicParticleDomainAdapter::primaryPositions()
            const {
        return primary().R;
    }

    PicParticleDomainAdapter::size_type PicParticleDomainAdapter::primaryLocalCount() const {
        return primary().getLocalNum();
    }

    PicParticleDomainAdapter::size_type PicParticleDomainAdapter::primaryTotalCount() const {
        return primary().getTotalNum();
    }

    bool PicParticleDomainAdapter::redistributionBlocked(const ParticleSetView& particles) const {
        validateParticleSet(particles);
        const std::span<ParticleContainerView> views = particles.containers();

        for (std::size_t i = 0; i < containers_m.size(); ++i) {
            if (i == primaryIndex_m) {
                continue;
            }
            if (views[i].trackingActive() || containers_m[i]->getTotalNum() > 0) {
                return true;
            }
        }
        return false;
    }

    bool PicParticleDomainAdapter::loadIsImbalanced(
            double threshold, std::span<int> rankFlags) const {
        const int rankCount = ippl::Comm->size();
        if (rankFlags.size() != static_cast<std::size_t>(rankCount)) {
            throw OpalException(
                    "PicParticleDomainAdapter::loadIsImbalanced",
                    "Rank-flag scratch must contain exactly one entry per communicator rank.");
        }

        std::fill(rankFlags.begin(), rankFlags.end(), 0);
        const size_type total = primaryTotalCount();
        if (total == 0 || rankCount < 2) {
            return false;
        }

        const double equalShare = static_cast<double>(total) / static_cast<double>(rankCount);
        const double deviation  = std::abs(static_cast<double>(primaryLocalCount()) - equalShare)
                                 / static_cast<double>(total);
        rankFlags[static_cast<std::size_t>(ippl::Comm->rank())] = deviation > threshold ? 1 : 0;
        ippl::Comm->allreduce(rankFlags.data(), rankFlags.size(), std::plus<int>());

        return std::any_of(rankFlags.begin(), rankFlags.end(), [](int flag) {
            return flag != 0;
        });
    }

    void PicParticleDomainAdapter::updateLayoutsAndMigrate(
            Layout& fieldLayout, Mesh& mesh, ParticleSetView& particles) {
        validateParticleSet(particles);

        for (ParticleContainer* container : containers_m) {
            container->getLayout().updateLayout(fieldLayout, mesh);
            container->update();
            container->markMomentsDirty();
        }

        advanceGeneration(particles);
    }

    void PicParticleDomainAdapter::updateMoments() {
        for (ParticleContainer* container : containers_m) {
            container->updateMoments();
        }
    }

    PicParticleDomainAdapter::ParticleContainer& PicParticleDomainAdapter::primary() {
        return *containers_m[primaryIndex_m];
    }

    const PicParticleDomainAdapter::ParticleContainer& PicParticleDomainAdapter::primary() const {
        return *containers_m[primaryIndex_m];
    }

    void PicParticleDomainAdapter::validateParticleSet(const ParticleSetView& particles) const {
        const std::span<ParticleContainerView> views = particles.containers();
        if (views.size() != containers_m.size() || particles.primaryIndex() != primaryIndex_m) {
            throw OpalException(
                    adapterName,
                    "ParticleSetView membership does not match the borrowed native containers.");
        }

        for (std::size_t i = 0; i < containers_m.size(); ++i) {
            const ParticleAttributeHandle* position = views[i].find(ParticleAttribute::Position);
            const bool wrongType =
                    position == nullptr
                    || position->nativeType() != std::type_index(typeid(ParticlePosition));
            if (wrongType
                || std::addressof(position->native<ParticlePosition>())
                           != std::addressof(containers_m[i]->R)) {
                throw OpalException(
                        adapterName,
                        "ParticleSetView order does not match the borrowed native containers.");
            }
        }
    }

    void PicParticleDomainAdapter::advanceGeneration(ParticleSetView& particles) {
        if (generation_m == std::numeric_limits<generation_type>::max()) {
            throw OpalException(adapterName, "Particle-storage generation counter overflowed.");
        }
        ++generation_m;
        particles.updateGeneration(generation_m);
    }

}  // namespace opalx::spacecharge
