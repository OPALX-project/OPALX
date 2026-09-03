/**
 * @file PicParticleDomainAdapter.cpp
 * @brief Implements the private Cartesian PIC particle-domain bridge.
 */

#include "SpaceCharge/Pic3D/PicParticleDomainAdapter.h"

#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>
#include <functional>

namespace opalx::spacecharge {

    PicParticleDomainAdapter::PicParticleDomainAdapter(
            std::span<const ParticleFieldBinding3d> bindings, std::size_t primaryIndex)
        : primaryIndex_m(primaryIndex) {
        containers_m.reserve(bindings.size());
        for (const ParticleFieldBinding3d& binding : bindings) {
            containers_m.push_back(binding.container);
        }
        if (containers_m.empty() || primaryIndex_m >= containers_m.size()) {
            throw OpalException("PicParticleDomainAdapter", "The particle binding set is invalid.");
        }
        if (std::any_of(
                    containers_m.begin(), containers_m.end(),
                    [](const ParticleContainer* container) {
                        return container == nullptr;
                    })) {
            throw OpalException(
                    "PicParticleDomainAdapter", "Null particle containers cannot be borrowed.");
        }
    }

    PicDomainBounds PicParticleDomainAdapter::computeBounds() {
        PicDomainBounds bounds;
        bool foundNonEmpty = false;
        for (ParticleContainer* container : containers_m) {
            if (container->getTotalNum() == 0) {
                continue;
            }
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
        const auto bindings = particles.bindings();
        for (std::size_t index = 0; index < containers_m.size(); ++index) {
            if (index == primaryIndex_m) {
                continue;
            }
            if (bindings[index].trackingActive || containers_m[index]->getTotalNum() > 0) {
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

    void PicParticleDomainAdapter::updateLayoutsAndMigrate(Layout& fieldLayout, Mesh& mesh) {
        for (ParticleContainer* container : containers_m) {
            container->updateLayout(fieldLayout, mesh);
            container->update();
            container->markMomentsDirty();
        }
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

}  // namespace opalx::spacecharge
