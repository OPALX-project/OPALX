/**
 * @file CartesianDomainUpdater.cpp
 * @brief Implements Cartesian PIC geometry, migration, and redistribution updates.
 */

#include "SpaceCharge/CartesianPIC/CartesianDomainUpdater.h"

#include "SpaceCharge/Poisson/PoissonSolver.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <utility>

namespace opalx::spacecharge {
    namespace {

        PoissonFieldBinding poissonFields(CartesianDomainUpdater::FieldStorage& storage) {
            return {&storage.chargeDensity(), &storage.electricField()};
        }

    }  // namespace

    CartesianDomainUpdater::CartesianDomainUpdater(
            CartesianPICConfig config, std::span<ParticleContainer* const> particles)
        : config_m(std::move(config)),
          particles_m(particles.begin(), particles.end()),
          rankFlags_m(static_cast<std::size_t>(ippl::Comm->size()), 0) {
        if (particles_m.empty()
            || std::any_of(
                    particles_m.begin(), particles_m.end(), [](const ParticleContainer* particles) {
                        return particles == nullptr;
                    })) {
            throw OpalException(
                    "CartesianDomainUpdater::CartesianDomainUpdater",
                    "Cartesian PIC requires non-null particle containers.");
        }
    }

    bool CartesianDomainUpdater::updateForSolve(
            DomainCoordinateFrame frame, const SpaceChargeSolveContext& context,
            const CorrectionConfig& correction, const FixedDomain* fixedDomain,
            FieldStorage& fieldStorage, PoissonSolver& poissonSolver) {
        const bool beamFrame = frame == DomainCoordinateFrame::Beam;
        const bool fixed     = fixedDomain != nullptr;
        const FieldStorage::Extents extents =
                fixed ? config_m.grid.meshSize : targetExtents(correction);

        const bool repeatFieldLayoutRefresh = poissonRebuildRequired_m;
        if (fieldStorage.layoutExtents() != extents) {
            poissonRebuildRequired_m = true;
        }
        const std::array<bool, 3>& decomposition =
                fixed ? config_m.grid.decomposition : config_m.layoutDecomposition();
        const bool layoutChanged =
                fieldStorage.domain().rebuildGlobalLayoutInPlace(extents, decomposition);
        if (layoutChanged || repeatFieldLayoutRefresh) {
            fieldStorage.updateFieldLayoutsAfterLayoutChange();
        }

        CartesianBounds bounds;
        if (fixed) {
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                bounds.lower[dimension] = fixedDomain->lower[dimension];
                bounds.upper[dimension] = fixedDomain->upper[dimension];
            }
        } else {
            bounds = computeBounds(beamFrame);
            extendImageBounds(bounds, correction);
            expandBounds(
                    bounds, beamFrame && context.stepState().emissionActive,
                    context.stepState().emittedFraction, fieldStorage.layoutExtents()[2]);
        }
        updateMeshGeometry(bounds, fieldStorage);

        updateLayoutsAndMigrate(fieldStorage, beamFrame);
        updateMoments(beamFrame);

        bool redistributed = false;
        if (!fixed && beamFrame && isRedistributionDue(context.stepState().step)
            && loadIsImbalanced(config_m.loadBalancingThreshold)) {
            redistributed = redistribute(context, fieldStorage);
        }

        if (poissonRebuildRequired_m) {
            poissonSolver.rebuildAfterLayoutChange(poissonFields(fieldStorage));
            poissonRebuildRequired_m = false;
        }
        if (redistributed) {
            ippl::Comm->barrier();
        }
        return redistributed;
    }

    CartesianBounds CartesianDomainUpdater::computeBounds(bool primaryOnly) {
        CartesianBounds bounds;
        bool foundNonEmpty     = false;
        const std::size_t last = primaryOnly ? 1 : particles_m.size();
        for (std::size_t index = 0; index < last; ++index) {
            ParticleContainer& particles = *particles_m[index];
            if (particles.getTotalNum() == 0) {
                continue;
            }
            particles.computeMinMaxR();
            const auto lower = particles.getMinR();
            const auto upper = particles.getMaxR();
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

    void CartesianDomainUpdater::updateLayoutsAndMigrate(
            FieldStorage& fieldStorage, bool primaryOnly) {
        const std::size_t last = primaryOnly ? 1 : particles_m.size();
        for (std::size_t index = 0; index < last; ++index) {
            ParticleContainer& particles = *particles_m[index];
            particles.updateLayout(fieldStorage.layout(), fieldStorage.mesh());
            particles.update();
            particles.markMomentsDirty();
        }
    }

    void CartesianDomainUpdater::updateMoments(bool primaryOnly) {
        const std::size_t last = primaryOnly ? 1 : particles_m.size();
        for (std::size_t index = 0; index < last; ++index) {
            particles_m[index]->updateMoments();
        }
    }

    bool CartesianDomainUpdater::isRedistributionBlocked(
            std::span<const std::uint8_t> trackingActive) const {
        if (trackingActive.size() != particles_m.size()) {
            throw OpalException(
                    "CartesianDomainUpdater::isRedistributionBlocked",
                    "The container activity set has the wrong size.");
        }
        for (std::size_t index = 1; index < particles_m.size(); ++index) {
            if (trackingActive[index] != 0 || particles_m[index]->getTotalNum() > 0) {
                return true;
            }
        }
        return false;
    }

    bool CartesianDomainUpdater::loadIsImbalanced(double threshold) {
        std::fill(rankFlags_m.begin(), rankFlags_m.end(), 0);
        const std::size_t total = primary().getTotalNum();
        const int rankCount     = ippl::Comm->size();
        if (total == 0 || rankCount < 2) {
            return false;
        }
        const double equalShare = static_cast<double>(total) / static_cast<double>(rankCount);
        const double deviation = std::abs(static_cast<double>(primary().getLocalNum()) - equalShare)
                                 / static_cast<double>(total);
        rankFlags_m[static_cast<std::size_t>(ippl::Comm->rank())] = deviation > threshold ? 1 : 0;
        ippl::Comm->allreduce(rankFlags_m.data(), rankFlags_m.size(), std::plus<int>());
        return std::any_of(rankFlags_m.begin(), rankFlags_m.end(), [](int flag) {
            return flag != 0;
        });
    }

    CartesianDomainUpdater::FieldStorage::Extents CartesianDomainUpdater::targetExtents(
            const CorrectionConfig& correction) const {
        FieldStorage::Extents extents = config_m.grid.meshSize;
        if (correction.kind == SpaceChargeCorrectionType::ImageCharge) {
            if (extents[2] > std::numeric_limits<std::size_t>::max() / 2) {
                throw OpalException(
                        "CartesianDomainUpdater::targetExtents",
                        "The image-charge longitudinal mesh extent overflows size_t.");
            }
            extents[2] *= 2;
        }
        return extents;
    }

    void CartesianDomainUpdater::extendImageBounds(
            CartesianBounds& bounds, const CorrectionConfig& correction) const {
        if (correction.kind != SpaceChargeCorrectionType::ImageCharge) {
            return;
        }
        const double mirroredMinZ = 2.0 * correction.planeZ - bounds.upper[2];
        const double mirroredMaxZ = 2.0 * correction.planeZ - bounds.lower[2];
        bounds.lower[2]           = std::min(bounds.lower[2], mirroredMinZ);
        bounds.upper[2]           = std::max(bounds.upper[2], mirroredMaxZ);
    }

    void CartesianDomainUpdater::expandBounds(
            CartesianBounds& bounds, bool applyEmissionStretch, double emittedFraction,
            std::size_t longitudinalExtent) const {
        ippl::Vector<double, 3> span = bounds.upper - bounds.lower;
        for (unsigned dimension = 0; dimension < 3; ++dimension) {
            span[dimension] = std::max(span[dimension], 1.0e-6);
        }

        const double relativeIncrement = config_m.grid.boundingBoxIncreasePercent / 100.0;
        bounds.lower                   = bounds.lower - span * relativeIncrement;
        bounds.upper                   = bounds.upper + span * relativeIncrement;
        if (!applyEmissionStretch || longitudinalExtent <= 1) {
            return;
        }

        const double percent =
                std::max(1.0 / static_cast<double>(longitudinalExtent - 1), emittedFraction);
        const double originalLength =
                std::abs(bounds.upper[2] - bounds.lower[2]) / (1.0 + 2.0 * relativeIncrement);
        if (percent >= 1.0 || percent <= 0.0 || originalLength <= 0.0) {
            return;
        }
        bounds.upper[2] -= relativeIncrement * originalLength;
        bounds.lower[2]              = bounds.upper[2] - originalLength / percent;
        const double stretchedLength = originalLength / percent;
        bounds.upper[2] += relativeIncrement * stretchedLength;
        bounds.lower[2] -= relativeIncrement * stretchedLength;
    }

    void CartesianDomainUpdater::updateMeshGeometry(
            const CartesianBounds& bounds, FieldStorage& fieldStorage) const {
        const FieldStorage::Extents extents = fieldStorage.layoutExtents();
        const FieldStorage::Vector span     = bounds.upper - bounds.lower;
        FieldStorage::Vector spacing(0.0);
        FieldStorage::Vector origin = bounds.lower;
        for (unsigned dimension = 0; dimension < 3; ++dimension) {
            spacing[dimension] =
                    extents[dimension] <= 1
                            ? span[dimension]
                            : span[dimension] / static_cast<double>(extents[dimension] - 1);
            origin[dimension] = bounds.lower[dimension] - 0.5 * spacing[dimension];
        }
        fieldStorage.domain().setGeometry(bounds.lower, bounds.upper, spacing, origin);
    }

    bool CartesianDomainUpdater::isRedistributionDue(std::size_t step) const {
        const std::size_t frequency = config_m.repartitionFrequency;
        return frequency > 0 && step % frequency + 1 == frequency;
    }

    bool CartesianDomainUpdater::redistribute(
            const SpaceChargeSolveContext& context, FieldStorage& fieldStorage) {
        Inform m("CartesianDomainUpdater::redistribute");
        const int ranks = context.stepState().mpiSize;
        if (ranks < 2 || primary().getTotalNum() == 0) {
            return false;
        }
        if ((ranks & (ranks - 1)) != 0) {
            m << level1 << "Skipping ORB load balancing because ORB requires a power-of-two MPI "
              << "rank count; current rank count is " << ranks << "." << endl;
            return false;
        }
        if (isRedistributionBlocked(context.trackingActive())) {
            m << level2
              << "Skipping ORB load balancing because a non-primary particle container is "
                 "tracking-active or non-empty."
              << endl;
            return false;
        }

        const std::size_t localBefore = primary().getLocalNum();
        orb_m                         = Orb();
        orb_m.initialize(fieldStorage.layout(), fieldStorage.mesh(), fieldStorage.chargeDensity());
        poissonRebuildRequired_m = true;
        const bool succeeded = orb_m.binaryRepartition(primary().R, fieldStorage.layout(), false);

        // ORB may mutate the layout even when it reports ordinary failure.
        fieldStorage.updateFieldLayoutsAfterLayoutChange();
        updateLayoutsAndMigrate(fieldStorage, true);
        if (!succeeded) {
            m << level2 << "ORB load balancing failed; retaining the resulting valid layout."
              << endl;
            return false;
        }

        m << level2 << "ORB load balancing done. Rank 0: " << localBefore << " -> "
          << primary().getLocalNum() << " particles in primary container." << endl;
        return true;
    }

}  // namespace opalx::spacecharge
