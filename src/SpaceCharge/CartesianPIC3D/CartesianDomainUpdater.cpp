/**
 * @file CartesianDomainUpdater.cpp
 * @brief Implements CartesianPIC3D geometry, migration, and redistribution updates.
 */

#include "SpaceCharge/CartesianPIC3D/CartesianDomainUpdater.h"

#include "SpaceCharge/Poisson/PoissonSolver.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <utility>

namespace opalx::spacecharge {
    namespace {

        PoissonFieldBinding poissonFields(CartesianDomainUpdater::FieldStorage& storage) {
            return {&storage.chargeDensity(), &storage.electricField()};
        }

        /**
         * @brief Replace the CIC particle counts n on the owned vertices by the mixed ORB weight.
         *
         * w = (1 - pairWeight) n / sum(n) + pairWeight n^2 / sum(n^2), with both sums over all
         * ranks, so the weights add up to one. The P3M pair stage meets, for each particle, the
         * particles within the cutoff; their number follows the local density, and all mesh cells
         * have the same volume, so the pair count of a cell scales with n^2. Halo values are left
         * unchanged because ORB sums owned vertices only. Collective.
         *
         * @return This rank's share of the total weight; zero when no particle was deposited.
         */
        double mixPairWeight(Field<double, 3>& weights, double pairWeight) {
            const auto view   = weights.getView();
            const auto policy = ippl::getRangePolicy(view, weights.getNghost());

            double sums[2] = {0.0, 0.0};
            Kokkos::parallel_reduce(
                    "ORB particle count", policy,
                    KOKKOS_LAMBDA(
                            const std::size_t i, const std::size_t j, const std::size_t k,
                            double& count) { count += view(i, j, k); },
                    sums[0]);
            Kokkos::parallel_reduce(
                    "ORB pair count", policy,
                    KOKKOS_LAMBDA(
                            const std::size_t i, const std::size_t j, const std::size_t k,
                            double& pairs) {
                        const double n = view(i, j, k);
                        pairs += n * n;
                    },
                    sums[1]);
            ippl::Comm->allreduce(sums, 2, std::plus<double>());
            if (!(sums[0] > 0.0)) {
                return 0.0;
            }

            const double countScale = (1.0 - pairWeight) / sums[0];
            const double pairScale  = pairWeight / sums[1];
            double share            = 0.0;
            Kokkos::parallel_reduce(
                    "ORB mixed weight", policy,
                    KOKKOS_LAMBDA(
                            const std::size_t i, const std::size_t j, const std::size_t k,
                            double& localShare) {
                        const double n = view(i, j, k);
                        const double w = n * (countScale + pairScale * n);
                        view(i, j, k)  = w;
                        localShare += w;
                    },
                    share);
            return share;
        }

        /** @brief Largest rank value over the mean rank value. Collective. */
        double maxOverMean(double local) {
            double extremes[2] = {local, local};
            ippl::Comm->allreduce(&extremes[0], 1, std::greater<double>());
            ippl::Comm->allreduce(&extremes[1], 1, std::plus<double>());
            const double mean = extremes[1] / static_cast<double>(ippl::Comm->size());
            return mean > 0.0 ? extremes[0] / mean : 0.0;
        }

        /**
         * @brief Axes ORB must not cut: those the layout does not distribute over MPI ranks.
         *
         * The layout carries the decomposition it was last built with (PARFFTX/Y/Z, or the
         * longitudinal rebuild decomposition after a mesh resize), so ORB keeps that decomposition.
         */
        std::array<bool, 3> orbForbiddenAxes(const ippl::FieldLayout<3>& layout) {
            const std::array<bool, 3> parallel = layout.isParallel();
            return {!parallel[0], !parallel[1], !parallel[2]};
        }

    }  // namespace

    CartesianDomainUpdater::CartesianDomainUpdater(
            CartesianPIC3DConfig config, std::span<ParticleContainer* const> particles)
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
                    "CartesianPIC3D requires non-null particle containers.");
        }
    }

    bool CartesianDomainUpdater::updateForSolve(
            DomainCoordinateFrame frame, const SpaceChargeSolveContext& context,
            const DirichletPlaneConfig& dirichletPlane, const FixedDomain* fixedDomain,
            FieldStorage& fieldStorage, PoissonSolver& poissonSolver) {
        const bool beamFrame = frame == DomainCoordinateFrame::Beam;
        const bool fixed     = fixedDomain != nullptr;
        const FieldStorage::Extents extents =
                fixed ? config_m.grid.meshSize : targetExtents(dirichletPlane);

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
            extendImageBounds(bounds, dirichletPlane);
            expandBounds(
                    bounds, beamFrame && context.stepState().emissionActive,
                    context.stepState().emittedFraction, fieldStorage.layoutExtents()[2]);
        }
        updateMeshGeometry(bounds, fieldStorage);

        updateLayoutsAndMigrate(fieldStorage, frame);
        updateMoments(beamFrame);

        bool redistributed = false;
        if (!fixed && beamFrame && isRedistributionDue(context.stepState().step)) {
            redistributed = config_m.loadBalancingPairWeight > 0.0
                                    ? redistributeByPairWeight(context, fieldStorage)
                                    : loadIsImbalanced(config_m.loadBalancingThreshold)
                                              && redistribute(context, fieldStorage);
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
            FieldStorage& fieldStorage, DomainCoordinateFrame frame) {
        const bool primaryOnly = frame == DomainCoordinateFrame::Beam;
        // A mesh-tied P3M cutoff follows the spacing just set by updateMeshGeometry(). Passing it
        // with the layout update checks the new cutoff against the new rank regions only.
        const bool meshTiedCutoff = config_m.p3mCutoffCells > 0.0;
        double cutoff             = 0.0;
        if (meshTiedCutoff) {
            const auto& spacing = fieldStorage.mesh().getMeshSpacing();
            const double requested =
                    primaryOnly ? resolveOverlapCutoff<3>(0.0, config_m.p3mCutoffCells, spacing)
                                : config_m.p3mCutoffCells
                                          * std::min({spacing[0], spacing[1], spacing[2]});
            // An anisotropic mesh, such as the stretched longitudinal domain during emission,
            // can make the requested cutoff wider than a rank region. A smaller cutoff keeps the
            // P3M split exact but lets mesh errors grow, as alpha * h_max exceeds 2 / RCUT_CELLS.
            cutoff = std::min(requested, fieldStorage.domain().largestFittingOverlapCutoff());
            if (primaryOnly && cutoff < requested && !cutoffReductionReported_m) {
                Inform m("CartesianDomainUpdater::updateLayoutsAndMigrate");
                m << level1 << "Reducing the mesh-tied P3M cutoff from " << requested << " m to "
                  << cutoff << " m, half the smallest rank region; mesh accuracy is reduced while "
                  << "this lasts. Further reductions are not reported." << endl;
                cutoffReductionReported_m = true;
            }
        }
        const std::size_t last = primaryOnly ? 1 : particles_m.size();
        for (std::size_t index = 0; index < last; ++index) {
            ParticleContainer& particles = *particles_m[index];
            if (meshTiedCutoff && particles.hasP3MLayout()) {
                particles.updateLayout(fieldStorage.layout(), fieldStorage.mesh(), cutoff);
            } else {
                particles.updateLayout(fieldStorage.layout(), fieldStorage.mesh());
            }
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
        const std::size_t total = primary().getTotalNum();
        const int rankCount     = ippl::Comm->size();
        if (total == 0 || rankCount < 2) {
            return false;
        }
        const double equalShare = static_cast<double>(total) / static_cast<double>(rankCount);
        const double deviation = std::abs(static_cast<double>(primary().getLocalNum()) - equalShare)
                                 / static_cast<double>(total);
        return anyRankExceeds(deviation, threshold);
    }

    bool CartesianDomainUpdater::anyRankExceeds(double deviation, double threshold) {
        std::fill(rankFlags_m.begin(), rankFlags_m.end(), 0);
        rankFlags_m[static_cast<std::size_t>(ippl::Comm->rank())] = deviation > threshold ? 1 : 0;
        ippl::Comm->allreduce(rankFlags_m.data(), rankFlags_m.size(), std::plus<int>());
        return std::any_of(rankFlags_m.begin(), rankFlags_m.end(), [](int flag) {
            return flag != 0;
        });
    }

    CartesianDomainUpdater::FieldStorage::Extents CartesianDomainUpdater::targetExtents(
            const DirichletPlaneConfig& dirichletPlane) const {
        FieldStorage::Extents extents = config_m.grid.meshSize;
        if (dirichletPlane.kind == DirichletPlaneType::ImageCharge) {
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
            CartesianBounds& bounds, const DirichletPlaneConfig& dirichletPlane) const {
        if (dirichletPlane.kind != DirichletPlaneType::ImageCharge) {
            return;
        }
        const double mirroredMinZ = 2.0 * dirichletPlane.planeZ - bounds.upper[2];
        const double mirroredMaxZ = 2.0 * dirichletPlane.planeZ - bounds.lower[2];
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

    bool CartesianDomainUpdater::redistributionAllowed(
            const SpaceChargeSolveContext& context) const {
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
        return true;
    }

    void CartesianDomainUpdater::applyRepartition(FieldStorage& fieldStorage) {
        poissonRebuildRequired_m = true;
        fieldStorage.updateFieldLayoutsAfterLayoutChange();
        updateLayoutsAndMigrate(fieldStorage, DomainCoordinateFrame::Beam);
    }

    bool CartesianDomainUpdater::redistribute(
            const SpaceChargeSolveContext& context, FieldStorage& fieldStorage) {
        Inform m("CartesianDomainUpdater::redistribute");
        if (!redistributionAllowed(context)) {
            return false;
        }

        const std::size_t localBefore = primary().getLocalNum();
        auto& layout                  = fieldStorage.layout();
        orb_m                         = Orb();
        orb_m.initialize(layout, fieldStorage.mesh(), fieldStorage.chargeDensity());
        // ORB changes the layout only after its domains pass validation, so a failed
        // repartition leaves fields, particles, and the Poisson solver valid.
        if (!orb_m.binaryRepartition(primary().R, layout, false, orbForbiddenAxes(layout))) {
            m << level2 << "ORB load balancing failed along the distributed axes "
              << "(PARFFTX/Y/Z); keeping the current layout." << endl;
            return false;
        }

        applyRepartition(fieldStorage);

        m << level2 << "ORB load balancing done. Rank 0: " << localBefore << " -> "
          << primary().getLocalNum() << " particles in primary container, local domain "
          << layout.getLocalNDIndex() << "." << endl;
        return true;
    }

    double CartesianDomainUpdater::depositPairWeight(FieldStorage& fieldStorage) {
        // The weight field shares the charge-density storage; depositCharge() zeroes it again.
        orb_m = Orb();
        orb_m.initialize(fieldStorage.layout(), fieldStorage.mesh(), fieldStorage.chargeDensity());
        orb_m.scatterR(primary().R);
        return mixPairWeight(orb_m.bf_m, config_m.loadBalancingPairWeight);
    }

    bool CartesianDomainUpdater::redistributeByPairWeight(
            const SpaceChargeSolveContext& context, FieldStorage& fieldStorage) {
        Inform m("CartesianDomainUpdater::redistribute");
        const int ranks = ippl::Comm->size();
        if (ranks < 2 || primary().getTotalNum() == 0) {
            return false;
        }
        const double shareBefore = depositPairWeight(fieldStorage);
        const double deviation   = std::abs(shareBefore - 1.0 / static_cast<double>(ranks));
        if (!anyRankExceeds(deviation, config_m.loadBalancingThreshold)
            || !redistributionAllowed(context)) {
            return false;
        }

        const std::size_t localBefore  = primary().getLocalNum();
        const double weightRatioBefore = maxOverMean(shareBefore);
        const double countRatioBefore  = maxOverMean(static_cast<double>(localBefore));

        auto& layout           = fieldStorage.layout();
        const auto hostDomains = layout.getHostLocalDomains();
        const std::vector<ippl::NDIndex<3>> previousDomains(
                hostDomains.data(), hostDomains.data() + hostDomains.extent(0));
        const double fitBefore = fieldStorage.domain().largestFittingOverlapCutoff();

        // ORB changes the layout only after its domains pass validation, so a failed
        // repartition leaves fields, particles, and the Poisson solver valid.
        if (!orb_m.binaryRepartition(
                    primary().R, layout, /*isFirstRepartition=*/true, orbForbiddenAxes(layout))) {
            m << level2 << "ORB load balancing failed along the distributed axes "
              << "(PARFFTX/Y/Z); keeping the current layout." << endl;
            return false;
        }

        // Pair weights shrink the regions in the dense core. A cut that no longer fits the P3M
        // cutoff would reduce a mesh-tied cutoff for every rank or abort with a fixed one, so it
        // is replaced by the particle-count cut, which today's runs would have used.
        const double requestedCutoff = resolveOverlapCutoff<3>(
                config_m.p3mCutoff, config_m.p3mCutoffCells, fieldStorage.mesh().getMeshSpacing());
        const double fitAfter = fieldStorage.domain().largestFittingOverlapCutoff();
        if (fitAfter < requestedCutoff && fitAfter < fitBefore) {
            layout.updateLayout(previousDomains);
            if (!pairWeightFallbackReported_m) {
                m << level1 << "The pair-weighted ORB cut leaves rank regions too small for "
                  << "the P3M cutoff " << requestedCutoff << " m (largest fitting cutoff "
                  << fitAfter << " m); using the particle-count cut instead. Further "
                  << "fallbacks are not reported." << endl;
                pairWeightFallbackReported_m = true;
            }
            return redistribute(context, fieldStorage);
        }

        applyRepartition(fieldStorage);

        const double weightRatioAfter = maxOverMean(depositPairWeight(fieldStorage));
        const double countRatioAfter  = maxOverMean(static_cast<double>(primary().getLocalNum()));
        m << level2 << "ORB load balancing done with pair weight "
          << config_m.loadBalancingPairWeight << ". Max/mean weight " << weightRatioBefore << " -> "
          << weightRatioAfter << ", max/mean particle count " << countRatioBefore << " -> "
          << countRatioAfter << ". Rank 0: " << localBefore << " -> " << primary().getLocalNum()
          << " particles in primary container." << endl;
        return true;
    }

}  // namespace opalx::spacecharge
