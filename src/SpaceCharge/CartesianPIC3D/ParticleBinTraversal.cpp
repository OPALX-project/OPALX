/**
 * @file ParticleBinTraversal.cpp
 * @brief Implements persistent adaptive or fixed PIC bin traversal.
 */

#include "SpaceCharge/CartesianPIC3D/ParticleBinTraversal.h"
#include "Utilities/OpalException.h"

#include <functional>
#include <limits>

namespace opalx::spacecharge {

    ParticleBinTraversal::ParticleBinTraversal(ParticleContainer& particles, BinningConfig config)
        : particles_m(particles), config_m(std::move(config)), bins_m(makeBins()) {}

    BinPreparationResult ParticleBinTraversal::prepareBins(bool captureSnapshots) {
        prepared_m = false;
        nextBin_m  = 0;

        // A single bin carries no information in its limits, histogram, or sort order. The
        // sort would also serialize every particle on one atomic counter. Snapshots read the
        // histogram, so a requested dump keeps the full path.
        singleBin_m = config_m.maximumBins == 1 && !captureSnapshots;
        if (singleBin_m) {
            Kokkos::deep_copy(particles_m.Bin.getView(), 0);
            prepared_m = true;
            return {.mergedBinCount = 1};
        }

        bins_m->doFullRebin(bins_m->getMaxBinCount());
        BinPreparationResult result;
        if (captureSnapshots) {
            // This snapshot represents the freshly rebuilt bins before adaptive merging.
            result.beforeMerge = hostSnapshot();
        }
        if (config_m.adaptive) {
            bins_m->genAdaptiveHistogram();
            if (captureSnapshots) {
                // Fixed binning has no merge stage, so only adaptive traversal has this snapshot.
                result.afterMerge = hostSnapshot();
            }
        }
        bins_m->sortContainerByBin();
        result.mergedBinCount = static_cast<std::size_t>(bins_m->getCurrentBinCount());
        prepared_m            = true;
        return result;
    }

    std::optional<ParticleBinTraversal::Unit> ParticleBinTraversal::nextNonemptyBin() {
        if (!prepared_m) {
            throw OpalException(
                    "ParticleBinTraversal::nextNonemptyBin",
                    "The particle bins must be prepared before traversal.");
        }

        using MomentumView = ParticleContainer::particle_position_type::view_type;
        const MomentumView::const_type momentum = particles_m.P.getView();
        ippl::Vector<double, 4> momentSums(0.0);

        if (singleBin_m) {
            // The container count is the same on every rank, so all ranks return together.
            const size_type globalCount = particles_m.getTotalNum();
            if (nextBin_m > 0 || globalCount == 0) {
                return std::nullopt;
            }
            nextBin_m                  = 1;
            const size_type localCount = particles_m.getLocalNum();
            Kokkos::parallel_reduce(
                    "ParticleBinTraversal::singleBinMoments", localCount,
                    KOKKOS_LAMBDA(const size_type index, ippl::Vector<double, 4>& sum) {
                        const ippl::Vector<double, 3> p = momentum(index);
                        for (unsigned dimension = 0; dimension < 3; ++dimension) {
                            sum[dimension] += p[dimension];
                        }
                        sum[3] += Kokkos::sqrt(1.0 + p.dot(p));
                    },
                    momentSums);
            return makeUnit(
                    0, localCount, globalCount, Selection::direct(0, localCount), true,
                    momentSums);
        }

        const bin_index_type mergedBinCount = bins_m->getCurrentBinCount();
        size_type globalCount               = 0;
        // getNPartInBin(..., true) is collective. Every rank advances through empty bins in the
        // same order and returns the same next globally nonempty bin.
        while (nextBin_m < mergedBinCount) {
            globalCount = bins_m->getNPartInBin(nextBin_m, true);
            if (globalCount > 0) {
                break;
            }
            ++nextBin_m;
        }
        if (nextBin_m >= mergedBinCount) {
            return std::nullopt;
        }

        const bin_index_type binIndex = nextBin_m++;
        const size_type localCount    = bins_m->getNPartInBin(binIndex, false);
        const auto policy             = bins_m->getBinIterationPolicy(binIndex);
        const auto hash               = bins_m->getHashArray();
        const ParticleMeshFieldTransfer::Hash::const_type readOnlyHash = hash;

        Kokkos::parallel_reduce(
                "ParticleBinTraversal::moments", policy,
                KOKKOS_LAMBDA(const size_type index, ippl::Vector<double, 4>& sum) {
                    const ippl::Vector<double, 3> p = momentum(readOnlyHash(index));
                    for (unsigned dimension = 0; dimension < 3; ++dimension) {
                        sum[dimension] += p[dimension];
                    }
                    sum[3] += Kokkos::sqrt(1.0 + p.dot(p));
                },
                momentSums);

        const size_type localParticleCount = particles_m.getLocalNum();
        const bool coversAllLocalParticles =
                static_cast<size_type>(policy.begin()) == 0
                && static_cast<size_type>(policy.end()) == localParticleCount;
        return makeUnit(
                static_cast<std::size_t>(binIndex), localCount, globalCount,
                Selection::indexed(policy, hash), coversAllLocalParticles, momentSums);
    }

    ParticleBinTraversal::Unit ParticleBinTraversal::makeUnit(
            std::size_t ordinal, size_type localCount, size_type globalCount, Selection selection,
            bool coversAllLocalParticles, ippl::Vector<double, 4> momentSums) {
        // Normalize only after the MPI reduction so every rank uses the same global mean
        // momentum and gamma for this bin's rest-frame solve.
        ippl::Comm->allreduce(momentSums.begin(), 4, std::plus<double>());

        std::array<double, 3> meanMomentum{};
        for (unsigned dimension = 0; dimension < 3; ++dimension) {
            meanMomentum[dimension] = momentSums[dimension] / static_cast<double>(globalCount);
        }
        const double gamma = momentSums[3] / static_cast<double>(globalCount);
        if (gamma <= 0.0) {
            throw OpalException(
                    "ParticleBinTraversal::nextNonemptyBin",
                    "Computed a non-positive mean gamma for a bin.");
        }
        return Unit{ordinal,      localCount,  globalCount, selection, coversAllLocalParticles,
                    meanMomentum, gamma};
    }

    BinConfigurationSnapshot ParticleBinTraversal::hostSnapshot() const {
        BinConfigurationSnapshot snapshot;
        snapshot.lowerBound = bins_m->getBinConfigHost(snapshot.particleCounts, snapshot.widths);
        return snapshot;
    }

    std::unique_ptr<ParticleBinTraversal::AdaptBins> ParticleBinTraversal::makeBins() {
        if (config_m.maximumBins
            > static_cast<std::size_t>(std::numeric_limits<bin_index_type>::max())) {
            throw OpalException(
                    "ParticleBinTraversal::ParticleBinTraversal",
                    "MAXBINS cannot be represented by the particle bin index type.");
        }
        const bin_index_type maximumBins = static_cast<bin_index_type>(config_m.maximumBins);

        using CoordinateSelector = ParticleBinning::CoordinateSelector<ParticleContainer>;
        using GammaSelector      = ParticleBinning::GammaSelector<ParticleContainer>;
        switch (config_m.parameter) {
            case BinningVariable::VelocityZ:
                return std::make_unique<
                        ParticleBinning::AdaptBins<ParticleContainer, CoordinateSelector>>(
                        particles_m, CoordinateSelector(2), maximumBins, config_m.alpha,
                        config_m.beta, config_m.desiredWidth, config_m.name);
            case BinningVariable::GammaZ:
                return std::make_unique<
                        ParticleBinning::AdaptBins<ParticleContainer, GammaSelector>>(
                        particles_m, GammaSelector(2), maximumBins, config_m.alpha, config_m.beta,
                        config_m.desiredWidth, config_m.name);
            case BinningVariable::PositionZ:
                throw OpalException(
                        "ParticleBinTraversal::ParticleBinTraversal",
                        "POSITIONZ binning is not supported; use VELOCITYZ or GAMMAZ.");
            case BinningVariable::MomentumZ:
                throw OpalException(
                        "ParticleBinTraversal::ParticleBinTraversal",
                        "MOMENTUMZ binning is not supported; use VELOCITYZ or GAMMAZ.");
        }
        throw OpalException(
                "ParticleBinTraversal::ParticleBinTraversal", "Unknown binning parameter.");
    }

}  // namespace opalx::spacecharge
