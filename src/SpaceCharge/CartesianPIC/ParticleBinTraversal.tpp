/**
 * @file ParticleBinTraversal.tpp
 * @brief Implements persistent adaptive or fixed PIC bin traversal.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC_PARTICLE_BIN_TRAVERSAL_TPP
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC_PARTICLE_BIN_TRAVERSAL_TPP

#include "Utilities/OpalException.h"

#include <functional>
#include <limits>

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    ParticleBinTraversal<T, Dim>::ParticleBinTraversal(
            ParticleContainer& particles, BinningConfig config)
        : particles_m(particles), config_m(std::move(config)), bins_m(makeBins()) {}

    template <typename T, unsigned Dim>
    BinPreparationResult ParticleBinTraversal<T, Dim>::prepareBins(bool captureSnapshots) {
        prepared_m = false;
        nextBin_m  = 0;

        bins_m->doFullRebin(bins_m->getMaxBinCount());
        BinPreparationResult result;
        if (captureSnapshots) {
            // This snapshot represents the freshly rebuilt bins before adaptive merging.
            result.beforeMerge = hostSnapshot();
        }
        if (config_m.adaptive()) {
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

    template <typename T, unsigned Dim>
    std::optional<typename ParticleBinTraversal<T, Dim>::Unit>
    ParticleBinTraversal<T, Dim>::nextNonemptyBin() {
        if (!prepared_m) {
            throw OpalException(
                    "ParticleBinTraversal::nextNonemptyBin",
                    "The particle bins must be prepared before traversal.");
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
        const auto momentum           = particles_m.P.getView();

        ippl::Vector<double, Dim> localMomentumSum(0.0);
        Kokkos::parallel_reduce(
                "ParticleBinTraversal::meanMomentum", policy,
                KOKKOS_LAMBDA(const size_type index, ippl::Vector<double, Dim>& sum) {
                    sum += momentum(hash(index));
                },
                localMomentumSum);

        double localGammaSum = 0.0;
        Kokkos::parallel_reduce(
                "ParticleBinTraversal::meanGamma", policy,
                KOKKOS_LAMBDA(const size_type index, double& sum) {
                    const ippl::Vector<double, Dim> p = momentum(hash(index));
                    sum += Kokkos::sqrt(1.0 + p.dot(p));
                },
                localGammaSum);

        // Normalize only after the MPI reductions so every rank uses the same global mean
        // momentum and gamma for this bin's rest-frame solve.
        ippl::Comm->allreduce(localMomentumSum, 1, std::plus<ippl::Vector<double, Dim>>());
        ippl::Vector<double, Dim> globalMomentumSum = localMomentumSum;
        ippl::Comm->allreduce(localGammaSum, 1, std::plus<double>());

        std::array<double, Dim> meanMomentum{};
        for (unsigned dimension = 0; dimension < Dim; ++dimension) {
            meanMomentum[dimension] =
                    globalMomentumSum[dimension] / static_cast<double>(globalCount);
        }
        const double gamma = localGammaSum / static_cast<double>(globalCount);
        if (gamma <= 0.0) {
            throw OpalException(
                    "ParticleBinTraversal::nextNonemptyBin",
                    "Computed a non-positive mean gamma for a bin.");
        }

        const size_type localParticleCount = particles_m.getLocalNum();
        const bool coversAllLocalParticles =
                static_cast<size_type>(policy.begin()) == 0
                && static_cast<size_type>(policy.end()) == localParticleCount;
        return Unit{
                static_cast<std::size_t>(binIndex),
                localCount,
                globalCount,
                Selection::indexed(policy, hash),
                coversAllLocalParticles,
                meanMomentum,
                gamma};
    }

    template <typename T, unsigned Dim>
    BinConfigurationSnapshot ParticleBinTraversal<T, Dim>::hostSnapshot() const {
        BinConfigurationSnapshot snapshot;
        snapshot.lowerBound = bins_m->getBinConfigHost(snapshot.particleCounts, snapshot.widths);
        return snapshot;
    }

    template <typename T, unsigned Dim>
    std::unique_ptr<typename ParticleBinTraversal<T, Dim>::AdaptBins>
    ParticleBinTraversal<T, Dim>::makeBins() {
        if (config_m.maximumBins()
            > static_cast<std::size_t>(std::numeric_limits<bin_index_type>::max())) {
            throw OpalException(
                    "ParticleBinTraversal::ParticleBinTraversal",
                    "MAXBINS cannot be represented by the particle bin index type.");
        }
        const bin_index_type maximumBins = static_cast<bin_index_type>(config_m.maximumBins());

        using CoordinateSelector = ParticleBinning::CoordinateSelector<ParticleContainer>;
        using GammaSelector      = ParticleBinning::GammaSelector<ParticleContainer>;
        switch (config_m.parameter()) {
            case BinningVariable::VelocityZ:
                return std::make_unique<
                        ParticleBinning::AdaptBins<ParticleContainer, CoordinateSelector>>(
                        particles_m, CoordinateSelector(2), maximumBins, config_m.alpha(),
                        config_m.beta(), config_m.desiredWidth(), config_m.name());
            case BinningVariable::GammaZ:
                return std::make_unique<
                        ParticleBinning::AdaptBins<ParticleContainer, GammaSelector>>(
                        particles_m, GammaSelector(2), maximumBins, config_m.alpha(),
                        config_m.beta(), config_m.desiredWidth(), config_m.name());
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

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC_PARTICLE_BIN_TRAVERSAL_TPP
