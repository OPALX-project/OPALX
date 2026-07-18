/**
 * @file BinningPlan.tpp
 * @brief Implements persistent adaptive/fixed PIC bin iteration.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_BINNING_PLAN_TPP
#define OPALX_SPACE_CHARGE_PIC3D_BINNING_PLAN_TPP

#include "Utilities/OpalException.h"

#include <array>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    BinningPlan<T, Dim>::BinningPlan(ParticleContainer& particles, BinningConfig config)
        : particles_m(particles), config_m(std::move(config)), bins_m(makeBins()) {}

    template <typename T, unsigned Dim>
    PreparedIteration BinningPlan<T, Dim>::prepare(
            generation_type particleGeneration, BinConfigurationObserver* observer) {
        advanceEpoch();
        prepared_m        = false;
        nextBin_m         = 0;
        nextGlobalCount_m = 0;
        nextBinReady_m    = false;

        bins_m->doFullRebin(bins_m->getMaxBinCount());
        recordIfRequested(observer, BinConfigurationPoint::BeforeMerge);
        if (config_m.adaptive()) {
            bins_m->genAdaptiveHistogram();
            recordIfRequested(observer, BinConfigurationPoint::AfterMerge);
        }
        bins_m->sortContainerByBin();

        const bin_index_type mergedBinCount = bins_m->getCurrentBinCount();
        prepared_m                          = true;
        return PreparedIteration{
                IterationKind::Binning, particleGeneration,
                static_cast<std::size_t>(mergedBinCount), preparationEpoch_m};
    }

    template <typename T, unsigned Dim>
    bool BinningPlan<T, Dim>::hasNext(
            const PreparedIteration& prepared, generation_type currentParticleGeneration) {
        validatePrepared(prepared, currentParticleGeneration);
        return stageNextNonemptyBin();
    }

    template <typename T, unsigned Dim>
    std::optional<typename BinningPlan<T, Dim>::Unit> BinningPlan<T, Dim>::next(
            const PreparedIteration& prepared, generation_type currentParticleGeneration) {
        validatePrepared(prepared, currentParticleGeneration);
        if (!stageNextNonemptyBin()) {
            return std::nullopt;
        }

        const bin_index_type binIndex = nextBin_m++;
        const size_type globalCount   = nextGlobalCount_m;
        nextGlobalCount_m             = 0;
        nextBinReady_m                = false;

        const size_type localCount = bins_m->getNPartInBin(binIndex, false);
        const auto policy          = bins_m->getBinIterationPolicy(binIndex);
        const auto hash            = bins_m->getHashArray();
        const auto momentum        = particles_m.P.getView();

        ippl::Vector<double, Dim> localMomentumSum(0.0);
        Kokkos::parallel_reduce(
                "BinningPlan::meanMomentum", policy,
                KOKKOS_LAMBDA(const size_type index, ippl::Vector<double, Dim>& sum) {
                    sum += momentum(hash(index));
                },
                localMomentumSum);

        double localGammaSum = 0.0;
        Kokkos::parallel_reduce(
                "BinningPlan::meanGamma", policy,
                KOKKOS_LAMBDA(const size_type index, double& sum) {
                    const ippl::Vector<double, Dim> p = momentum(hash(index));
                    sum += Kokkos::sqrt(1.0 + p.dot(p));
                },
                localGammaSum);

        ippl::Vector<double, Dim> globalMomentumSum(0.0);
        ippl::Comm->allreduce(localMomentumSum, 1, std::plus<ippl::Vector<double, Dim>>());
        globalMomentumSum = localMomentumSum;
        ippl::Comm->allreduce(localGammaSum, 1, std::plus<double>());

        std::array<double, Dim> meanMomentum{};
        for (unsigned dimension = 0; dimension < Dim; ++dimension) {
            meanMomentum[dimension] =
                    globalMomentumSum[dimension] / static_cast<double>(globalCount);
        }
        const double gamma = localGammaSum / static_cast<double>(globalCount);
        if (gamma <= 0.0) {
            throw OpalException(
                    "BinningPlan::next", "Computed a non-positive mean gamma for a bin.");
        }

        const size_type localParticleCount = particles_m.getLocalNum();
        const bool coversAllLocalParticles =
                static_cast<size_type>(policy.begin()) == 0
                && static_cast<size_type>(policy.end()) == localParticleCount;
        return Unit(
                IterationKind::Binning, SolveUnitFieldMode::LorentzTransformed,
                static_cast<std::size_t>(binIndex), localCount, globalCount,
                Selection::indexed(policy, hash), coversAllLocalParticles, meanMomentum, gamma);
    }

    template <typename T, unsigned Dim>
    BinConfigurationSnapshot BinningPlan<T, Dim>::hostSnapshot() const {
        BinConfigurationSnapshot snapshot;
        snapshot.diagnosticName = config_m.name();
        snapshot.lowerBound = bins_m->getBinConfigHost(snapshot.particleCounts, snapshot.widths);
        return snapshot;
    }

    template <typename T, unsigned Dim>
    std::unique_ptr<typename BinningPlan<T, Dim>::AdaptBins> BinningPlan<T, Dim>::makeBins() {
        if (config_m.maximumBins()
            > static_cast<std::size_t>(std::numeric_limits<bin_index_type>::max())) {
            throw OpalException(
                    "BinningPlan::BinningPlan",
                    "MAXBINS cannot be represented by the particle bin index type.");
        }
        const bin_index_type maximumBins = static_cast<bin_index_type>(config_m.maximumBins());

        using CoordinateSelector = ParticleBinning::CoordinateSelector<ParticleContainer>;
        using GammaSelector      = ParticleBinning::GammaSelector<ParticleContainer>;
        switch (config_m.parameter()) {
            case BinningParameterKind::VelocityZ:
                return std::make_unique<
                        ParticleBinning::AdaptBins<ParticleContainer, CoordinateSelector>>(
                        particles_m, CoordinateSelector(2), maximumBins, config_m.alpha(),
                        config_m.beta(), config_m.desiredWidth(), config_m.name());
            case BinningParameterKind::GammaZ:
                return std::make_unique<
                        ParticleBinning::AdaptBins<ParticleContainer, GammaSelector>>(
                        particles_m, GammaSelector(2), maximumBins, config_m.alpha(),
                        config_m.beta(), config_m.desiredWidth(), config_m.name());
            case BinningParameterKind::PositionZ:
                throw OpalException(
                        "BinningPlan::BinningPlan",
                        "POSITIONZ binning is not supported; use VELOCITYZ or GAMMAZ.");
            case BinningParameterKind::MomentumZ:
                throw OpalException(
                        "BinningPlan::BinningPlan",
                        "MOMENTUMZ binning is not supported; use VELOCITYZ or GAMMAZ.");
        }
        throw OpalException("BinningPlan::BinningPlan", "Unknown binning parameter.");
    }

    template <typename T, unsigned Dim>
    void BinningPlan<T, Dim>::recordIfRequested(
            BinConfigurationObserver* observer, BinConfigurationPoint point) const {
        if (observer == nullptr || !observer->wants(point)) {
            return;
        }
        const BinConfigurationSnapshot snapshot = hostSnapshot();
        observer->record(point, snapshot);
    }

    template <typename T, unsigned Dim>
    void BinningPlan<T, Dim>::validatePrepared(
            const PreparedIteration& prepared, generation_type currentParticleGeneration) const {
        if (!prepared_m || prepared.kind != IterationKind::Binning
            || prepared.preparationEpoch != preparationEpoch_m) {
            throw OpalException(
                    "BinningPlan::next",
                    "The prepared iteration does not belong to the current binning epoch.");
        }
        if (prepared.particleGeneration != currentParticleGeneration) {
            throw OpalException(
                    "BinningPlan::next",
                    "Particle storage changed after the bin iteration was prepared.");
        }
        if (prepared.mergedBinCount != static_cast<std::size_t>(bins_m->getCurrentBinCount())) {
            throw OpalException("BinningPlan::next", "The prepared merged-bin count is stale.");
        }
    }

    template <typename T, unsigned Dim>
    bool BinningPlan<T, Dim>::stageNextNonemptyBin() {
        if (nextBinReady_m) {
            return true;
        }

        const bin_index_type mergedBinCount = bins_m->getCurrentBinCount();
        while (nextBin_m < mergedBinCount) {
            nextGlobalCount_m = bins_m->getNPartInBin(nextBin_m, true);
            if (nextGlobalCount_m > 0) {
                nextBinReady_m = true;
                return true;
            }
            ++nextBin_m;
        }
        nextGlobalCount_m = 0;
        return false;
    }

    template <typename T, unsigned Dim>
    void BinningPlan<T, Dim>::advanceEpoch() {
        if (preparationEpoch_m == std::numeric_limits<std::uint64_t>::max()) {
            throw OpalException(
                    "BinningPlan::prepare", "The iteration preparation epoch overflowed.");
        }
        ++preparationEpoch_m;
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_BINNING_PLAN_TPP
