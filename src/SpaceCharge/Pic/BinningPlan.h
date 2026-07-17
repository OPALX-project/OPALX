/**
 * @file BinningPlan.h
 * @brief Declares persistent adaptive/fixed PIC bin iteration.
 */

#ifndef OPALX_SPACE_CHARGE_PIC_BINNING_PLAN_H
#define OPALX_SPACE_CHARGE_PIC_BINNING_PLAN_H

// Complete ParticleContainer field aliases before loading the legacy AdaptBins header.
#include "SpaceCharge/Pic/IterationPlan.h"

#include "PartBunch/Binning/AdaptBins.h"
#include "SpaceCharge/SelfFieldConfig.h"

#include <cstdint>
#include <memory>

namespace opalx::spacecharge {

    /**
     * @brief Reuses AdaptBins to prepare and lazily traverse globally nonempty merged bins.
     *
     * The plan persistently owns the concrete AdaptBins selected by immutable configuration and
     * borrows one particle container for its lifetime. The plan adds no particle-view cache;
     * AdaptBins reacquires its selector view during every full rebin.
     */
    template <typename T, unsigned Dim>
    class BinningPlan final : public IterationPlan<T, Dim> {
        static_assert(Dim == 3, "BinningPlan currently supports Dim == 3 only.");

    public:
        using Base              = IterationPlan<T, Dim>;
        using ParticleContainer = typename Base::ParticleContainer;
        using generation_type   = typename Base::generation_type;
        using Unit              = typename Base::Unit;
        using Selection         = typename Unit::Selection;
        using size_type         = typename Unit::size_type;
        using AdaptBins         = ParticleBinning::AdaptBinsBase<ParticleContainer>;
        using bin_index_type    = typename AdaptBins::bin_index_type;

        BinningPlan(ParticleContainer& particles, BinningConfig config);

        BinningPlan(const BinningPlan&)            = delete;
        BinningPlan& operator=(const BinningPlan&) = delete;
        BinningPlan(BinningPlan&&)                 = delete;
        BinningPlan& operator=(BinningPlan&&)      = delete;

        [[nodiscard]] IterationKind kind() const override { return IterationKind::Binning; }
        [[nodiscard]] const std::string& diagnosticName() const override { return config_m.name(); }
        [[nodiscard]] std::size_t maximumBinCount() const override {
            return config_m.maximumBins();
        }

        PreparedIteration prepare(
                generation_type particleGeneration,
                BinConfigurationObserver* observer = nullptr) override;

        [[nodiscard]] bool hasNext(
                const PreparedIteration& prepared,
                generation_type currentParticleGeneration) override;

        [[nodiscard]] std::optional<Unit> next(
                const PreparedIteration& prepared,
                generation_type currentParticleGeneration) override;

        /** @brief Copy the current global histogram into host-owned diagnostic storage. */
        [[nodiscard]] BinConfigurationSnapshot hostSnapshot() const;

    private:
        [[nodiscard]] std::unique_ptr<AdaptBins> makeBins();
        void recordIfRequested(
                BinConfigurationObserver* observer, BinConfigurationPoint point) const;
        void validatePrepared(
                const PreparedIteration& prepared, generation_type currentParticleGeneration) const;
        [[nodiscard]] bool stageNextNonemptyBin();
        void advanceEpoch();

        ParticleContainer& particles_m;
        BinningConfig config_m;
        std::unique_ptr<AdaptBins> bins_m;
        std::uint64_t preparationEpoch_m = 0;
        bin_index_type nextBin_m         = 0;
        size_type nextGlobalCount_m      = 0;
        bool nextBinReady_m              = false;
        bool prepared_m                  = false;
    };

    extern template class BinningPlan<double, 3>;

}  // namespace opalx::spacecharge

#include "SpaceCharge/Pic/BinningPlan.tpp"

#endif  // OPALX_SPACE_CHARGE_PIC_BINNING_PLAN_H
