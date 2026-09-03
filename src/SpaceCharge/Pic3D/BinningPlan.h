/**
 * @file BinningPlan.h
 * @brief Declares persistent adaptive or fixed PIC bin traversal.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_BINNING_PLAN_H
#define OPALX_SPACE_CHARGE_PIC3D_BINNING_PLAN_H

// Complete ParticleContainer field aliases before loading the legacy AdaptBins header.
#include "SpaceCharge/Pic3D/PicScatterGather.h"

#include "PartBunch/Binning/AdaptBins.h"
#include "SpaceCharge/SelfFieldConfig.h"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace opalx::spacecharge {

    struct BinConfigurationSnapshot final {
        using size_type = ippl::detail::size_type;

        double lowerBound = 0.0;
        std::vector<size_type> particleCounts;
        std::vector<double> widths;
    };

    struct BinningPreparation final {
        std::size_t mergedBinCount = 0;
        std::optional<BinConfigurationSnapshot> beforeMerge;
        std::optional<BinConfigurationSnapshot> afterMerge;
    };

    /** @brief One globally nonempty bin and its current particle selection. */
    template <typename T, unsigned Dim>
    struct BinnedSolveUnit final {
        using ScatterGather = PicScatterGather<T, Dim>;
        using Selection     = typename ScatterGather::Selection;
        using size_type     = typename ScatterGather::size_type;

        [[nodiscard]] Selection depositSelection() const {
            if (coversAllLocalParticles) {
                return Selection::direct(0, localParticleCount);
            }
            return indexedSelection;
        }

        std::size_t ordinal           = 0;
        size_type localParticleCount  = 0;
        size_type globalParticleCount = 0;
        Selection indexedSelection;
        bool coversAllLocalParticles = false;
        std::array<double, Dim> meanMomentum{};
        double gamma = 1.0;
    };

    /**
     * @brief Reuses AdaptBins to prepare and lazily traverse globally nonempty merged bins.
     *
     * prepare() owns the full rebin and optional host snapshots. next() performs one
     * rank-synchronous empty-bin query at a time and returns no external epoch or lifetime token.
     */
    template <typename T, unsigned Dim>
    class BinningPlan final {
        static_assert(Dim == 3, "BinningPlan currently supports Dim == 3 only.");

    public:
        using ParticleContainer = ::ParticleContainer<T, Dim>;
        using Unit              = BinnedSolveUnit<T, Dim>;
        using Selection         = typename Unit::Selection;
        using size_type         = typename Unit::size_type;
        using AdaptBins         = ParticleBinning::AdaptBinsBase<ParticleContainer>;
        using bin_index_type    = typename AdaptBins::bin_index_type;

        BinningPlan(ParticleContainer& particles, BinningConfig config);

        BinningPlan(const BinningPlan&)            = delete;
        BinningPlan& operator=(const BinningPlan&) = delete;
        BinningPlan(BinningPlan&&)                 = delete;
        BinningPlan& operator=(BinningPlan&&)      = delete;

        [[nodiscard]] const std::string& diagnosticName() const { return config_m.name(); }
        [[nodiscard]] std::size_t maximumBinCount() const { return config_m.maximumBins(); }

        [[nodiscard]] BinningPreparation prepare(bool captureSnapshots);
        [[nodiscard]] std::optional<Unit> next();

    private:
        [[nodiscard]] BinConfigurationSnapshot hostSnapshot() const;
        [[nodiscard]] std::unique_ptr<AdaptBins> makeBins();

        ParticleContainer& particles_m;
        BinningConfig config_m;
        std::unique_ptr<AdaptBins> bins_m;
        bin_index_type nextBin_m = 0;
        bool prepared_m          = false;
    };

    extern template class BinningPlan<double, 3>;

}  // namespace opalx::spacecharge

#include "SpaceCharge/Pic3D/BinningPlan.tpp"

#endif  // OPALX_SPACE_CHARGE_PIC3D_BINNING_PLAN_H
