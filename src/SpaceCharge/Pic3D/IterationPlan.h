/**
 * @file IterationPlan.h
 * @brief Declares particle-selection plans for Cartesian PIC solve units.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_ITERATION_PLAN_H
#define OPALX_SPACE_CHARGE_PIC3D_ITERATION_PLAN_H

#include "SpaceCharge/ParticleSetView.h"
#include "SpaceCharge/Pic3D/PicScatterGather.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace opalx::spacecharge {

    /** @brief Stable kind of particle traversal prepared for one solver call. */
    enum class IterationKind : std::uint8_t { WholeBunch, Binning };

    /** @brief How the executor composes the field produced for one solve unit. */
    enum class SolveUnitFieldMode : std::uint8_t { Direct, LorentzTransformed };

    /** @brief Point at which a host bin-configuration snapshot was captured. */
    enum class BinConfigurationPoint : std::uint8_t { BeforeMerge, AfterMerge };

    /** @brief Host-owned global bin data suitable for diagnostics and file output. */
    struct BinConfigurationSnapshot final {
        using size_type = ippl::detail::size_type;

        std::string diagnosticName;
        double lowerBound = 0.0;
        std::vector<size_type> particleCounts;
        std::vector<double> widths;
    };

    /** @brief Optional host callback for observing pre- and post-merge bin configurations. */
    class BinConfigurationObserver {
    public:
        virtual ~BinConfigurationObserver() = default;

        /** @brief Return whether a snapshot is required at @p point. */
        [[nodiscard]] virtual bool wants(BinConfigurationPoint point) const = 0;

        /** @brief Consume one requested host snapshot synchronously. */
        virtual void record(
                BinConfigurationPoint point, const BinConfigurationSnapshot& snapshot) = 0;
    };

    /** @brief Token proving that an iteration was prepared for one particle generation. */
    struct PreparedIteration final {
        using generation_type = ParticleSetView::generation_type;

        IterationKind kind                 = IterationKind::WholeBunch;
        generation_type particleGeneration = 0;
        std::size_t mergedBinCount         = 0;
        std::uint64_t preparationEpoch     = 0;
    };

    /**
     * @brief Host-side description of one particle selection and its global kinematics.
     *
     * The indexed selection borrows the current bin hash allocation. It remains usable only
     * until the owning plan is prepared again. No particle attribute or field view is retained.
     */
    template <typename T, unsigned Dim>
    struct SolveUnit final {
        using ScatterGather = PicScatterGather<T, Dim>;
        using Selection     = typename ScatterGather::Selection;
        using size_type     = typename ScatterGather::size_type;

        SolveUnit(
                IterationKind unitKind, SolveUnitFieldMode unitFieldMode, std::size_t unitOrdinal,
                size_type unitLocalCount, size_type unitGlobalCount, Selection unitIndexedSelection,
                bool unitCoversAllLocalParticles, std::array<double, Dim> unitMeanMomentum,
                double unitGamma)
            : kind(unitKind),
              fieldMode(unitFieldMode),
              ordinal(unitOrdinal),
              localParticleCount(unitLocalCount),
              globalParticleCount(unitGlobalCount),
              indexedSelection(std::move(unitIndexedSelection)),
              coversAllLocalParticles(unitCoversAllLocalParticles),
              meanMomentum(unitMeanMomentum),
              gamma(unitGamma) {}

        /** @brief Select direct deposition when this unit covers every local particle. */
        [[nodiscard]] Selection depositSelection() const {
            if (coversAllLocalParticles) {
                return Selection::direct(0, localParticleCount);
            }
            return indexedSelection;
        }

        IterationKind kind;
        SolveUnitFieldMode fieldMode;
        std::size_t ordinal;
        size_type localParticleCount;
        size_type globalParticleCount;
        Selection indexedSelection;
        bool coversAllLocalParticles;
        std::array<double, Dim> meanMomentum;
        double gamma;
    };

    /**
     * @brief Stateful host plan that lazily emits solve units for one prepared iteration.
     *
     * Concrete plans borrow a particle container but never own particle storage. prepare() must
     * follow every storage-generation change. hasNext() and next() reject a stale token before
     * accessing a retained bin hash or acquiring a current particle view.
     */
    template <typename T, unsigned Dim>
    class IterationPlan {
        static_assert(Dim == 3, "IterationPlan currently supports Dim == 3 only.");

    public:
        using ParticleContainer = ::ParticleContainer<T, Dim>;
        using generation_type   = ParticleSetView::generation_type;
        using Unit              = SolveUnit<T, Dim>;

        virtual ~IterationPlan() = default;

        /** @brief Return the immutable traversal selected at construction. */
        [[nodiscard]] virtual IterationKind kind() const = 0;

        /** @brief Stable diagnostic name; empty for whole-bunch traversal. */
        [[nodiscard]] virtual const std::string& diagnosticName() const = 0;

        /** @brief Configured maximum bin count; zero for whole-bunch traversal. */
        [[nodiscard]] virtual std::size_t maximumBinCount() const = 0;

        /** @brief Prepare selections for @p particleGeneration and return a validation token. */
        virtual PreparedIteration prepare(
                generation_type particleGeneration,
                BinConfigurationObserver* observer = nullptr) = 0;

        /**
         * @brief Return whether another solve unit is available for this preparation.
         *
         * Calls are rank-synchronous and may stage a global empty-bin query. Repeated calls before
         * next() are idempotent.
         */
        [[nodiscard]] virtual bool hasNext(
                const PreparedIteration& prepared, generation_type currentParticleGeneration) = 0;

        /** @brief Lazily emit the next solve unit, rejecting stale preparation tokens. */
        [[nodiscard]] virtual std::optional<Unit> next(
                const PreparedIteration& prepared, generation_type currentParticleGeneration) = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_ITERATION_PLAN_H
