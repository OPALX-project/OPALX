/**
 * @file ParticleBinTraversal.h
 * @brief Adaptive and fixed PIC bin traversal.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC3D_PARTICLE_BIN_TRAVERSAL_H
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC3D_PARTICLE_BIN_TRAVERSAL_H

// Complete ParticleContainer field aliases before loading the legacy AdaptBins header.
#include "SpaceCharge/CartesianPIC3D/ParticleMeshFieldTransfer.h"

#include "PartBunch/Binning/AdaptBins.h"
#include "SpaceCharge/SpaceChargeConfig.h"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace opalx::spacecharge {

    /** @brief Host snapshot of bin bounds, widths, and particle counts. */
    struct BinConfigurationSnapshot final {
        using size_type = ippl::detail::size_type;

        double lowerBound = 0.0;
        std::vector<size_type> particleCounts;
        std::vector<double> widths;
    };

    /** @brief Prepared bin count and optional snapshots around adaptive merging. */
    struct BinPreparationResult final {
        std::size_t mergedBinCount = 0;
        std::optional<BinConfigurationSnapshot> beforeMerge;
        std::optional<BinConfigurationSnapshot> afterMerge;
    };

    /** @brief One globally nonempty bin; mean momentum uses Cartesian solve axes and beta-gamma
     * units. */
    struct ParticleBin final {
        using ScatterGather = ParticleMeshFieldTransfer;
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
        std::array<double, 3> meanMomentum{};
        double gamma = 1.0;
    };

    /**
     * @brief Prepares and traverses globally nonempty bins through AdaptBins.
     *
     * Call prepareBins() once per solve, then call nextNonemptyBin() collectively until it returns
     * @c std::nullopt.
     */
    class ParticleBinTraversal final {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;
        using Unit              = ParticleBin;
        using Selection         = typename Unit::Selection;
        using size_type         = typename Unit::size_type;
        using AdaptBins         = ParticleBinning::AdaptBinsBase<ParticleContainer>;
        using bin_index_type    = typename AdaptBins::bin_index_type;

        ParticleBinTraversal(ParticleContainer& particles, BinningConfig config);

        ParticleBinTraversal(const ParticleBinTraversal&)            = delete;
        ParticleBinTraversal& operator=(const ParticleBinTraversal&) = delete;
        ParticleBinTraversal(ParticleBinTraversal&&)                 = delete;
        ParticleBinTraversal& operator=(ParticleBinTraversal&&)      = delete;

        [[nodiscard]] const std::string& diagnosticName() const { return config_m.name; }
        [[nodiscard]] std::size_t maximumBinCount() const { return config_m.maximumBins; }

        /** @brief Rebuild bins, optionally capture host snapshots, and reset traversal. */
        [[nodiscard]] BinPreparationResult prepareBins(bool captureSnapshots);
        /** @brief Collectively return the next globally nonempty bin. */
        [[nodiscard]] std::optional<Unit> nextNonemptyBin();

    private:
        [[nodiscard]] BinConfigurationSnapshot hostSnapshot() const;
        [[nodiscard]] std::unique_ptr<AdaptBins> makeBins();

        ParticleContainer& particles_m;
        BinningConfig config_m;
        std::unique_ptr<AdaptBins> bins_m;
        bin_index_type nextBin_m = 0;
        bool prepared_m          = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC3D_PARTICLE_BIN_TRAVERSAL_H
