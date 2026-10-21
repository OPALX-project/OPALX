/**
 * @file CartesianDomainConfig.h
 * @brief Defines algorithm-neutral Cartesian particle-storage setup.
 */

#ifndef OPALX_PART_BUNCH_CARTESIAN_DOMAIN_CONFIG_H
#define OPALX_PART_BUNCH_CARTESIAN_DOMAIN_CONFIG_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>

namespace opalx::spacecharge {

    enum class ParticleLayoutType : std::uint8_t { Spatial, SpatialOverlap };

    /**
     * @brief Immutable construction values for a Cartesian domain and particle layouts.
     *
     * The overlap cutoff is meaningful only for SpatialOverlap. It is either fixed in metres
     * (overlapCutoff) or a multiple of the largest mesh spacing (overlapCutoffCells), see
     * resolveOverlapCutoff(). The periodic flag controls particle and field-layout wrapping,
     * independently of later Poisson backend dispatch.
     */
    template <typename T, unsigned Dim>
    struct CartesianDomainConfig {
        std::array<std::size_t, Dim> meshSize = [] {
            std::array<std::size_t, Dim> result{};
            result.fill(8);
            return result;
        }();
        std::array<bool, Dim> decomposition = [] {
            std::array<bool, Dim> result{};
            result.fill(true);
            return result;
        }();
        ParticleLayoutType layoutType = ParticleLayoutType::Spatial;
        T overlapCutoff               = T(0);
        T overlapCutoffCells          = T(0);
        bool periodicParticleBoundary = false;
        T boundingBoxIncreasePercent  = T(2);
    };

    using CartesianDomainConfig3D = CartesianDomainConfig<double, 3>;

    /**
     * @brief Overlap cutoff in metres for a mesh spacing.
     * @return cutoffCells times the largest spacing if cutoffCells is positive, otherwise
     * fixedCutoff.
     */
    template <unsigned Dim, typename T, typename Spacing>
    [[nodiscard]] T resolveOverlapCutoff(T fixedCutoff, T cutoffCells, const Spacing& spacing) {
        if (!(cutoffCells > T(0))) {
            return fixedCutoff;
        }
        T largestSpacing = T(0);
        for (unsigned dimension = 0; dimension < Dim; ++dimension) {
            largestSpacing = std::max(largestSpacing, static_cast<T>(spacing[dimension]));
        }
        return cutoffCells * largestSpacing;
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_PART_BUNCH_CARTESIAN_DOMAIN_CONFIG_H
