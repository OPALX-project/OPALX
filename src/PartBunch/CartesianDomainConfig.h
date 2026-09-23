/**
 * @file CartesianDomainConfig.h
 * @brief Defines algorithm-neutral Cartesian particle-storage setup.
 */

#ifndef OPALX_PART_BUNCH_CARTESIAN_DOMAIN_CONFIG_H
#define OPALX_PART_BUNCH_CARTESIAN_DOMAIN_CONFIG_H

#include <array>
#include <cstddef>
#include <cstdint>

namespace opalx::spacecharge {

    enum class ParticleLayoutType : std::uint8_t { Spatial, SpatialOverlap };

    /**
     * @brief Immutable construction values for a Cartesian domain and particle layouts.
     *
     * The overlap cutoff is meaningful only for SpatialOverlap. The periodic flag controls
     * particle and field-layout wrapping, independently of later Poisson backend dispatch.
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
        bool periodicParticleBoundary = false;
        T boundingBoxIncreasePercent  = T(2);
    };

    using CartesianDomainConfig3D = CartesianDomainConfig<double, 3>;

}  // namespace opalx::spacecharge

#endif  // OPALX_PART_BUNCH_CARTESIAN_DOMAIN_CONFIG_H
