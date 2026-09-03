/**
 * @file ParticleStorageConfig.h
 * @brief Defines algorithm-neutral Cartesian particle-storage setup.
 */

#ifndef OPALX_SPACE_CHARGE_PARTICLE_STORAGE_CONFIG_H
#define OPALX_SPACE_CHARGE_PARTICLE_STORAGE_CONFIG_H

#include <array>
#include <cstddef>
#include <cstdint>

namespace opalx::spacecharge {

    enum class ParticleLayoutKind : std::uint8_t { Spatial, SpatialOverlap };

    /**
     * @brief Immutable construction values for a Cartesian domain and particle layouts.
     *
     * The overlap cutoff is meaningful only for SpatialOverlap. The periodic flag controls
     * particle and field-layout wrapping, independently of later Poisson backend dispatch.
     */
    template <typename T, unsigned Dim>
    struct ParticleStorageConfig {
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
        ParticleLayoutKind layoutKind = ParticleLayoutKind::Spatial;
        T overlapCutoff               = T(0);
        bool periodicParticleBoundary = false;
        T boundingBoxIncreasePercent  = T(2);
    };

    using ParticleStorageConfig3d = ParticleStorageConfig<double, 3>;

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PARTICLE_STORAGE_CONFIG_H
