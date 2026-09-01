/**
 * @file ParticleLayoutConfig.h
 * @brief Algorithm-neutral setup contract for particle spatial layouts.
 */

#ifndef OPALX_SPACE_CHARGE_PARTICLE_LAYOUT_CONFIG_H
#define OPALX_SPACE_CHARGE_PARTICLE_LAYOUT_CONFIG_H

#include <cstdint>

namespace opalx::spacecharge {

    /** @brief Particle ownership layout selected before runtime solver construction. */
    enum class ParticleLayoutKind : std::uint8_t { Spatial, SpatialOverlap };

    /**
     * @brief Immutable setup values needed while constructing particle containers.
     *
     * The record deliberately contains no parser, mesh, field, or concrete solver type. The
     * overlap cutoff is meaningful only for SpatialOverlap layouts.
     */
    struct ParticleLayoutConfig {
        ParticleLayoutKind kind = ParticleLayoutKind::Spatial;
        double overlapCutoff    = 0.0;
        bool periodic           = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PARTICLE_LAYOUT_CONFIG_H
