/**
 * @file SpaceChargeFrames.h
 * @brief Explicit particle transformations at the space-charge boundary.
 */

#ifndef OPALX_SPACE_CHARGE_FRAMES_H
#define OPALX_SPACE_CHARGE_FRAMES_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/SpaceChargeSolveContext.h"

namespace opalx::spacecharge {

    using SpaceChargeParticleContainer = ::ParticleContainer<double, 3>;

    /** @brief Overwrite existing self-field storage without allocating particle scratch. */
    inline void clearSelfFields(SpaceChargeParticleContainer& particles) {
        particles.E = Vector_t<double, 3>(0.0);
        particles.B = Vector_t<double, 3>(0.0);
    }

    /** @brief Move R, P and Pol into the solve frame; E and B are cleared, so not rotated. */
    inline void enterSolveFrame(
            const ParticleFrame& frame, SpaceChargeParticleContainer& particles) {
        particles.enterFrame(frame, true);
    }

    /** @brief Return to the reference frame, rotating the computed E/B; not a Lorentz boost. */
    inline void leaveSolveFrame(SpaceChargeParticleContainer& particles) {
        particles.enterReferenceFrame();
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_FRAMES_H
