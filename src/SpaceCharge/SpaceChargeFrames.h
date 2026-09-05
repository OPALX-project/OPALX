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

    /** @brief Apply a spatial transform to R and rotate normalized momentum P into the same axes.
     */
    inline void enterSolveFrame(
            const CoordinateFrameTransforms& frames, SpaceChargeParticleContainer& particles) {
        frames.trackerToSolve.transformBunchTo(particles.R.getView(), particles.getLocalNum());
        frames.trackerToSolve.rotateBunchTo(particles.P.getView(), particles.getLocalNum());
        particles.markMomentsDirty();
    }

    /** @brief Restore R/P and rotate computed E/B to tracker axes; this is not a Lorentz boost. */
    inline void leaveSolveFrame(
            const CoordinateFrameTransforms& frames, SpaceChargeParticleContainer& particles) {
        particles.transformBunch(frames.solveToTracker);
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_FRAMES_H
