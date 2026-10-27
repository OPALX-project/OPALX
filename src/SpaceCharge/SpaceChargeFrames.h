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
     *  @note Pol is rotated too, because leaveSolveFrame() rotates it back via transformBunch().
     */
    inline void enterSolveFrame(
            const CoordinateFrameTransforms& frames, SpaceChargeParticleContainer& particles) {
        frames.trackerToSolve.transformBunchTo(particles.R.getView(), particles.getLocalNum());
        frames.trackerToSolve.rotateBunchTo(particles.P.getView(), particles.getLocalNum());
        if (particles.hasSpin()) {
            frames.trackerToSolve.rotateBunchTo(particles.Pol.getView(), particles.getLocalNum());
        }
        particles.markMomentsDirty();
    }

    /**
     * @brief Restore R/P to tracker axes; this is not a Lorentz boost.
     * @param rotateFields Also rotate E/B. False when the solve already gathered them in tracker
     * axes, which saves two passes over the container.
     */
    inline void leaveSolveFrame(
            const CoordinateFrameTransforms& frames, SpaceChargeParticleContainer& particles,
            bool rotateFields = true) {
        particles.transformBunch(frames.solveToTracker, rotateFields);
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_FRAMES_H
