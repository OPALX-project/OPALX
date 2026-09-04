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

    inline void enterSolveFrame(
            const CoordinateFrameTransforms& frames, SpaceChargeParticleContainer& particles) {
        frames.trackerToSolve.transformBunchTo(particles.R.getView(), particles.getLocalNum());
    }

    /** @brief Restore R, E, and B in the established compatibility order. */
    inline void leaveSolveFrame(
            const CoordinateFrameTransforms& frames, SpaceChargeParticleContainer& particles) {
        frames.solveToTracker.transformBunchTo(particles.R.getView(), particles.getLocalNum());
        frames.solveToTracker.rotateBunchTo(particles.E.getView(), particles.getLocalNum());
        frames.solveToTracker.rotateBunchTo(particles.B.getView(), particles.getLocalNum());
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_FRAMES_H
