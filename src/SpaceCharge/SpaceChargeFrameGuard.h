/**
 * @file SpaceChargeFrameGuard.h
 * @brief Defines exception-safe tracker and solve-frame restoration.
 */

#ifndef OPALX_SPACE_CHARGE_FRAME_GUARD_H
#define OPALX_SPACE_CHARGE_FRAME_GUARD_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/SpaceChargeSolveContext.h"
#include "Utilities/OpalException.h"

namespace opalx::spacecharge {

    /**
     * @brief Restores primary R, E, and B in the exact compatibility order after a solve.
     *
     * Views are reacquired for each transform, so layout migration between enter() and leave()
     * is safe. The destructor suppresses cleanup failures to preserve an active solver exception.
     */
    template <typename T, unsigned Dim>
    class SpaceChargeFrameGuard final {
        static_assert(Dim == 3, "SpaceChargeFrameGuard currently supports Dim == 3 only.");

    public:
        using ParticleContainer = ::ParticleContainer<T, Dim>;

        SpaceChargeFrameGuard(const CoordinateFrameTransforms& frames, ParticleContainer& particles)
            : frames_m(frames), particles_m(particles) {}

        ~SpaceChargeFrameGuard() { restoreNoThrow(); }

        SpaceChargeFrameGuard(const SpaceChargeFrameGuard&)            = delete;
        SpaceChargeFrameGuard& operator=(const SpaceChargeFrameGuard&) = delete;
        SpaceChargeFrameGuard(SpaceChargeFrameGuard&&)                 = delete;
        SpaceChargeFrameGuard& operator=(SpaceChargeFrameGuard&&)      = delete;

        void enter() {
            if (positionsInSolveFrame_m || electricInSolveFrame_m || magneticInSolveFrame_m) {
                throw OpalException(
                        "SpaceChargeFrameGuard::enter",
                        "The particle frame guard is already active.");
            }
            // Set the state before launching so exception cleanup conservatively restores a
            // partially submitted device transformation.
            positionsInSolveFrame_m = true;
            frames_m.trackerToSolve.transformBunchTo(
                    particles_m.R.getView(), particles_m.getLocalNum());
        }

        void markComputedFields() noexcept {
            // A failing algorithm may have written either field before throwing. Treat both as
            // solve-frame data so cleanup rotates any partial result consistently.
            electricInSolveFrame_m = true;
            magneticInSolveFrame_m = true;
        }

        void leave() {
            // Preserve the established compatibility order: restore R first, then rotate E and B.
            // Each view is acquired here because domain migration may have replaced its storage.
            if (positionsInSolveFrame_m) {
                frames_m.solveToTracker.transformBunchTo(
                        particles_m.R.getView(), particles_m.getLocalNum());
                positionsInSolveFrame_m = false;
            }
            if (electricInSolveFrame_m) {
                frames_m.solveToTracker.rotateBunchTo(
                        particles_m.E.getView(), particles_m.getLocalNum());
                electricInSolveFrame_m = false;
            }
            if (magneticInSolveFrame_m) {
                frames_m.solveToTracker.rotateBunchTo(
                        particles_m.B.getView(), particles_m.getLocalNum());
                magneticInSolveFrame_m = false;
            }
        }

        void restoreNoThrow() noexcept {
            try {
                leave();
            } catch (...) {
                // Preserve the original algorithm or migration exception during stack unwinding.
            }
        }

        [[nodiscard]] bool positionsInSolveFrame() const noexcept {
            return positionsInSolveFrame_m;
        }

    private:
        const CoordinateFrameTransforms& frames_m;
        ParticleContainer& particles_m;
        bool positionsInSolveFrame_m = false;
        bool electricInSolveFrame_m  = false;
        bool magneticInSolveFrame_m  = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_FRAME_GUARD_H
