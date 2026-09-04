/**
 * @file SpaceChargeFrameTransform.h
 * @brief Defines explicit tracker and solve-frame transformations.
 */

#ifndef OPALX_SPACE_CHARGE_FRAME_TRANSFORM_H
#define OPALX_SPACE_CHARGE_FRAME_TRANSFORM_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/SpaceChargeSolveContext.h"
#include "Utilities/OpalException.h"

namespace opalx::spacecharge {

    /**
     * @brief Transforms primary R, E, and B in the exact compatibility order.
     *
     * Views are reacquired for each operation, so layout migration between enter() and leave() is
     * safe. Callers invoke leave() explicitly after a successful solve. An exception is terminal
     * for the current run, and this object does not attempt rollback during stack unwinding.
     */
    template <typename T, unsigned Dim>
    class SpaceChargeFrameTransform final {
        static_assert(Dim == 3, "SpaceChargeFrameTransform currently supports Dim == 3 only.");

    public:
        using ParticleContainer = ::ParticleContainer<T, Dim>;

        SpaceChargeFrameTransform(
                const CoordinateFrameTransforms& frames, ParticleContainer& particles)
            : frames_m(frames), particles_m(particles) {}

        SpaceChargeFrameTransform(const SpaceChargeFrameTransform&)            = delete;
        SpaceChargeFrameTransform& operator=(const SpaceChargeFrameTransform&) = delete;
        SpaceChargeFrameTransform(SpaceChargeFrameTransform&&)                 = delete;
        SpaceChargeFrameTransform& operator=(SpaceChargeFrameTransform&&)      = delete;

        void enter() {
            if (positionsInSolveFrame_m || electricInSolveFrame_m || magneticInSolveFrame_m) {
                throw OpalException(
                        "SpaceChargeFrameTransform::enter",
                        "The particle frame transform is already active.");
            }
            frames_m.trackerToSolve.transformBunchTo(
                    particles_m.R.getView(), particles_m.getLocalNum());
            positionsInSolveFrame_m = true;
        }

        void markComputedFields() noexcept {
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

    private:
        const CoordinateFrameTransforms& frames_m;
        ParticleContainer& particles_m;
        bool positionsInSolveFrame_m = false;
        bool electricInSolveFrame_m  = false;
        bool magneticInSolveFrame_m  = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_FRAME_TRANSFORM_H
