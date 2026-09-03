/**
 * @file ParticleFrameGuard.h
 * @brief Defines exception-safe tracker and solve-frame restoration.
 */

#ifndef OPALX_SPACE_CHARGE_PARTICLE_FRAME_GUARD_H
#define OPALX_SPACE_CHARGE_PARTICLE_FRAME_GUARD_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/SolveContext.h"
#include "Utilities/OpalException.h"

namespace opalx::spacecharge {

    /**
     * @brief Restores primary R, E, and B in the exact compatibility order after a solve.
     *
     * Views are reacquired for each transform, so layout migration between enter() and leave()
     * is safe. The destructor suppresses cleanup failures to preserve an active solver exception.
     */
    template <typename T, unsigned Dim>
    class ParticleFrameGuard final {
        static_assert(Dim == 3, "ParticleFrameGuard currently supports Dim == 3 only.");

    public:
        using ParticleContainer = ::ParticleContainer<T, Dim>;

        ParticleFrameGuard(const FrameState& frames, ParticleContainer& particles)
            : frames_m(frames), particles_m(particles) {}

        ~ParticleFrameGuard() { restoreNoThrow(); }

        ParticleFrameGuard(const ParticleFrameGuard&)            = delete;
        ParticleFrameGuard& operator=(const ParticleFrameGuard&) = delete;
        ParticleFrameGuard(ParticleFrameGuard&&)                 = delete;
        ParticleFrameGuard& operator=(ParticleFrameGuard&&)      = delete;

        void enter() {
            if (positionsInSolveFrame_m || electricInSolveFrame_m || magneticInSolveFrame_m) {
                throw OpalException(
                        "ParticleFrameGuard::enter", "The particle frame guard is already active.");
            }
            positionsInSolveFrame_m = true;
            frames_m.trackerToSolve.transformBunchTo(
                    particles_m.R.getView(), particles_m.getLocalNum());
        }

        void markComputedFields() noexcept {
            electricInSolveFrame_m = true;
            magneticInSolveFrame_m = true;
        }

        void leave() {
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
            }
        }

        [[nodiscard]] bool positionsInSolveFrame() const noexcept {
            return positionsInSolveFrame_m;
        }

    private:
        const FrameState& frames_m;
        ParticleContainer& particles_m;
        bool positionsInSolveFrame_m = false;
        bool electricInSolveFrame_m  = false;
        bool magneticInSolveFrame_m  = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PARTICLE_FRAME_GUARD_H
