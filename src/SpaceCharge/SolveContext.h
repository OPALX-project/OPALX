/**
 * @file SolveContext.h
 * @brief Defines all borrowed and per-call state needed for one self-field solve.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H
#define OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H

#include "Algorithms/CoordinateSystemTrafo.h"
#include "SpaceCharge/ParticleSetView.h"

#include <cstddef>

namespace opalx::spacecharge {

    /** @brief Coordinate transforms owned by the current solve context. */
    struct FrameState {
        CoordinateSystemTrafo trackerToSolve;
        CoordinateSystemTrafo solveToTracker;
    };

    /** @brief Tracker state captured for one self-field solve. */
    struct StepState {
        std::size_t step       = 0;
        double time            = 0.0;
        double timeStep        = 0.0;
        bool emissionActive    = false;
        double emittedFraction = 1.0;
        int mpiSize            = 1;
        FrameState frames;
    };

    struct CorrectionRequest {
        CorrectionKind kind = CorrectionKind::None;
        double planeZ       = 0.0;
    };

    struct RequestedPhysics {
        bool useBinning     = false;
        bool writePotential = false;
        CorrectionRequest correction;
    };

    /**
     * @brief Borrowed particle bindings and immutable state for one self-field call.
     *
     * The context must not be retained after execute() returns. Persistent solver objects use the
     * native storage registered at construction and the system validates those identities before
     * every dispatch.
     */
    class SolveContext {
    public:
        SolveContext(
                ParticleSetView particles, StepState stepState,
                RequestedPhysics requestedPhysics = {});

        [[nodiscard]] ParticleSetView& particles() { return particles_m; }
        [[nodiscard]] const ParticleSetView& particles() const { return particles_m; }
        [[nodiscard]] const StepState& stepState() const { return stepState_m; }
        [[nodiscard]] const RequestedPhysics& requestedPhysics() const {
            return requestedPhysics_m;
        }

    private:
        ParticleSetView particles_m;
        StepState stepState_m;
        RequestedPhysics requestedPhysics_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H
