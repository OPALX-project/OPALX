/**
 * @file SpaceChargeSolveContext.h
 * @brief Defines all borrowed and per-call state needed for one space-charge solve.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H
#define OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H

#include "Algorithms/CoordinateSystemTrafo.h"
#include "SpaceCharge/ParticleFieldSet.h"

#include <cstddef>

namespace opalx::spacecharge {

    /** @brief Coordinate transforms owned by the current solve context. */
    struct CoordinateFrameTransforms {
        CoordinateSystemTrafo trackerToSolve;
        CoordinateSystemTrafo solveToTracker;
    };

    /** @brief Tracker state captured for one space-charge solve. */
    struct SpaceChargeStepState {
        std::size_t step       = 0;
        double time            = 0.0;
        double timeStep        = 0.0;
        bool emissionActive    = false;
        double emittedFraction = 1.0;
        int mpiSize            = 1;
        CoordinateFrameTransforms frames;
    };

    struct SpaceChargeCorrectionRequest {
        SpaceChargeCorrectionType kind = SpaceChargeCorrectionType::None;
        double planeZ                  = 0.0;
    };

    struct SpaceChargeRequest {
        bool useBinning     = false;
        bool writePotential = false;
        SpaceChargeCorrectionRequest correction;
    };

    /**
     * @brief Borrowed particle bindings and immutable state for one space-charge call.
     *
     * The context must not be retained after solve() returns. Persistent algorithm objects use the
     * native storage registered at construction and SpaceChargeSolver validates those identities
     * before every dispatch.
     */
    class SpaceChargeSolveContext {
    public:
        SpaceChargeSolveContext(
                ParticleFieldSet particles, SpaceChargeStepState stepState,
                SpaceChargeRequest request = {});

        [[nodiscard]] ParticleFieldSet& particles() { return particles_m; }
        [[nodiscard]] const ParticleFieldSet& particles() const { return particles_m; }
        [[nodiscard]] const SpaceChargeStepState& stepState() const { return stepState_m; }
        [[nodiscard]] const SpaceChargeRequest& request() const { return request_m; }

    private:
        ParticleFieldSet particles_m;
        SpaceChargeStepState stepState_m;
        SpaceChargeRequest request_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H
