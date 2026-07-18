/**
 * @file FieldFramePolicy.cpp
 * @brief Implements tracker-to-beam frame handling for Cartesian 3D PIC.
 */

#include "SpaceCharge/Pic3D/FieldFramePolicy.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Utilities/OpalException.h"

#include <typeindex>

namespace opalx::spacecharge {

    void FieldFramePolicy::validate(const SolveContext& context) const {
        const FrameState& frames = context.stepState().frames;
        if (!frames.trackerToSolve.has_value() || !frames.solveToTracker.has_value()) {
            throw OpalException(
                    "FieldFramePolicy::validate",
                    "Cartesian 3D PIC requires both per-call coordinate transforms.");
        }

        const std::type_index coordinateTransformType(typeid(CoordinateSystemTrafo));
        if (frames.trackerToSolve->nativeType() != coordinateTransformType
            || frames.solveToTracker->nativeType() != coordinateTransformType) {
            throw OpalException(
                    "FieldFramePolicy::validate",
                    "Cartesian 3D PIC coordinate handles must contain CoordinateSystemTrafo "
                    "objects.");
        }
    }

    void FieldFramePolicy::enter(
            const SolveContext& context, ParticleContainer& particles, State& state) const {
        validate(context);
        if (state.positionsInBeam || state.electricInBeam || state.magneticInBeam) {
            throw OpalException(
                    "FieldFramePolicy::enter", "The field-frame state is already active.");
        }

        const auto& trackerToBeam =
                context.stepState().frames.trackerToSolve->native<CoordinateSystemTrafo>();
        state.positionsInBeam = true;
        trackerToBeam.transformBunchTo(particles.R.getView(), particles.getLocalNum());
    }

    void FieldFramePolicy::markComputedFieldsInBeam(State& state) const noexcept {
        state.electricInBeam = true;
        state.magneticInBeam = true;
    }

    void FieldFramePolicy::leave(
            const SolveContext& context, ParticleContainer& particles, State& state) const {
        validate(context);
        const auto& beamToTracker =
                context.stepState().frames.solveToTracker->native<CoordinateSystemTrafo>();

        if (state.positionsInBeam) {
            beamToTracker.transformBunchTo(particles.R.getView(), particles.getLocalNum());
            state.positionsInBeam = false;
        }
        if (state.electricInBeam) {
            beamToTracker.rotateBunchTo(particles.E.getView(), particles.getLocalNum());
            state.electricInBeam = false;
        }
        if (state.magneticInBeam) {
            beamToTracker.rotateBunchTo(particles.B.getView(), particles.getLocalNum());
            state.magneticInBeam = false;
        }
    }

    void FieldFramePolicy::restore(
            const SolveContext& context, ParticleContainer& particles,
            State& state) const noexcept {
        try {
            leave(context, particles, state);
        } catch (...) {
        }
    }

}  // namespace opalx::spacecharge
