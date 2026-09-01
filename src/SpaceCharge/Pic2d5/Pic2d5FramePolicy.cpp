/**
 * @file Pic2d5FramePolicy.cpp
 * @brief Implements exception-safe FFT2D5 compatibility-frame restoration.
 */

#include "SpaceCharge/Pic2d5/Pic2d5FramePolicy.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Utilities/OpalException.h"

#include <typeindex>

namespace opalx::spacecharge {

    void Pic2d5FramePolicy::validate(const SolveContext& context) const {
        const FrameState& frames = context.stepState().frames;
        if (!frames.trackerToSolve.has_value() || !frames.solveToTracker.has_value()) {
            throw OpalException(
                    "Pic2d5FramePolicy::validate",
                    "FFT2D5 requires both per-call coordinate transforms.");
        }
        const std::type_index transformType(typeid(CoordinateSystemTrafo));
        if (frames.trackerToSolve->nativeType() != transformType
            || frames.solveToTracker->nativeType() != transformType) {
            throw OpalException(
                    "Pic2d5FramePolicy::validate",
                    "FFT2D5 coordinate handles must contain CoordinateSystemTrafo objects.");
        }
    }

    void Pic2d5FramePolicy::enter(
            const SolveContext& context, ParticleContainer& primary, State& state) const {
        validate(context);
        if (state.primaryPositionsInSolveFrame || state.primaryElectricInSolveFrame
            || state.primaryMagneticInSolveFrame) {
            throw OpalException(
                    "Pic2d5FramePolicy::enter", "The FFT2D5 frame state is already active.");
        }
        const auto& trackerToSolve =
                context.stepState().frames.trackerToSolve->native<CoordinateSystemTrafo>();
        state.primaryPositionsInSolveFrame = true;
        trackerToSolve.transformBunchTo(primary.R.getView(), primary.getLocalNum());
    }

    void Pic2d5FramePolicy::markComputedFields(State& state) const noexcept {
        state.primaryElectricInSolveFrame = true;
        state.primaryMagneticInSolveFrame = true;
    }

    void Pic2d5FramePolicy::leave(
            const SolveContext& context, ParticleContainer& primary, State& state) const {
        validate(context);
        const auto& solveToTracker =
                context.stepState().frames.solveToTracker->native<CoordinateSystemTrafo>();
        if (state.primaryPositionsInSolveFrame) {
            solveToTracker.transformBunchTo(primary.R.getView(), primary.getLocalNum());
            state.primaryPositionsInSolveFrame = false;
        }
        if (state.primaryElectricInSolveFrame) {
            solveToTracker.rotateBunchTo(primary.E.getView(), primary.getLocalNum());
            state.primaryElectricInSolveFrame = false;
        }
        if (state.primaryMagneticInSolveFrame) {
            solveToTracker.rotateBunchTo(primary.B.getView(), primary.getLocalNum());
            state.primaryMagneticInSolveFrame = false;
        }
    }

    void Pic2d5FramePolicy::restore(
            const SolveContext& context, ParticleContainer& primary, State& state) const noexcept {
        try {
            leave(context, primary, state);
        } catch (...) {
        }
    }

}  // namespace opalx::spacecharge
