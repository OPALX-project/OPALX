/**
 * @file SpaceChargeSolveContext.cpp
 * @brief Implements construction of one-call space-charge contexts.
 */

#include "SpaceCharge/SpaceChargeSolveContext.h"

#include "Algorithms/Quaternion.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace opalx::spacecharge {

    ParticleFrame makeSolveFrame(
            SolveFrame request, const ippl::Vector<double, 3>& meanR, ippl::Vector<double, 3> meanP,
            const ippl::Vector<double, 3>& fallbackP, double pathOrigin) {
        if (request == SolveFrame::Reference) return ParticleFrame();
        // getQuaternion needs a direction; the mean vanishes with zero particles.
        constexpr double tolerance = 1e-12;
        const auto vanishes        = [](const ippl::Vector<double, 3>& p) {
            return p[0] * p[0] + p[1] * p[1] + p[2] * p[2] < tolerance * tolerance;
        };
        if (vanishes(meanP)) meanP = fallbackP;
        if (vanishes(meanP)) meanP = ippl::Vector<double, 3>(0, 0, 1);
        const Quaternion alignment = getQuaternion(meanP, ippl::Vector<double, 3>(0, 0, 1));
        if (request == SolveFrame::Centroid) {
            return ParticleFrame(
                    ParticleFrame::Kind::SpaceChargeSolve, CoordinateSystemTrafo(meanR, alignment));
        }
        const CoordinateSystemTrafo solveToTracker(
                ippl::Vector<double, 3>(0, 0, pathOrigin), alignment.conjugate());
        return ParticleFrame(
                ParticleFrame::Kind::SpaceChargeSolve, solveToTracker.inverted(), solveToTracker);
    }

    SpaceChargeSolveContext::SpaceChargeSolveContext(
            std::span<const std::uint8_t> trackingActive, SpaceChargeStepState stepState)
        : trackingActive_m(trackingActive), stepState_m(std::move(stepState)) {
        if (!std::isfinite(stepState_m.time) || !std::isfinite(stepState_m.timeStep)) {
            throw std::invalid_argument(
                    "SpaceChargeSolveContext requires finite time and time step");
        }
        if (stepState_m.mpiSize < 1) {
            throw std::invalid_argument(
                    "SpaceChargeSolveContext communicator size must be positive");
        }
        if (!(stepState_m.emittedFraction >= 0.0 && stepState_m.emittedFraction <= 1.0)) {
            throw std::invalid_argument(
                    "SpaceChargeSolveContext emitted fraction must be in [0, 1]");
        }
    }

}  // namespace opalx::spacecharge
