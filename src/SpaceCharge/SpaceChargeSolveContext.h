/**
 * @file SpaceChargeSolveContext.h
 * @brief Defines all borrowed and per-call state needed for one space-charge solve.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H
#define OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H

#include "PartBunch/ParticleFrame.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace opalx::spacecharge {

    /**
     * @brief Frame an algorithm solves in, relative to the primary container's reference frame.
     *
     * A rigid (Galilean) frame, not a Lorentz boost. Axes and origin follow the primary
     * container.
     */
    enum class SolveFrame {
        Reference,  ///< The reference frame itself: no transform.
        Centroid,   ///< Origin at the centroid, z along the mean momentum.
        PathLength  ///< z along the mean momentum, z = path length at the reference particle.
    };

    /**
     * @brief Build the requested solve frame from the primary container's moments.
     *
     * @param meanR, meanP Centroid [m] and mean momentum (beta-gamma) in the reference frame.
     * @param fallbackP Direction used when the mean momentum vanishes; +z if it vanishes too.
     * @param pathOrigin z of the reference particle in a PathLength frame [m].
     */
    ParticleFrame makeSolveFrame(
            SolveFrame request, const ippl::Vector<double, 3>& meanR, ippl::Vector<double, 3> meanP,
            const ippl::Vector<double, 3>& fallbackP, double pathOrigin);

    /** @brief Tracker state captured for one space-charge solve. */
    struct SpaceChargeStepState {
        std::size_t step       = 0;    ///< Global tracker step.
        double time            = 0.0;  ///< Tracker time in seconds.
        double timeStep        = 0.0;  ///< Tracker time step in seconds.
        bool emissionActive    = false;
        double emittedFraction = 1.0;  ///< Least-complete active source fraction in [0, 1].
        int mpiSize            = 1;    ///< Active communicator size.
        /// Frame solve() enters; fixed bounds and planeZ are read in it. The tracker passes
        /// makeSolveFrame() of the algorithm's solveFrame() request; unit tests may pass any
        /// rigid frame.
        ParticleFrame solveFrame;
    };

    /**
     * @brief Per-container activity and tracker state for one space-charge call.
     *
     * @c trackingActive[i] corresponds to the ith particle container supplied at construction;
     * zero means inactive. The activity span is borrowed for the call.
     */
    class SpaceChargeSolveContext {
    public:
        SpaceChargeSolveContext(
                std::span<const std::uint8_t> trackingActive, SpaceChargeStepState stepState);

        [[nodiscard]] std::span<const std::uint8_t> trackingActive() const {
            return trackingActive_m;
        }
        [[nodiscard]] const SpaceChargeStepState& stepState() const { return stepState_m; }

    private:
        std::span<const std::uint8_t> trackingActive_m;
        SpaceChargeStepState stepState_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H
