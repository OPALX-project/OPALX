// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPALX_BORIS_STEP_CONTROL_H
#define OPALX_BORIS_STEP_CONTROL_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

#include "Algorithms/CompensatedSum.h"

namespace boris_step {
    /**
     * @brief Time [s] required to move 64 coordinate roundoff scales.
     *
     * positionScale is the largest coordinate/translation scale [m] used by
     * tracker and solver frames, clamped to at least one metre; speed is the
     * positive physical speed [m/s]. The controlled ring PIC frame has zero
     * origin translation, so its scale is set by laboratory ring coordinates. Using
     * the particle/reference speed rather than c prevents low-beta boundary
     * substeps whose movement is lost in repeated coordinate transformations.
     */
    inline double positionTimeFloor(double positionScale, double speed) {
        if (!std::isfinite(positionScale) || positionScale < 0 || !std::isfinite(speed)
            || speed <= 0)
            throw std::invalid_argument("Invalid coordinate scale or speed for Boris time floor");
        const double floor =
                (64 * std::numeric_limits<double>::epsilon() * std::max(1., positionScale)) / speed;
        if (!std::isfinite(floor))
            throw std::overflow_error("Boris position time floor is not representable");
        return floor;
    }

    /**
     * @brief Internal, host-side schedule of dyadic Boris intervals [s].
     *
     * Every accepted interval retains both children of each rejected interval;
     * completion depends on the pending stack, never on a rounded remaining-time
     * subtraction. The clock is relative to this nominal step, avoiding a floor
     * that grows with the absolute simulation time. The refinement floor is
     * max(1e-12 * nominalDt, minimumDt), with an additional depth limit of 64.
     *
     * This class does not operate on particles or call MPI. The tracker must make
     * one collective accept/split decision and apply it identically on every rank
     * so that collective space-charge solves see a common physical midpoint.
     * Trial side effects must be discarded by the caller before split(). A work
     * budget bounds attempted intervals, including rejected parents.
     */
    class Control {
    public:
        explicit Control(
                double nominalDt, double minimumDt = 0, std::size_t maximumTrials = 1000000)
            : nominalDt_m(nominalDt), maximumTrials_m(maximumTrials) {
            if (!std::isfinite(nominalDt) || nominalDt <= 0 || !std::isfinite(minimumDt)
                || minimumDt < 0 || maximumTrials == 0)
                throw std::invalid_argument("Invalid Boris substep schedule");
            timeFloor_m  = std::max(1e-12 * nominalDt, minimumDt);
            pending_m[0] = 0;
        }

        bool done() const { return count_m == 0; }

        /// Duration [s] of the current trial interval; requires !done().
        double step() const {
            requirePending();
            return std::ldexp(nominalDt_m, -static_cast<int>(pending_m[count_m - 1]));
        }

        /// Sum [s] of accepted intervals, relative to the nominal step start.
        double elapsed() const { return sum_m - correction_m; }

        /// Whether a further bisection is numerically resolvable (budget excluded).
        bool canSplit() const {
            if (done() || pending_m[count_m - 1] >= maximumDepth) return false;
            const double interval = step(), half = std::ldexp(interval, -1);
            // Odd subnormal intervals cannot be bisected without changing their
            // duration. Retain that interval instead of rounding both children.
            return interval > timeFloor_m && half > 0 && half + half == interval;
        }

        /** Reject the current trial and visit its first half next.
         * Throws before changing state if refinement or the work budget is exhausted.
         */
        void split() {
            requirePending();
            if (!canSplit())
                throw std::runtime_error("Boris boundary refinement reached its time floor");
            requireNextTrial();
            const auto depth       = static_cast<unsigned char>(pending_m[count_m - 1] + 1);
            pending_m[count_m - 1] = depth;  // Retain the second half.
            pending_m[count_m++]   = depth;  // Visit the first half.
            ++trials_m;
        }

        /** Accept the current interval and visit the next pending interval.
         * When no intervals remain, elapsed() is the nominal duration exactly.
         */
        void accept() {
            requirePending();
            if (count_m > 1) requireNextTrial();
            compensated::add(step(), sum_m, correction_m);
            --count_m;
            if (done()) {
                sum_m        = nominalDt_m;
                correction_m = 0;
            } else {
                ++trials_m;
            }
        }

    private:
        void requirePending() const {
            if (done()) throw std::logic_error("Boris substep schedule is already complete");
        }
        void requireNextTrial() const {
            if (trials_m >= maximumTrials_m)
                throw std::runtime_error("Boris boundary refinement exceeded its trial budget");
        }

        static constexpr unsigned maximumDepth = 64;
        std::array<unsigned char, maximumDepth + 1> pending_m{};
        std::size_t count_m = 1, trials_m = 1;
        double nominalDt_m, timeFloor_m = 0, sum_m = 0, correction_m = 0;
        std::size_t maximumTrials_m;
    };
}  // namespace boris_step

#endif
