// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPALX_PASSIVE_PROBE_H
#define OPALX_PASSIVE_PROBE_H

#include <Kokkos_Core.hpp>
#include <cstdint>
#include <limits>
#include "OPALTypes.h"

/** @brief Passive, directed plane observations of accepted particle endpoints.
 *
 * All endpoints and planes use one fixed Cartesian frame. Position and the
 * arming tolerance are in metres, time in seconds, and mechanical momentum is
 * normalized by mc. Endpoints must represent the same physical time for both
 * position and momentum, rather than staggered leapfrog quantities.
 *
 * Linear interpolation is exact for constant drift and has O(dt^2) event error
 * for a smooth, transverse crossing. It cannot repair transport error or detect
 * multiple crossings hidden inside one accepted interval. Place observation
 * planes away from discontinuous fields and check timestep convergence.
 *
 * Only diagnostic state changes. Neither particles nor their step schedule are
 * writable through this interface. No allocation, reduction or exception occurs
 * on the device. Stable particle IDs and migration are the caller's concern.
 */
namespace passive_probe {
    using Vector = Vector_t<double, 3>;

    struct Endpoint {
        Vector position = Vector(0.0), momentum = Vector(0.0);
        double time = 0.0;
    };

    struct Plane {
        Vector origin = Vector(0.0), normal = Vector(0.0, 0.0, 1.0);
        double tolerance = 1e-9;
    };

    struct State {
        Endpoint previous;
        std::uint64_t turns = 0;
        bool initialized = false, armed = false;
    };

    struct Sample {
        Endpoint crossing;
        std::uint64_t turn = 0;
        double fraction    = 0.0;
    };

    enum class Status {
        Disabled,
        Initialized,
        Advanced,
        Crossing,
        InvalidInput,
        DecreasingTime,
        TurnOverflow
    };

    namespace detail {
        KOKKOS_INLINE_FUNCTION bool finite(const Endpoint& endpoint) {
            if (!Kokkos::isfinite(endpoint.time)) return false;
            for (unsigned d = 0; d < 3; ++d)
                if (!Kokkos::isfinite(endpoint.position(d))
                    || !Kokkos::isfinite(endpoint.momentum(d)))
                    return false;
            return true;
        }
    }  // namespace detail

    /** @brief Observe one accepted endpoint, optionally recording a forward crossing.
     *
     * The first call initializes history without recording the initial plane point.
     * A visit farther than tolerance onto the negative side arms the next return.
     * A crossing requires a negative-to-nonnegative position bracket and positive
     * interpolated momentum along the plane normal. The normal need not have unit
     * length. Exact endpoints on the plane are recorded once. Opposite crossings
     * and near-plane jitter cannot create duplicate turns.
     *
     * Disabled calls and errors leave state and sample unchanged. Successful calls
     * always advance state; sample changes only for Status::Crossing. Successive
     * accepted positive integration intervals can be smaller than one ULP of the
     * absolute clock, giving identical represented endpoint times. Such endpoints
     * still update position/momentum history. A bracketed crossing remains defined
     * by its spatial fraction and is recorded at that common rounded timestamp;
     * no division by the time interval is needed. Strictly decreasing time,
     * invalid/overflowing arithmetic, and turn-count overflow are reported, never
     * thrown. Missing-particle intervals must be handled by resetting history
     * externally, rather than interpolating across an unobserved interval.
     */
    KOKKOS_INLINE_FUNCTION Status
    update(const Plane& plane, State& state, const Endpoint& accepted, Sample& sample,
           bool enabled = true) {
        if (!enabled) return Status::Disabled;
        if (!detail::finite(accepted) || !Kokkos::isfinite(plane.tolerance) || plane.tolerance <= 0)
            return Status::InvalidInput;
        for (unsigned d = 0; d < 3; ++d)
            if (!Kokkos::isfinite(plane.origin(d)) || !Kokkos::isfinite(plane.normal(d)))
                return Status::InvalidInput;
        const double norm2 = dot(plane.normal, plane.normal);
        if (!Kokkos::isfinite(norm2) || norm2 <= 0) return Status::InvalidInput;
        const Vector normal       = plane.normal / Kokkos::sqrt(norm2);
        const Vector displacement = accepted.position - plane.origin;
        const double after        = dot(displacement, normal);
        if (!Kokkos::isfinite(after)) return Status::InvalidInput;
        if (!state.initialized) {
            state.previous    = accepted;
            state.initialized = true;
            state.armed       = after < -plane.tolerance;
            return Status::Initialized;
        }
        if (!detail::finite(state.previous)) return Status::InvalidInput;
        const double dt = accepted.time - state.previous.time;
        if (!Kokkos::isfinite(dt)) return Status::InvalidInput;
        if (dt < 0) return Status::DecreasingTime;
        const Vector previousDisplacement = state.previous.position - plane.origin;
        const double before               = dot(previousDisplacement, normal);
        if (!Kokkos::isfinite(before)) return Status::InvalidInput;

        if (state.armed && before < 0 && after >= 0) {
            const double denominator = after - before;
            if (!Kokkos::isfinite(denominator) || denominator <= 0) return Status::InvalidInput;
            Sample candidate;
            candidate.fraction      = -before / denominator;
            candidate.crossing.time = state.previous.time + candidate.fraction * dt;
            for (unsigned d = 0; d < 3; ++d) {
                candidate.crossing.position(d) =
                        state.previous.position(d)
                        + candidate.fraction * (accepted.position(d) - state.previous.position(d));
                candidate.crossing.momentum(d) =
                        state.previous.momentum(d)
                        + candidate.fraction * (accepted.momentum(d) - state.previous.momentum(d));
            }
            if (!Kokkos::isfinite(candidate.fraction) || candidate.fraction < 0
                || candidate.fraction > 1 || !detail::finite(candidate.crossing))
                return Status::InvalidInput;
            const double direction = dot(candidate.crossing.momentum, normal);
            if (!Kokkos::isfinite(direction)) return Status::InvalidInput;
            if (direction > 0) {
                if (state.turns == std::numeric_limits<std::uint64_t>::max())
                    return Status::TurnOverflow;
                candidate.turn = state.turns + 1;
                sample         = candidate;
                state.turns    = candidate.turn;
                state.previous = accepted;
                state.armed    = false;
                return Status::Crossing;
            }
        }
        state.previous = accepted;
        if (after < -plane.tolerance)
            state.armed = true;
        else if (after >= 0)
            state.armed = false;
        return Status::Advanced;
    }
}  // namespace passive_probe

#endif
