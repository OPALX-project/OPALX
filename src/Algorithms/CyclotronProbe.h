// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPALX_CYCLOTRON_PROBE_H
#define OPALX_CYCLOTRON_PROBE_H

#include <climits>
#include "Algorithms/PassiveProbe.h"

/** @brief Cyclotron segment crossings in the fixed laboratory X-Z plane.
 *
 * Legacy (XSTART, YSTART) denotes laboratory (x,z), in metres after parsing;
 * laboratory y is vertical and is unrestricted. Both crossing directions count.
 * Two directed observers suppress initial-plane hits and endpoint duplicates,
 * with a 1 nm arming distance. A clipped hit increments the particle's passage
 * number (one-based), independently of reference-turn or local array indices.
 * Linear interpolation has O(dt^2) event error for smooth transverse motion;
 * it cannot resolve hidden intra-step recrossings or RF discontinuities.
 */
namespace cyclotron_probe {
    using Vector   = passive_probe::Vector;
    using Endpoint = passive_probe::Endpoint;
    using Sample   = passive_probe::Sample;
    using Status   = passive_probe::Status;

    struct Segment {
        Vector start = Vector(0.0), end = Vector(0.0);
    };

    /// Trivially copyable particle attribute: migrated/sorted with its owner.
    struct State {
        passive_probe::State forward, reverse;
        std::int64_t id = -1;
        int passages    = 0;
    };

    KOKKOS_INLINE_FUNCTION Status
    update(const Segment& segment, State& state, std::int64_t id, const Endpoint& endpoint,
           Sample& sample) {
        const Vector tangent = segment.end - segment.start;
        const double length2 = dot(tangent, tangent);
        if (id < 0 || !Kokkos::isfinite(length2) || length2 <= 0 || segment.start(1) != 0
            || segment.end(1) != 0)
            return Status::InvalidInput;
        passive_probe::Plane plane;
        plane.origin    = segment.start;
        plane.normal    = Vector(-tangent(2), 0.0, tangent(0));
        State candidate = state.id == id ? state : State{};
        candidate.id    = id;
        Sample forward, reverse;
        const auto a = passive_probe::update(plane, candidate.forward, endpoint, forward);
        plane.normal = -plane.normal;
        const auto b = passive_probe::update(plane, candidate.reverse, endpoint, reverse);
        if (a == Status::InvalidInput || a == Status::DecreasingTime || a == Status::TurnOverflow)
            return a;
        if (b == Status::InvalidInput || b == Status::DecreasingTime || b == Status::TurnOverflow)
            return b;
        if (a == Status::Crossing || b == Status::Crossing) {
            Sample hit            = a == Status::Crossing ? forward : reverse;
            const Vector offset   = hit.crossing.position - segment.start;
            const double fraction = dot(offset, tangent) / length2;
            // Include segment endpoints; tolerate only floating-point projection rounding.
            constexpr double eps = 16 * Kokkos::epsilon_v<double>;
            if (fraction >= -eps && fraction <= 1 + eps) {
                if (candidate.passages == INT_MAX) return Status::TurnOverflow;
                hit.turn = ++candidate.passages;
                state    = candidate;
                sample   = hit;
                return Status::Crossing;
            }
        }
        state = candidate;
        return a == Status::Initialized ? a : Status::Advanced;
    }
}  // namespace cyclotron_probe
#endif
