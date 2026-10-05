#ifndef OPALX_DIRECTED_TURN_COUNTER_H
#define OPALX_DIRECTED_TURN_COUNTER_H
#include <cmath>
#include <stdexcept>
#include "OPALTypes.h"
#include "Physics/Physics.h"

/** @brief Counts directed reference-trajectory returns through a fixed launch plane.
 * The fixed plane normal is the initial momentum direction. A negative-side
 * excursion greater than tolerance metres arms the next positive crossing.
 * The initial point and the opposite-direction crossing cannot count as turns.
 * update() counts supplied endpoints. terminalStep() can locate the final requested
 * crossing before the caller commits a step; earlier turns remain endpoint-counted.
 *
 * All input states must use the same Cartesian frame. No geometry circumference
 * or assumed orbit radius enters the count. This host-side object is maintained
 * independently for each particle-container reference, not for each macro-particle.
 * Its origin, normal, armed flag, and count are not currently checkpointed.
 * @see ParallelTracker::setRequestedTurns
 */
class DirectedTurnCounter {
public:
    /** @param origin Fixed launch position [m].
     * @param momentum Initial mechanical momentum [beta*gamma]; normalized internally.
     * @param tolerance Positive hysteresis distance [m], default 1 nm.
     * @pre Origin and tolerance are finite, and tolerance is positive.
     * @throws std::invalid_argument If momentum has zero or nonfinite norm.
     */
    DirectedTurnCounter(
            const Vector_t<double, 3>& origin, const Vector_t<double, 3>& momentum,
            double tolerance = 1e-9)
        : origin(origin), tolerance(tolerance) {
        double norm = std::sqrt(dot(momentum, momentum));
        if (!(norm > 0) || !std::isfinite(norm))
            throw std::invalid_argument("Invalid turn-counter momentum");
        normal = momentum / norm;
    }
    /** @brief Observe the next endpoint in forward tracking order.
     * @param position Current reference position [m].
     * @param momentum Current reference momentum [beta*gamma].
     * @return True only when the armed count increments: the point is on the
     * nonnegative side and its momentum has positive projection on the plane normal.
     * @pre Sample successive timesteps finely enough to resolve each negative-side
     * excursion. This method does not interpolate the crossing or shorten a step.
     */
    bool update(const Vector_t<double, 3>& position, const Vector_t<double, 3>& momentum) {
        const Vector_t<double, 3> displacement = position - origin;
        const double distance                  = dot(displacement, normal);
        if (distance < -tolerance) armed = true;
        if (armed && distance >= 0 && dot(momentum, normal) > 0) {
            ++turns;
            armed = false;
            return true;
        }
        return false;
    }
    /// Number of accepted crossings since construction (initially zero).
    unsigned long long count() const { return turns; }

    /** Locate only the final requested forward return by reintegration from the
     * unchanged start. advance(h) returns position, momentum and hitMaterial.
     * No counter state changes here: update() observes the subsequently committed
     * endpoint. Static magnetic tracking and a crossing-resolving DT are assumed.
     *
     * A light-speed bound avoids field trials away from the plane. Bisection keeps
     * the nonnegative endpoint, with relative time resolution 1e-12 of nominal DT
     * (as in OneTurnMap) and a 1 nm section acceptance equal to the arming distance.
     * A discontinuous/unresolved event fails rather than projecting the position.
     */
    template <class Advance>
    double terminalStep(
            unsigned long long requested, const Vector_t<double, 3>& position, double dt,
            Advance&& advance) const {
        if (!(dt > 0) || !std::isfinite(dt))
            throw std::invalid_argument("Terminal turn requires finite positive DT");
        const Vector_t<double, 3> displacement = position - origin;
        const double before                    = dot(displacement, normal);
        if (!armed || requested == 0 || turns != requested - 1 || before >= 0
            || -before > Physics::c * dt)
            return dt;
        auto end            = advance(dt);
        const auto distance = [&](const auto& state) {
            const Vector_t<double, 3> displacement = state.position - origin;
            const double d                         = dot(displacement, normal);
            for (unsigned i = 0; i < 3; ++i)
                if (!std::isfinite(state.position(i)) || !std::isfinite(state.momentum(i)))
                    throw std::runtime_error("Nonfinite terminal-turn reference trial");
            return d;
        };
        if (distance(end) < 0 || dot(end.momentum, normal) <= 0) return dt;
        double lo = 0, hi = dt;
        for (unsigned iteration = 0; iteration < 80; ++iteration) {
            const double mid = lo + (hi - lo) / 2;
            if (mid == lo || mid == hi || hi - lo <= 1e-12 * dt) break;
            auto trial = advance(mid);
            if (distance(trial) >= 0) {
                hi  = mid;
                end = trial;
            } else
                lo = mid;
        }
        if (distance(end) > tolerance || dot(end.momentum, normal) <= 0 || end.hitMaterial)
            throw std::runtime_error(
                    "Terminal turn is unresolved or intercepts material; reduce DT");
        return hi;
    }

private:
    Vector_t<double, 3> origin, normal;
    double tolerance;
    bool armed               = false;
    unsigned long long turns = 0;
};
#endif
