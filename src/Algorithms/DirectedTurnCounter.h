#ifndef OPALX_DIRECTED_TURN_COUNTER_H
#define OPALX_DIRECTED_TURN_COUNTER_H
#include <cmath>
#include <stdexcept>
#include "OPALTypes.h"

/** @brief Counts directed reference-trajectory returns through a fixed launch plane.
 * The fixed plane normal is the initial momentum direction. A negative-side
 * excursion greater than tolerance metres arms the next positive crossing.
 * The initial point and the opposite-direction crossing cannot count as turns.
 * Completion is reported at the first integration endpoint beyond the plane.
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

private:
    Vector_t<double, 3> origin, normal;
    double tolerance;
    bool armed               = false;
    unsigned long long turns = 0;
};
#endif
