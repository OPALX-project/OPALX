#ifndef OPALX_DIRECTED_TURN_COUNTER_H
#define OPALX_DIRECTED_TURN_COUNTER_H
#include <cmath>
#include <stdexcept>
#include "OPALTypes.h"

/** Counts returns of a reference trajectory through its launch plane.
 * The fixed plane normal is the initial momentum direction. A negative-side
 * excursion of at least tolerance metres arms the next positive crossing.
 * The initial point and the opposite-direction crossing cannot count as turns.
 * Completion is reported at the first integration endpoint beyond the plane.
 */
class DirectedTurnCounter {
public:
    DirectedTurnCounter(
            const Vector_t<double, 3>& origin, const Vector_t<double, 3>& momentum,
            double tolerance = 1e-9)
        : origin(origin), tolerance(tolerance) {
        double norm = std::sqrt(dot(momentum, momentum));
        if (!(norm > 0) || !std::isfinite(norm))
            throw std::invalid_argument("Invalid turn-counter momentum");
        normal = momentum / norm;
    }
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
    unsigned long long count() const { return turns; }

private:
    Vector_t<double, 3> origin, normal;
    double tolerance;
    bool armed               = false;
    unsigned long long turns = 0;
};
#endif
