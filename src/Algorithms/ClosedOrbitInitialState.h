#ifndef OPAL_CLOSED_ORBIT_INITIAL_STATE_H
#define OPAL_CLOSED_ORBIT_INITIAL_STATE_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/PartData.h"
#include "Utilities/OpalException.h"

/** Immutable-by-value COF launch snapshot, in laboratory metres and p/(mc).
 * section maps laboratory coordinates to the COF Poincare section. The bunch
 * frame is the shortest rotation from section +z to the solved momentum; local
 * particle momenta include their full longitudinal momentum, not a delta-p.
 * This rigid placement changes neither momentum norms nor internal distances.
 * It supplies an external-field orbit, not a space-charge matched distribution.
 */
struct ClosedOrbitInitialState {
    ippl::Vector<double, 3> position{0.0, 0.0, 0.0}, momentum{0.0, 0.0, 0.0};
    CoordinateSystemTrafo section;
    double time = 0, massEV = 0, chargeE = 0;
    std::string line, species, json;

    /// Validate species, charge, mass and reference momentum without changing BEAM.
    void validate(
            const std::string& ring, const std::string& particle, const PartData& reference) const {
        const auto equal = [](double a, double b) {
            return std::isfinite(a) && std::isfinite(b)
                   && std::abs(a - b) <= 64 * std::numeric_limits<double>::epsilon()
                                                 * std::max(std::abs(a), std::abs(b));
        };
        if (ring != line || particle != species || !equal(massEV, reference.getM())
            || !equal(chargeE, reference.getQ())
            || !equal(
                    std::hypot(momentum[0], momentum[1], momentum[2]),
                    reference.getP() / reference.getM()))
            throw OpalException(
                    "INITIALORBIT", "COF lattice, species or energy does not match TRACK.");
    }

    /// Laboratory-to-orbit rigid frame; forward momentum fixes its transverse roll.
    CoordinateSystemTrafo frame() const {
        const ippl::Vector<double, 3> p = section.rotateTo(momentum);
        const double norm               = std::hypot(p[0], p[1], p[2]);
        if (!std::isfinite(norm) || norm <= 0 || p[2] <= 0 || !std::isfinite(time))
            throw OpalException("INITIALORBIT", "Invalid forward COF launch state.");
        for (unsigned d = 0; d < 3; ++d)
            if (!std::isfinite(position[d]))
                throw OpalException("INITIALORBIT", "Nonfinite COF launch position.");
        Quaternion align(1 + p[2] / norm, -p[1] / norm, p[0] / norm, 0);
        align.normalize();
        return CoordinateSystemTrafo(position, align.conjugate() * section.getRotation());
    }
};
#endif
