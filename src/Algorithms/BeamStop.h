//
// Header BeamStop
//   When a beam stops being tracked: its ZSTOP plane and its UNTILMAXSTEPS flag.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_BEAM_STOP_H
#define OPALX_BEAM_STOP_H

#include "OPALTypes.h"

#include <Kokkos_Core.hpp>

#include <cmath>
#include <limits>
#include <optional>
#include <vector>

/**
 * @brief Stopping state of one beam (one particle container) on a LINE.
 *
 * A beam with a ZSTOP is tracked until its last particle has passed the plane through the
 * reference particle at path length ZSTOP, normal to the reference momentum there. The primary
 * beam takes its ZSTOP from the last TRACK ZSTOP entry, every other beam from its BEAM ZSTOP.
 * A beam with UNTILMAXSTEPS keeps the run going until MAXSTEPS is used up. A beam with neither
 * is tracked while the run lasts.
 */
struct BeamStop {
    std::optional<double> zStop;  ///< Path length [m] of the stop plane; empty without ZSTOP.
    bool untilMaxSteps = false;   ///< Keep the run going until MAXSTEPS for this beam.
    bool planeSet      = false;   ///< The reference reached zStop; the plane below is set.
    Vector_t<double, 3> planeOrigin = Vector_t<double, 3>(0.0);  ///< Lab position [m].
    Vector_t<double, 3> planeNormal = Vector_t<double, 3>(0.0);  ///< Lab unit vector.
    bool done = false;  ///< Every particle passed the plane; the beam is no longer tracked.
};

namespace beam_stop {

    /// Whether a reference at path length @p sPos has reached @p zStop, to half a step @p drift.
    inline bool reached(const double sPos, const double zStop, const double drift) {
        return sPos > zStop - 0.5 * drift;
    }

    /**
     * @brief Put the stop plane at path length zStop.
     *
     * The reference is within half a step of zStop; the plane origin is moved along the
     * reference momentum by the remaining path length.
     * @param stop Beam whose zStop is set.
     * @param position Reference position in the lab [m].
     * @param momentum Reference momentum in the lab (beta gamma), nonzero.
     * @param sPos Reference path length [m].
     */
    inline void setPlane(
            BeamStop& stop, const Vector_t<double, 3>& position,
            const Vector_t<double, 3>& momentum, const double sPos) {
        const Vector_t<double, 3> normal = momentum / euclidean_norm(momentum);
        stop.planeNormal                 = normal;
        stop.planeOrigin                 = position + (*stop.zStop - sPos) * normal;
        stop.planeSet                    = true;
    }

    /**
     * @brief Smallest signed distance (r - origin) . normal of the first @p count positions.
     * @return Infinity when @p count is zero, so an empty beam counts as passed.
     */
    template <class View>
    double minDistanceToPlane(
            const View& positions, const std::size_t count, const Vector_t<double, 3>& origin,
            const Vector_t<double, 3>& normal) {
        double minimum = std::numeric_limits<double>::infinity();
        Kokkos::parallel_reduce(
                "beam_stop::minDistanceToPlane", count,
                KOKKOS_LAMBDA(const std::size_t i, double& local) {
                    const Vector_t<double, 3> d = positions(i) - origin;
                    const double distance = d(0) * normal(0) + d(1) * normal(1) + d(2) * normal(2);
                    local                 = distance < local ? distance : local;
                },
                Kokkos::Min<double>(minimum));
        return minimum;
    }

    /**
     * @brief Whether the run ends before MAXSTEPS is used up.
     *
     * True when at least one beam has a ZSTOP, every beam with a ZSTOP is done, and no beam has
     * UNTILMAXSTEPS. Without any ZSTOP the run ends by MAXSTEPS.
     */
    inline bool runFinished(const std::vector<BeamStop>& stops) {
        bool anyZStop = false;
        for (const BeamStop& stop : stops) {
            if (stop.untilMaxSteps) return false;
            if (stop.zStop) {
                anyZStop = true;
                if (!stop.done) return false;
            }
        }
        return anyZStop;
    }

}  // namespace beam_stop

#endif  // OPALX_BEAM_STOP_H
