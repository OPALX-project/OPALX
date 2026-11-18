//
// Header ElementRegions
//   Device-callable tests of where an element applies its field and where it deletes particles.
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
#ifndef OPALX_ELEMENT_REGIONS_H
#define OPALX_ELEMENT_REGIONS_H

#include "BeamlineGeometry/Geometry.h"
#include "OPALTypes.h"

#include <Kokkos_Core.hpp>

#include <limits>
#include <utility>
#include <vector>

enum class ApertureType : unsigned short { RECTANGULAR, ELLIPTICAL };

/**
 * @brief Stateless, device-callable aperture functions.
 */
namespace ApertureHelper {

    /**
     * @brief Transverse-aperture containment test.
     *
     * Limits are half-apertures; the boundary itself counts as outside. The
     * aperture is constant along the element.
     *
     * @param x transverse horizontal position [m]
     * @param y transverse vertical position [m]
     * @param type aperture shape
     * @param xLimit horizontal half-aperture [m]
     * @param yLimit vertical half-aperture [m]
     */
    KOKKOS_INLINE_FUNCTION bool isInsideAperture(
            const double x, const double y, const ApertureType type, const double xLimit,
            const double yLimit) {
        switch (type) {
            case ApertureType::RECTANGULAR:
                return (Kokkos::abs(x) < xLimit && Kokkos::abs(y) < yLimit);
            case ApertureType::ELLIPTICAL: {
                const double u = x / xLimit;
                const double v = y / yLimit;
                return (u * u + v * v < 1.0);
            }
            default:
                return false;
        }
    }

    /**
     * @brief Host convenience overload
     *
     * @param r position in the element's coordinatesystem
     * @param aperture vector of aperture coefficients (usually aperture_m)
     */
    inline bool isInsideAperture(
            const Vector_t<double, 3>& r,
            const std::pair<ApertureType, std::vector<double>>& aperture) {
        return isInsideAperture(r[0], r[1], aperture.first, aperture.second[0], aperture.second[1]);
    }

}  // namespace ApertureHelper

/// A rectangle or ellipse in the transverse plane, centred at (centreX, centreY) [m].
struct ApertureShape {
    ApertureType type = ApertureType::RECTANGULAR;
    double centreX    = 0.0;
    double centreY    = 0.0;
    double halfX      = std::numeric_limits<double>::max();
    double halfY      = std::numeric_limits<double>::max();

    /// Whether (x, y) lies inside; the boundary counts as outside.
    KOKKOS_INLINE_FUNCTION bool contains(const double x, const double y) const {
        return ApertureHelper::isInsideAperture(x - centreX, y - centreY, type, halfX, halfY);
    }

    /// Whether the shape lies inside |x| <= limitX, |y| <= limitY.
    bool fitsInside(const double limitX, const double limitY) const {
        return Kokkos::abs(centreX) + halfX <= limitX && Kokkos::abs(centreY) + halfY <= limitY;
    }
};

/**
 * @brief Where an element applies its field and where it deletes particles.
 *
 * All positions are in the element's local frame. A sector bend (curvature != 0) measures them
 * in the arc coordinates of GeometryHelper::toBendArcCoords(): x radial, y vertical, z the arc
 * length. For a particle at (x, y, z) in those coordinates:
 * - the field applies if fieldBegin <= z < fieldEnd and (x, y) is inside the aperture;
 * - the particle is deleted if deletesParticles is set, bodyBegin <= z < bodyEnd,
 *   |x| < reachX, |y| < reachY and (x, y) is outside the aperture;
 * - otherwise the element does nothing to it.
 *
 * A default-constructed object applies the field everywhere and deletes nothing.
 */
struct ElementRegions {
    double fieldBegin = -std::numeric_limits<double>::max();  ///< Field window start [m].
    double fieldEnd   = std::numeric_limits<double>::max();   ///< Field window end [m].
    double bodyBegin  = 0.0;                                  ///< Body start [m].
    double bodyEnd    = 0.0;                                  ///< Body end [m].
    ApertureShape aperture;                                   ///< Inside: field; outside: material.
    double reachX     = 0.0;  ///< Horizontal half-width of the element [m].
    double reachY     = 0.0;  ///< Vertical half-width of the element [m].
    double curvature  = 0.0;  ///< Nonzero for a sector bend [1/m].
    double bodyLength = 0.0;  ///< Arc length of a sector bend [m].
    /// False when nothing around the aperture is material, e.g. a field map without APERTURE.
    bool deletesParticles = false;

    /// The local position in the coordinates the windows are measured in.
    KOKKOS_INLINE_FUNCTION Vector_t<double, 3> toArcCoords(const Vector_t<double, 3>& r) const {
        return curvature != 0.0 ? GeometryHelper::toBendArcCoords(r, curvature, bodyLength) : r;
    }

    /// Field test for a position already in toArcCoords() coordinates.
    KOKKOS_INLINE_FUNCTION bool fieldAppliesAt(const Vector_t<double, 3>& p) const {
        return p(2) >= fieldBegin && p(2) < fieldEnd && aperture.contains(p(0), p(1));
    }

    /// Loss test for a position already in toArcCoords() coordinates.
    KOKKOS_INLINE_FUNCTION bool particleLostAt(const Vector_t<double, 3>& p) const {
        return deletesParticles && p(2) >= bodyBegin && p(2) < bodyEnd && Kokkos::abs(p(0)) < reachX
               && Kokkos::abs(p(1)) < reachY && !aperture.contains(p(0), p(1));
    }

    /// Field test for a local Cartesian position.
    KOKKOS_INLINE_FUNCTION bool inFieldRegion(const Vector_t<double, 3>& r) const {
        return fieldAppliesAt(toArcCoords(r));
    }

    /// Loss test for a local Cartesian position.
    KOKKOS_INLINE_FUNCTION bool inLossRegion(const Vector_t<double, 3>& r) const {
        return particleLostAt(toArcCoords(r));
    }
};

#endif  // OPALX_ELEMENT_REGIONS_H
