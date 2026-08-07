#ifndef OPALX_Geometry_HH
#define OPALX_Geometry_HH

//
// Class Geometry
//   Single geometry class for all beamline elements. Replaces the former
//   geometry class hierarchy (StraightGeometry, PlanarArcGeometry, RBendGeometry,
//   NullGeometry).
//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "Algorithms/CoordinateSystemTrafo.h"
#include "OPALTypes.h"

#include <Kokkos_Core.hpp>

#include <cstddef>
#include <vector>

/// The kind of body geometry represented by a Geometry object.
enum class GeometryKind : unsigned char { Null, Straight, SBend, RBend };

/**
 * @class Geometry
 * @brief Single geometry class for all beamline elements.
 *
 * Carries the body length plus the bend parameters (curvature/angle and
 * entrance/exit pole-face angles) and provides:
 *   - the body length and arc length,
 *   - the bend angle, curvature and chord length,
 *   - the entrance/exit edge transforms as CoordinateSystemTrafo (the single
 *     source for placement), and
 *   - the sampled design path.
 *
 * @note A GeometryKind tag selects straight / planar-arc / rectangular-bend /
 *       null behaviour. For straight and null geometries the edge transforms are
 *       the identity at the entrance and a +z shift to the exit. For bends they
 *       reproduce the planar-arc and rectangular-bend pole-face frames of the
 *       former geometry subclasses, expressed directly as CoordinateSystemTrafo.
 */
class Geometry {
public:
    Geometry() = default;

    /// @name Factories
    ///@{
    static Geometry makeNull();
    static Geometry makeStraight(double length);
    /// Sector-bend (planar-arc) geometry of arc length and curvature h.
    static Geometry makeSBend(double length, double curvature);
    /// Rectangular-bend geometry of straight body length and full bend angle.
    static Geometry makeRBend(double bodyLength, double angle);
    ///@}

    GeometryKind kind() const { return kind_m; }

    /// @name Lengths
    ///@{
    /// Length measured along the design arc.
    double getArcLength() const;
    /// Design / body length (straight body length for rectangular bends).
    double getElementLength() const { return len_m; }
    void setElementLength(double length);
    ///@}

    /// @name Bend parameters
    ///@{
    double getBendAngle() const;
    void setBendAngle(double angle);
    double getCurvature() const { return (kind_m == GeometryKind::SBend) ? h_m : 0.0; }
    void setCurvature(double curvature);
    /// Straight-line distance between the entrance and exit frames.
    double getChordLength() const;
    double getEntranceAngle() const { return entranceAngle_m; }
    void setEntranceAngle(double angle) { entranceAngle_m = angle; }
    double getExitAngle() const { return exitAngle_m; }
    void setExitAngle(double angle) { exitAngle_m = angle; }
    /// Sample the local design path from entrance to exit.
    std::vector<Vector_t<double, 3>> getDesignPath(std::size_t minSamples = 32) const;
    ///@}

    /// @name Edge transforms (single source of truth for placement)
    ///@{
    /// Entrance-frame to entrance-edge transform. The stored frame IS the entrance face for
    /// every kind, so this is always the identity.
    CoordinateSystemTrafo getEdgeToBegin() const;
    /// Entrance-frame to exit-edge transform: a pure +z shift for a straight body (straight
    /// element, rectangular-bend box, null); the arc-end frame turned by the full bend angle
    /// for a sector bend.
    CoordinateSystemTrafo getEdgeToEnd() const;
    ///@}

private:
    double halfAngle() const { return angle_m / 2.0; }
    /// Position (x, 0, z) of the centred body chart at arc-position s.
    Vector_t<double, 3> framePosition(double s) const;

    GeometryKind kind_m    = GeometryKind::Straight;
    double len_m           = 0.0;  ///< design / body length
    double h_m             = 0.0;  ///< curvature (SBend)
    double angle_m         = 0.0;  ///< bend angle (SBend: h*len; RBend: full angle)
    double entranceAngle_m = 0.0;  ///< entrance pole-face angle
    double exitAngle_m     = 0.0;  ///< exit pole-face angle
};

/**
 * @namespace GeometryHelper
 * @brief Stateless, device-callable geometry functions.
 *
 * A Geometry object cannot be used inside a Kokkos lambda, so the coordinate
 * conversions needed by the field kernels live here as pure functions of plain
 * scalars (curvature, body length) that the caller reads from its Geometry on the
 * host and captures by value. The Geometry class stays the single stateful source
 * of those scalars.
 */
namespace GeometryHelper {

    /// |cos(angle)| floored to a small positive value, for pole-face projections.
    KOKKOS_INLINE_FUNCTION double safeAbsCos(const double angle) {
        const double c = Kokkos::abs(Kokkos::cos(angle));
        return (c > 1.0e-6) ? c : 1.0e-6;
    }

    /// Rotate a vector by @p angle about the y-axis (the bend plane is x-z).
    KOKKOS_INLINE_FUNCTION Vector_t<double, 3> rotateAboutY(
            const Vector_t<double, 3>& v, const double angle) {
        const double c = Kokkos::cos(angle);
        const double s = Kokkos::sin(angle);
        return Vector_t<double, 3>(c * v(0) + s * v(2), v(1), -s * v(0) + c * v(2));
    }

    /// Signed design-arc phase of an entrance-frame point on a sector arc of the
    /// given curvature. Kept branch-consistent for negative curvature.
    KOKKOS_INLINE_FUNCTION double referencePhase(
            const Vector_t<double, 3>& entry, const double curvature) {
        const double radius = 1.0 / curvature;
        const double sign   = (radius >= 0.0) ? 1.0 : -1.0;
        return Kokkos::atan2(sign * entry(2), sign * (entry(0) + radius));
    }

    /**
     * @brief Map an entrance-frame Cartesian point to arc coordinates.
     *
     * Returns (radial offset x, vertical y, arc length s measured from the
     * entrance) for the sector design arc of the given curvature and body length —
     * the inverse of the design-arc map anchored at the entrance frame. Upstream
     * of the entrance and downstream of the exit the reference path continues as
     * the straight entrance/exit tangent, so s continues linearly there. Zero
     * curvature (a straight body) returns the point unchanged.
     */
    KOKKOS_INLINE_FUNCTION Vector_t<double, 3> toBendArcCoords(
            const Vector_t<double, 3>& entry, const double curvature, const double bodyLength) {
        if (Kokkos::abs(curvature) <= 1.0e-15 || entry(2) <= 0.0) {
            return entry;  // straight body, or upstream entrance tangent (s = z)
        }

        // Past the exit face: continue along the straight exit tangent.
        const double exitPhi    = curvature * bodyLength;
        const double cosExit    = Kokkos::cos(exitPhi);
        const double sinExit    = Kokkos::sin(exitPhi);
        const double dx         = entry(0) - (cosExit - 1.0) / curvature;
        const double dz         = entry(2) - sinExit / curvature;
        const double exitLocalX = cosExit * dx + sinExit * dz;
        const double exitLocalZ = -sinExit * dx + cosExit * dz;
        if (exitLocalZ >= 0.0) {
            return Vector_t<double, 3>(exitLocalX, entry(1), bodyLength + exitLocalZ);
        }

        // Inside the curved body: radial offset and arc length from the phase.
        const double radius = 1.0 / curvature;
        const double s      = referencePhase(entry, curvature) / curvature;
        const double radial = Kokkos::hypot(entry(0) + radius, entry(2)) - Kokkos::abs(radius);
        return Vector_t<double, 3>(radial, entry(1), s);
    }

    /**
     * @brief Rotate a vector from the local tangent basis at arc length s into the
     * entrance frame.
     *
     * A field evaluated in the basis tangent to the design arc at s (where x is
     * radial) must be rotated by R_y(curvature*s) to be accumulated in the element
     * entrance frame — the tangent has turned by that angle relative to the
     * entrance. Upstream of the entrance no rotation; downstream of the exit the
     * angle is frozen at the exit. The vertical component is invariant.
     */
    KOKKOS_INLINE_FUNCTION Vector_t<double, 3> rotateArcFieldToEntry(
            const Vector_t<double, 3>& B, const double s, const double curvature,
            const double bodyLength) {
        if (Kokkos::abs(curvature) <= 1.0e-15 || s <= 0.0) {
            return B;
        }
        const double phi = curvature * ((s >= bodyLength) ? bodyLength : s);
        const double c   = Kokkos::cos(phi);
        const double sn  = Kokkos::sin(phi);
        return Vector_t<double, 3>(c * B(0) - sn * B(2), B(1), sn * B(0) + c * B(2));
    }

}  // namespace GeometryHelper

#endif  // OPALX_Geometry_HH
