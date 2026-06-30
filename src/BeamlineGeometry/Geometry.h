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

#include <cstddef>
#include <vector>

/// The kind of body geometry represented by a Geometry object.
enum class GeometryKind : unsigned char { Null, Straight, Arc, RBend };

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
    static Geometry makeArc(double length, double curvature);
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

    /// @name Reference points along the body chart (centred at the origin)
    ///@{
    double getOrigin() const { return len_m / 2.0; }
    double getEntrance() const { return -len_m / 2.0; }
    double getExit() const { return len_m / 2.0; }
    ///@}

    /// @name Bend parameters
    ///@{
    double getBendAngle() const;
    void setBendAngle(double angle);
    double getCurvature() const { return (kind_m == GeometryKind::Arc) ? h_m : 0.0; }
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
    /// Body-origin to entrance-edge transform.
    CoordinateSystemTrafo getEdgeToBegin() const;
    /// Body-origin to exit-edge transform.
    CoordinateSystemTrafo getEdgeToEnd() const;
    ///@}

private:
    double halfAngle() const { return angle_m / 2.0; }
    /// Position (x, 0, z) of the centred body chart at arc-position s.
    Vector_t<double, 3> framePosition(double s) const;

    GeometryKind kind_m    = GeometryKind::Straight;
    double len_m           = 0.0;  ///< design / body length
    double h_m             = 0.0;  ///< curvature (Arc)
    double angle_m         = 0.0;  ///< bend angle (Arc: h*len; RBend: full angle)
    double entranceAngle_m = 0.0;  ///< entrance pole-face angle
    double exitAngle_m     = 0.0;  ///< exit pole-face angle
};

#endif  // OPALX_Geometry_HH
