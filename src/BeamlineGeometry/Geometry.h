#ifndef OPALX_BGeometryBase_HH
#define OPALX_BGeometryBase_HH

// ------------------------------------------------------------------------
// $RCSfile: BGeometryBase.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BGeometryBase
//
// ------------------------------------------------------------------------
// Class category: BeamlineBGeometryBase
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

// Euclid3D represents an aribitrary 3-d rotation and displacement.
#include "BeamlineGeometry/Euclid3D.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "OPALTypes.h"

#include <cstddef>
#include <vector>

// Class BGeometryBase
// ------------------------------------------------------------------------
/// Abstract base class for accelerator geometry classes.
//  A BGeometryBase can be considered a 3-dimensional space line parameterised by
//  the distance along the line (arc length) s. All Geometries have an exit
//  and entrance plane and an origin. At any position s, a BGeometryBase can define
//  a unique 3-d rectilinear coordinate frame whose origin is on the geometry
//  at s, and whose local z-axis is tangential to the geometry at s. The
//  orientation of the local x- and y-axes are arbitrarilly specified by
//  the BGeometryBase. A special frame, referred to as the BGeometryBase Local Frame
//  (or Local Frame when it is unambiguous) is specified at s = origin. The
//  Local Frame is is used to define that frame about which translations and
//  rotations can be applied to the BGeometryBase. The entrance and exit planes
//  are defined as those x-y planes (z=0, s=constant) in the frames defined
//  at s=entrance and s=exit.

class BGeometryBase {
public:
    BGeometryBase();
    BGeometryBase(const BGeometryBase& right);
    virtual ~BGeometryBase();
    const BGeometryBase& operator=(const BGeometryBase& right);

    /// Get arc length.
    //  Return the length of the geometry, measured along the design arc.
    virtual double getArcLength() const = 0;

    /// Get geometry length.
    //  Return or the design length of the geometry.
    //  Depending on the element this may be the arc length or the
    //  straight length.
    virtual double getElementLength() const = 0;

    /// Set geometry length.
    //  Assign the design length of the geometry.
    //  Depending on the element this may be the arc length or the
    //  straight length.
    virtual void setElementLength(double length);

    /// Get origin position.
    //  Return the arc length from the entrance to the origin of the
    //  geometry (non-negative).
    virtual double getOrigin() const;

    /// Get entrance position.
    //  Return the arc length from the origin to the entrance of the
    //  geometry (non-positive).
    virtual double getEntrance() const;

    /// Get exit position.
    //  Return the arc length from the origin to the exit of the
    //  geometry (non-negative).
    virtual double getExit() const;

    /// Get transform.
    //  Return the transform of the local coordinate system from the
    //  position [b]fromS[/b] to the position [b]toS[/b].
    virtual Euclid3D getTransform(double fromS, double toS) const = 0;

    /// Get transform.
    //  Equivalent to getTransform(0.0, s).
    //  Return the transform of the local coordinate system from the
    //  origin and [b]s[/b].
    virtual Euclid3D getTransform(double s) const;

    /// Get transform.
    //  Equivalent to getTransform(getEntrance(), getExit()).
    //  Return the transform of the local coordinate system from the
    //  entrance to the exit of the element.
    virtual Euclid3D getTotalTransform() const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, getEntrance()).
    //  Return the transform of the local coordinate system from the
    //  origin to the entrance of the element.
    virtual Euclid3D getEntranceFrame() const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, getExit()).
    //  Return the transform of the local coordinate system from the
    //  origin to the exit of the element.
    virtual Euclid3D getExitFrame() const;

    /// Get patch.
    //  Returns the entrance patch (transformation) which is used to transform
    //  the global geometry to the local geometry for a misaligned element
    //  at its entrance. The default behaviour returns identity transformation.
    //  This function should be overidden by derived concrete classes which
    //  model complex geometries.
    virtual Euclid3D getEntrancePatch() const;

    /// Get patch.
    //  Returns the entrance patch (transformation) which is used to transform
    //  the local geometry to the global geometry for a misaligned element
    //  at its exit. The default behaviour returns identity transformation.
    //  This function should be overidden by derived concrete classes which
    //  model complex geometries.
    virtual Euclid3D getExitPatch() const;
};

// inlined (trivial) member functions
inline BGeometryBase::BGeometryBase() {}

inline BGeometryBase::BGeometryBase(const BGeometryBase&) {}

inline const BGeometryBase& BGeometryBase::operator=(const BGeometryBase&) { return *this; }

// Class Geometry
// ------------------------------------------------------------------------
/// The kind of body geometry represented by a Geometry object.
enum class GeometryKind : unsigned char { Null, Straight, Arc, RBend };

/**
 * @class Geometry
 * @brief Single geometry class for all beamline elements.
 *
 * Replaces the BGeometryBase class tree (StraightGeometry, PlanarArcGeometry,
 * RBendGeometry, NullGeometry) with one concrete, mode-tagged class. It carries
 * the body length plus the bend parameters (curvature/angle and entrance/exit
 * pole-face angles) and provides:
 *   - the body length and arc length,
 *   - the bend angle, curvature and chord length,
 *   - the entrance/exit edge transforms as CoordinateSystemTrafo (the single
 *     source for placement), and
 *   - the sampled design path.
 *
 * @note A GeometryKind tag selects straight / planar-arc / rectangular-bend /
 *       null behaviour. For straight and null geometries the edge transforms are
 *       the identity at the entrance and a +z shift to the exit, matching the
 *       former ElementBase defaults. For bends they reproduce the planar-arc and
 *       rectangular-bend pole-face frames of the former geometry subclasses.
 */
class Geometry : public BGeometryBase {
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
    double getArcLength() const override;
    /// Design / body length (straight body length for rectangular bends).
    double getElementLength() const override { return len_m; }
    void setElementLength(double length) override;
    ///@}

    /// @name Reference points along the body chart (centred at the origin)
    ///@{
    double getOrigin() const override { return len_m / 2.0; }
    double getEntrance() const override { return -len_m / 2.0; }
    double getExit() const override { return len_m / 2.0; }
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

    /// @name Legacy Euclid3D frame API (transitional, removed once bends are
    /// folded; reproduces the former geometry subclasses for BendBase/TBeamline).
    ///@{
    Euclid3D getTransform(double fromS, double toS) const override;
    Euclid3D getTransform(double s) const override;
    Euclid3D getTotalTransform() const override;
    Euclid3D getEntranceFrame() const override;
    Euclid3D getExitFrame() const override;
    ///@}

private:
    double halfAngle() const { return angle_m / 2.0; }

    // Default-constructed geometry is a zero-length straight body (the common
    // case for newly created representations); markers use makeNull().
    GeometryKind kind_m = GeometryKind::Straight;
    double len_m           = 0.0;  ///< design / body length
    double h_m             = 0.0;  ///< curvature (Arc)
    double angle_m         = 0.0;  ///< bend angle (Arc: h*len; RBend: full angle)
    double entranceAngle_m = 0.0;  ///< entrance pole-face angle
    double exitAngle_m     = 0.0;  ///< exit pole-face angle
};

#endif  // OPALX_BGeometryBase_HH
