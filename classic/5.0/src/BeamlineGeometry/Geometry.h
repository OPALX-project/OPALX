#ifndef CLASSIC_Geometry_HH
#define CLASSIC_Geometry_HH

// ------------------------------------------------------------------------
// $RCSfile: Geometry.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Geometry
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

// Euclid3D represents an aribitrary 3-d rotation and displacement.
class Euclid3D;


// Class Geometry
// ------------------------------------------------------------------------
/// Abstract base class for accelerator geometry classes.
//  A Geometry can be considered a 3-dimensional space line parameterised by
//  the distance along the line (arc length) s. All Geometries have an exit
//  and entrance plane and an origin. At any position s, a Geometry can define
//  a unique 3-d rectilinear coordinate frame whose origin is on the geometry
//  at s, and whose local z-axis is tangential to the geometry at s. The
//  orientation of the local x- and y-axes are arbitrarilly specified by
//  the Geometry. A special frame, referred to as the Geometry Local Frame
//  (or Local Frame when it is unambiguous) is specified at s = origin. The
//  Local Frame is is used to define that frame about which translations and
//  rotations can be applied to the Geometry. The entrance and exit planes
//  are defined as those x-y planes (z=0, s=constant) in the frames defined
//  at s=entrance and s=exit.

class Geometry {
public:

    Geometry();
    Geometry(const Geometry &right);
    virtual ~Geometry();
    const Geometry &operator=(const Geometry &right);

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
    virtual Euclid3D getExitPatch()     const;
};

// inlined (trivial) member functions
inline Geometry::Geometry()
{ }

inline Geometry::Geometry(const Geometry &)
{ }

inline const Geometry &Geometry::operator=(const Geometry &)
{ return *this; }

#endif // CLASSIC_Geometry_HH
