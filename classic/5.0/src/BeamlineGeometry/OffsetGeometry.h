#ifndef CLASSIC_OffsetGeometry_HH
#define CLASSIC_OffsetGeometry_HH

// ------------------------------------------------------------------------
// $RCSfile: OffsetGeometry.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OffsetGeometry
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/Geometry.h"
#include "BeamlineGeometry/Euclid3D.h"

// Class OffsetGeometry
// ------------------------------------------------------------------------
/// A geometry which offset with respect to some global geometry.
//  OffsetGeometry effectively acts as a bridge between two geometries,
//  designated as the global geometry and the local geometry. The global
//  geometry is used to define a segment of the true global geometry, and
//  so most functions are delegated to the global geometry. The local and
//  global geometries are ``offset'' with respect to each other by a single
//  transformation (Euclid3D object), which defines the transformation from
//  the global geometry's Local Frame to the local geometry's Local Frame.
//  The two patch functions return the necessary transformations to/from the
//  global/local geometries. OffsetGeometry primary task is to calculate
//  these two patches. Note that when OffsetGeometry is being used to model
//  an alignment error, the local and global geometries can both refer to
//  the same geometry.

class OffsetGeometry : public BGeometryBase {
public:

    /// Constructor.
    //  Assign the [b]global[/b] and [b]local[/b] geometries,
    //  and the displacement [b]euclid[/b] between their origins.
    OffsetGeometry
    (const BGeometryBase &global, const BGeometryBase &local, const Euclid3D &euclid);

    /// Constructor.
    //  Both global and local geometries are set to [b]geom[/b],
    //  and the displacement is [b]euclid[/b].
    OffsetGeometry(const BGeometryBase &geom, const Euclid3D &euclid);

    OffsetGeometry(const OffsetGeometry &);
    virtual ~OffsetGeometry();
    const OffsetGeometry &operator=(const OffsetGeometry &);

    /// Get displacement.
    //  Return the displacement from the global to the local origin.
    Euclid3D getGtoL() const;

    /// Set displacement.
    //  Assign the displacement from the global to the local origin.
    void setGtoL(const Euclid3D &);

    /// Get arc length.
    //  Return the length of the global geometry,
    //  measured along its design arc.
    virtual double getArcLength() const;

    /// Get design length.
    //  Return the design length of the global geometry,
    //  measured along its design polygone.
    virtual double getElementLength() const;

    /// Get origin.
    //  Return the arc length from the entrance to the origin of the
    //  global geometry (non-negative).
    virtual double getOrigin() const;

    /// Get entrance position.
    //  Return the arc length from the origin to the entrance of the
    //  global geometry (non-positive).
    virtual double getEntrance() const;

    /// Get exit position.
    //  Return the arc length from the origin to the exit of the
    // global geometry (non-negative).
    virtual double getExit() const;

    /// Get transform.
    //  Return the transform of the global coordinate system from the
    //  position [b]fromS[/b] to the position [b]toS[/b].
    virtual Euclid3D getTransform(double fromS, double toS) const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, s).
    //  Return the transform of the local coordinate system from the
    //  global origin at [b]s[/b].
    virtual Euclid3D getTransform(double s) const;

    /// Get transform.
    //  Equivalent to getTransform(getEntrance(), getExit()).
    //  Return the transform of the global coordinate system from the
    //  entrance to the exit of the element.
    virtual Euclid3D getTotalTransform() const;


    /// Get transform.
    //  Equivalent to getTransform(0.0, getEntrance()).
    //  Return the transform of the local coordinate system from the
    //  global origin to the entrance of the element.
    virtual Euclid3D getEntranceFrame() const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, getExit()).
    //  Return the transform of the local coordinate system from the
    //  global origin to the exit of the element.
    virtual Euclid3D getExitFrame() const;

    /// Get patch.
    //  Return the entrance patch (transformation) which is used to transform
    //  the global geometry to the local geometry at entrance.
    virtual Euclid3D getEntrancePatch() const;

    /// Get patch.
    //  Returns the entrance patch (transformation) which is used to transform
    //  the local geometry to the global geometry at exit.
    virtual Euclid3D getExitPatch() const;

    /// Transform global to local.
    //  Special OffsetGeometry function which calculates the transformation
    //  from position [b]globalS[/b] on the global geometry to position
    //  [b]localS[/b] on the local geometry.
    Euclid3D getGlobalToLocalTransform(double globalS, double localS) const;

private:

    const BGeometryBase &global;
    const BGeometryBase &local;
    Euclid3D g2l;
};

#endif // CLASSIC_OffsetGeometry_HH

