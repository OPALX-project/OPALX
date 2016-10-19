#ifndef CLASSIC_SRotatedGeometry_HH
#define CLASSIC_SRotatedGeometry_HH

// ------------------------------------------------------------------------
// $RCSfile: SRotatedGeometry.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SRotatedGeometry
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

// Class SRotatedGeometry
// ------------------------------------------------------------------------
/// A Geometry which wraps an arbitrary geometry in two s-rotations.
//  An SRotatedGeometry object is a Geometry wrapper which adds two
//  arbitrary s-rotations (local z-axis rotations) to the entrance and
//  exit planes of an arbitrary geometry. Any Geometry object can be
//  wrapped, including another SRotatedGeometry. The s-rotations become
//  part of the global geometry definition. Functions for setting the
//  two in- and out-rotations using certain constraints are provided.
//  [P]
//  NOTE: in general the transformations returned include the effects of
//  the s-rotations when the required distance parameter specified is
//  either the entrance or exit point. Requests for transformations
//  within the geometry (i.e. from s1 to s2, where s1 and/or s2 are not
//  the entrance or exit planes) do not contain the s-rotations.

class SRotatedGeometry : public BGeometryBase {
public:

    /// Balance mode.
    enum BalanceMode {

        /// srotOut is set to -srotIn.
        tilt,

        /// srotOut calculated.
        //  srotOut is calculated to be the vlaue which returns the local
        //  x-axis after rotation to be in a plane parallel with the local
        //  x-axis before application of srotIn.
        balanceX

    };

    /// Constructor.
    //  Use the wrapped geometry [b]geom[/b],
    //  and the two angles [b]srotIn[/b] and [b]srotOut[/b].
    SRotatedGeometry(const BGeometryBase &geom, double srotIn, double srotOut);

    /// Constructor.
    //  Use the wrapped geometry [b]geom[/b],
    //  the entrance rotation [b]srotIn[/b],
    //  and the balanc mode [b]mode[/b].
    SRotatedGeometry(const BGeometryBase &geom, double srotIn, BalanceMode mode = tilt);

    SRotatedGeometry(const SRotatedGeometry &);
    virtual ~SRotatedGeometry();
    const SRotatedGeometry &operator=(const SRotatedGeometry &);

    /// Get arc length.
    //  Return or the design length of the embedded geometry.
    //  This is the length of the straight line connecting entrance and exit.
    virtual double getArcLength() const;

    /// Get design length.
    //  Return or the design length of the embedded geometry.
    //  Depending on the element this may be the arc length or the straight
    //  length.
    virtual double getElementLength() const;

    /// Get entrance rotation.
    double getSrotIn() const;

    /// Get exit rotation.
    double getSrotOut() const;

    /// Set entrance rotation.
    void setSrotIn(double);

    /// Set exit rotation.
    void setSrotOut(double);

    /// Balance rotations.
    //  Set the exit rotation in one of two modes:
    //  [UL]
    //  [LI]if [tt]mode==tilt[/tt], then srotOut is set to -srotIn.
    //  [LI]if [tt]mode==balanceX[/tt], then srotOut is calculated to be
    //  the value which returns the local x-axis after rotation to be in
    //  a parallel plane with the local x-axis before the application of
    //  [tt]srotIn[/tt].
    //  [/UL]
    void balanceSrots(BalanceMode mode = tilt);

    /// Get origin.
    //  Return the arc length from the entrance to the origin of the
    //  geometry (non-negative).
    double getOrigin() const;

    /// Get entrance.
    //  Return the arc length from the origin to the entrance of the
    //  geometry (non-positive)
    double getEntrance() const;

    /// Get exit.
    //  Return the arc length from the origin to the exit of the
    //  geometry (non-negative)
    double getExit() const;

    /// Get transform.
    //  Return the transform of the local coordinate system from the
    //  position [b]fromS[/b] to the position [b]toS[/b].
    Euclid3D getTransform(double fromS, double toS) const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, s).
    //  Return the transform of the local coordinate system from the
    //  origin and [b]s[/b].
    Euclid3D getTransform(double s) const;

    /// Get transform.
    //  Equivalent to getTransform(getEntrance(), getExit()).
    //  Return the transform of the local coordinate system from the
    //  entrance to the exit of the element.
    Euclid3D getTotalTransform() const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, getEntrance()).
    //  Return the transform of the local coordinate system from the
    //  origin to the entrance of the element.
    Euclid3D getEntranceFrame() const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, getExit()).
    //  Return the transform of the local coordinate system from the
    //  origin to the exit of the element.
    Euclid3D getExitFrame()     const;

    /// Get patch.
    //  Returns the entrance patch (transformation) which is used to
    //  transform the global geometry to the local geometry at entrance.
    Euclid3D getEntrancePatch() const;

    /// Get patch.
    //  Returns the entrance patch (transformation) which is used to
    //  transform the local geometry to the global geometry at exit.
    Euclid3D getExitPatch()     const;

private:

    double srotIn;
    double srotOut;
    const BGeometryBase &geom;
};

#endif // CLASSIC_SRotatedGeometry_HH

