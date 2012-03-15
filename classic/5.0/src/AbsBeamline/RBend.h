#ifndef CLASSIC_RBend_HH
#define CLASSIC_RBend_HH

// ------------------------------------------------------------------------
// $RCSfile: RBend.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBend
//   Defines the abstract interface for a rectangular bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Fields/BMultipoleField.h"


// Class RBend
// ------------------------------------------------------------------------
/// Interface for rectangular bend.
//  Class RBend defines the abstract interface for rectangular bend magnets.
//  A rectangular bend magnet has a rectilinear geometry about which its
//  multipole components are specified.

class RBend: public Component {

public:

  /// Constructor with given name.
  explicit RBend(const string &name);

  RBend();
  RBend(const RBend &);
  virtual ~RBend();

  /// Apply visitor to RBend.
  virtual void accept(BeamlineVisitor &) const;

  /// Get dipole field of RBend.
  virtual double getB() const = 0;

  /// Get RBend geometry.
  //  Version for non-constant object.
  virtual RBendGeometry &getGeometry() = 0;

  /// Get RBend geometry
  //  Version for constant object.
  virtual const RBendGeometry &getGeometry() const = 0;

  /// Get multipole expansion of field.
  //  Version for non-constant object.
  virtual BMultipoleField &getField() = 0;

  /// Get multipole expansion of field.
  //  Version for constant object.
  virtual const BMultipoleField &getField() const = 0;

  /// Get normal component.
  //  Return the normal component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the return value is zero.
  double getNormalComponent(int) const;

  /// Get skew component.
  //  Return the skew component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the return value is zero.
  double getSkewComponent(int) const;

  /// Set normal component.
  //  Set the normal component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the component is created.
  void setNormalComponent(int, double);

  /// Set skew component.
  //  Set the skew component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the component is created.
  void setSkewComponent(int, double);

  /// Get pole entry face rotation.
  //  Return the rotation of the entry pole face with respect to the x-axis.
  //  A positive angle rotates the pole face normal away from the centre
  //  of the machine.
  virtual double getEntryFaceRotation() const = 0;

  /// Get exit pole face rotation.
  //  Return the rotation of the exit pole face with respect to the x-axis.
  //  A positive angle rotates the pole face normal away from the centre
  //  of the machine.
  virtual double getExitFaceRotation() const = 0;

  /// Get entry pole face curvature.
  //  Return the curvature of the entry pole face.
  //  A positive curvature creates a convex pole face.
  virtual double getEntryFaceCurvature() const = 0;
  
  /// Get exit pole face curvature.
  //  Return the curvature of the exit pole face.
  //  A positive curvature creates a convex pole face.
  virtual double getExitFaceCurvature() const = 0;

  /// Get number of slices.
  //  Slices and stepsize used to determine integration step.
  virtual double getSlices() const = 0;

  /// Get stepsize.
  //  Slices and stepsize used to determine integration step.
  virtual double getStepsize() const = 0;

private:

  // Not implemented.
  void operator=(const RBend &);
};

#endif // CLASSIC_RBend_HH
