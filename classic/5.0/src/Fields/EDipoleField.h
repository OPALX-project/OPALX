#ifndef CLASSIC_EDipoleField_HH
#define CLASSIC_EDipoleField_HH

// ------------------------------------------------------------------------
// $RCSfile: EDipoleField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EDipoleField
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/ConstEField.h"


// Class EDipoleField
// ------------------------------------------------------------------------
/// A static homogeneous electrostatic field in the (x,y)-plane.

class EDipoleField: public ConstEField {

public:

  /// Default constructor.
  //  Constructs a null field.
  EDipoleField();

  virtual ~EDipoleField();

  /// Get field.
  //  Return the time-independent part of the electric field in point [b]P[/b].
  virtual EVector Efield(const Point3D &P) const;

  /// Get field.
  //  Return the electric field at time [b]t[/b] in point [b]P[/b].
  virtual EVector Efield(const Point3D &P, double t) const;

  /// Get component.
  //  Return the x-component of the electric field in A/m.
  virtual double getEx() const;

  /// Get component.
  //  Return the y-component of the electric field in A/m.
  virtual double getEy() const;

  /// Set component.
  //  Assign the x-component of the electric field in A/m.
  virtual void setEx(double);

  /// Set component.
  //  Assign the y-component of the electric field in A/m.
  virtual void setEy(double);

  /// Scale the field.
  //  Multiply the field by [b]scalar[/b].
  virtual void scale(double scalar);

private:

  // The field components.
  double Ex, Ey;
};

#endif // CLASSIC_EDipoleField_HH
