#ifndef CLASSIC_ConstEzField_HH
#define CLASSIC_ConstEzField_HH

// ------------------------------------------------------------------------
// $RCSfile: ConstEzField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstEzField
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


// Class ConstEzField
// ------------------------------------------------------------------------
/// A homogeneous electrostatic field in z-direction.

class ConstEzField: public StaticElectricField {

public:

  /// Default constructor.
  //  Constructs a null field.
  ConstEzField();

  virtual ~ConstEzField();

  /// Get field.
  //  Return the time-independent part of the electric field in point [b]P[/b].
  virtual EVector Efield(const Point3D &P) const;

  /// Get field.
  //  Return the electric field at time [b]t[/b] in point [b]P[/b].
  virtual EVector Efield(const Point3D &P, double t) const;

  /// Get component.
  //  Return the x-component of the electric field in A/m.
  virtual double getEz() const;

  /// Set component.
  //  Assign the z-component of the electric field in A/m.
  virtual void setEz(double);

  /// Scale the field.
  //  Multiply the field by [b]scalar[/b].
  virtual void scale(double scalar);

private:

  // The field components.
  double Ez;
};

#endif // CLASSIC_ConstEzField_HH
