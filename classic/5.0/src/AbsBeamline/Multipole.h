#ifndef CLASSIC_Multipole_HH
#define CLASSIC_Multipole_HH

// ------------------------------------------------------------------------
// $RCSfile: Multipole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/BMultipoleField.h"


// Class Multipole
// ------------------------------------------------------------------------
/// Interface for general multipole.
//  Class Multipole defines the abstract interface for magnetic multipoles.
//  The order n of multipole components runs from 1 to N and is dynamically
//  adjusted. It is connected with the number of poles by the table
//
//  [tab 2 b]
//  [ROW]1[&]dipole[/ROW]
//  [ROW]2[&]quadrupole[/ROW]
//  [ROW]3[&]sextupole[/ROW]
//  [ROW]4[&]octupole[/ROW]
//  [ROW]5[&]decapole[/ROW]
//  [ROW]n[&]multipole with 2*n poles[/ROW]
//  [/TAB]
//  Units for multipole strengths are Teslas / m**(n-1).

class Multipole: public Component {

public:

  /// Constructor with given name.
  explicit Multipole(const string &name);

  Multipole();
  Multipole(const Multipole &);
  virtual ~Multipole();

  /// Apply visitor to Multipole.
  virtual void accept(BeamlineVisitor &) const;


  /// Get multipole field.
  virtual BMultipoleField &getField() = 0;

  /// Get multipole field. Version for const object.
  virtual const BMultipoleField &getField() const = 0;

  /// Get normal component.
  //  Return the normal component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the return value is zero.
  double getNormalComponent(int n) const;

  /// Get skew component.
  //  Return the skew component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the return value is zero.
  double getSkewComponent(int n) const;

  /// Set normal component.
  //  Set the normal component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the component is created.
  void setNormalComponent(int, double);

  /// Set skew component.
  //  Set the skew component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the component is created.
  void setSkewComponent(int, double);

  /// Get geometry.
  virtual StraightGeometry &getGeometry() = 0;

  /// Get geometry.
  virtual const StraightGeometry &getGeometry() const = 0;

private:

  // Not implemented.
  void operator=(const Multipole &);
};

#endif // CLASSIC_Multipole_HH
