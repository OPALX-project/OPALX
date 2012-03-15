#ifndef OPAL_BOUNDARY_GEOMETRY_HH
#define OPAL_BOUNDARY_GEOMETRY_HH

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
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
class ElementBase;

// Class BoundaryGeometry
// ------------------------------------------------------------------------
/// The GEOMETRY definition.
//  A GEOMETRY definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.
//
//  i.e:
//  G1: Geometry, FILE="input.h5"
//  G2: Geometry, L=1.0, A=0.0025, B=0.0001
//

class BoundaryGeometry: public Definition {

public:

  /// Exemplar constructor.
  BoundaryGeometry();

  virtual ~BoundaryGeometry();

  /// Test if replacement is allowed.
  //  Can replace only by another GEOMETRY.
  virtual bool canReplaceBy(Object *object);

  /// Make clone.
  virtual BoundaryGeometry *clone(const string &name);

  /// Check the GEOMETRY data.
  virtual void execute();

  /// Find named GEOMETRY.
  static BoundaryGeometry *find(const string &name);

  /// Update the GEOMETRY data.
  virtual void update();

  Inform &print(Inform &os) const; 
  string getTopology();
  string getFilename();
  double getA();
  double getB();

private:

  // Not implemented.
  BoundaryGeometry(const BoundaryGeometry &);
  void operator=(const BoundaryGeometry &);

  // Clone constructor.
  BoundaryGeometry(const string &name, BoundaryGeometry *parent);

  // The particle reference data.
  PartData reference;

};

inline Inform &operator<<(Inform &os, const BoundaryGeometry &b)
{
  return b.print(os);
}

#endif // OPAL_BOUNDARY_GEOMETRY_HH
