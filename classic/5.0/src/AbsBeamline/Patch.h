#ifndef CLASSIC_Patch_HH
#define CLASSIC_Patch_HH

// ------------------------------------------------------------------------
// $RCSfile: Patch.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Patch
//   Defines the abstract interface for a geometry patch.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"


// Class Patch
// ------------------------------------------------------------------------
/// Interface for a geometric patch.
//  Class Patch defines the abstract interface for geometric patches,
//  i.e. a geometric transform representing a misalignment of local 
//  coordinate systems of two subsequent element with each other.

class Patch: public Component {

public:

  /// Constructor with given name.
  explicit Patch(const string &name);

  Patch();
  Patch(const Patch &);
  virtual ~Patch();

  /// Apply visitor to patch.
  virtual void accept(BeamlineVisitor &) const;

  /// Get patch transform.
  virtual const Euclid3D &getPatch() const = 0;

  virtual bool apply(const int &i, const double &t, double E[], double B[]);

  virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);
  
  virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);
  
  virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

  virtual void finalise();

  virtual bool bends() const;

  virtual  const string& getType() const;

  virtual void getDimensions(double &zBegin, double &zEnd) const;

private:

  // Not implemented.
  void operator=(const Patch &);
};

#endif // CLASSIC_Patch_HH
