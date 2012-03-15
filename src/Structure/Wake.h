#ifndef OPAL_Wake_HH
#define OPAL_Wake_HH

// ------------------------------------------------------------------------
// $RCSfile: Wake.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Wake
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Solvers/WakeFunction.hh"
class ElementBase;

// Class Wake
// ------------------------------------------------------------------------
/// The WAKE definition.
//  A WAKE definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.

class Wake: public Definition {

public:

  /// Exemplar constructor.
  Wake();

  virtual ~Wake();

  /// Test if replacement is allowed.
  //  Can replace only by another WAKE.
  virtual bool canReplaceBy(Object *object);

  /// Make clone.
  virtual Wake *clone(const string &name);

  /// Check the WAKE data.
  virtual void execute();

  /// Find named WAKE.
  static Wake *find(const string &name);

  /// Update the WAKE data.
  virtual void update();

  /// Print the TFS descriptors for the wake
  void tfsDescriptors(std::ostream &os) const;

  Inform &print(Inform &os) const; 

  int getNumberOfBins();

  void initWakefunction(ElementBase &element);

  WakeFunction *wf_m;

private:

  // Not implemented.
  Wake(const Wake &);
  void operator=(const Wake &);

  // Clone constructor.
  Wake(const string &name, Wake *parent);

  // The particle reference data.
  PartData reference;

  // The conversion from GeV to eV.
  static const double energy_scale;

  // the element the wake is attached to
  ElementBase *itsElement_m;

};

inline Inform &operator<<(Inform &os, const Wake &b)
{
  return b.print(os);
}

#endif // OPAL_Wake_HH
