#ifndef OPAL_OpalDrift_HH
#define OPAL_OpalDrift_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalDrift.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalDrift
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalDrift
// ------------------------------------------------------------------------
/// The DRIFT element.

class OpalDrift: public OpalElement {

public:

  /// Exemplar constructor.
  OpalDrift();

  virtual ~OpalDrift();

  /// Make clone.
  virtual OpalDrift *clone(const string &name);

  /// Test for drift.
  //  Return true.
  virtual bool isDrift() const;

  /// Update the embedded CLASSIC drift.
  virtual void update();

private:

  // Not implemented.
  OpalDrift(const OpalDrift &);
  void operator=(const OpalDrift &);

  // Clone constructor.
  OpalDrift(const string &name, OpalDrift *parent);
};

#endif // OPAL_OpalDrift_HH
