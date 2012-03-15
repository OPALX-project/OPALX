#ifndef OPAL_OpalRBend_HH
#define OPAL_OpalRBend_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalRBend.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalRBend
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/24 19:35:09 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalBend.h"
 

// Class OpalRBend
// ------------------------------------------------------------------------
/// The RBEND element.

class OpalRBend: public OpalBend {

public:

  /// Exemplar constructor.
  OpalRBend();

  virtual ~OpalRBend();

  /// Make clone.
  virtual OpalRBend *clone(const string &name);

  /// Fill in all registered attributes.
  virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

  /// Update the embedded CLASSIC bend.
  virtual void update();

private:

  // Not implemented.
  OpalRBend(const OpalRBend &);
  void operator=(const OpalRBend &);

  // Clone constructor.
  OpalRBend(const string &name, OpalRBend *parent);
};

#endif // OPAL_OpalRBend_HH
