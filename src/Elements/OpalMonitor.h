#ifndef OPAL_OpalMonitor_HH
#define OPAL_OpalMonitor_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalMonitor.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalMonitor
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalMonitor
// ------------------------------------------------------------------------
/// The MONITOR element.

class OpalMonitor: public OpalElement {

public:

  /// Exemplar constructor.
  OpalMonitor();

  virtual ~OpalMonitor();

  /// Make clone.
  virtual OpalMonitor *clone(const string &name);

  /// Update the embedded CLASSIC monitor.
  virtual void update();

private:

  // Not implemented.
  OpalMonitor(const OpalMonitor &);
  void operator=(const OpalMonitor &);

  // Clone constructor.
  OpalMonitor(const string &name, OpalMonitor *parent);
};

#endif // OPAL_OpalMonitor_HH
