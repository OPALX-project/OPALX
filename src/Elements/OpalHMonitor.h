#ifndef OPAL_OpalHMonitor_HH
#define OPAL_OpalHMonitor_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalHMonitor.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalHMonitor
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalHMonitor
// ------------------------------------------------------------------------
/// The HMONITOR element.

class OpalHMonitor: public OpalElement {

public:

  /// Exemplar constructor.
  OpalHMonitor();

  virtual ~OpalHMonitor();

  /// Make clone.
  virtual OpalHMonitor *clone(const string &name);

  /// Update the embedded CLASSIC monitor.
  virtual void update();

private:

  // Not implemented.
  OpalHMonitor(const OpalHMonitor &);
  void operator=(const OpalHMonitor &);

  // Clone constructor.
  OpalHMonitor(const string &name, OpalHMonitor *parent);
};

#endif // OPAL_OpalHMonitor_HH
