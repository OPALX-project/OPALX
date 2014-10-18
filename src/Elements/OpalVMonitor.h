#ifndef OPAL_OpalVMonitor_HH
#define OPAL_OpalVMonitor_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalVMonitor.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalVMonitor
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalVMonitor
// ------------------------------------------------------------------------
/// The VMONITOR element.

class OpalVMonitor: public OpalElement {

public:

    /// Exemplar constructor.
    OpalVMonitor();

    virtual ~OpalVMonitor();

    /// Make clone.
    virtual OpalVMonitor *clone(const std::string &name);

    /// Update the embedded CLASSIC monitor.
    virtual void update();

private:

    // Not implemented.
    OpalVMonitor(const OpalVMonitor &);
    void operator=(const OpalVMonitor &);

    // Clone constructor.
    OpalVMonitor(const std::string &name, OpalVMonitor *parent);
};

#endif // OPAL_OpalVMonitor_HH
