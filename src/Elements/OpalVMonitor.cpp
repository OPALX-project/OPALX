// ------------------------------------------------------------------------
// $RCSfile: OpalVMonitor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalVMonitor
//   The class of OPAL monitors for the vertical plane.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalVMonitor.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/YMonitorRep.h"


// Class OpalVMonitor
// ------------------------------------------------------------------------

OpalVMonitor::OpalVMonitor():
    OpalElement(COMMON, "VMONITOR",
                "The \"VMONITOR\" element defines a monitor "
                "for the vertical plane.") {
    setElement((new YMonitorRep("VMONITOR"))->makeAlignWrapper());
}


OpalVMonitor::OpalVMonitor(const string &name, OpalVMonitor *parent):
    OpalElement(name, parent) {
    setElement((new YMonitorRep(name))->makeAlignWrapper());
}


OpalVMonitor::~OpalVMonitor()
{}


OpalVMonitor *OpalVMonitor::clone(const string &name) {
    return new OpalVMonitor(name, this);
}


void OpalVMonitor::update() {
    YMonitorRep *mon =
        dynamic_cast<YMonitorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    mon->setElementLength(length);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mon);
}
