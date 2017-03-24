// ------------------------------------------------------------------------
// $RCSfile: OpalHMonitor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalHMonitor
//   The class of OPAL monitors for the horizontal plane.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalHMonitor.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/XMonitorRep.h"


// Class OpalHMonitor
// ------------------------------------------------------------------------


OpalHMonitor::OpalHMonitor():
    OpalElement(COMMON, "HMONITOR",
                "The \"HMONITOR\" element defines a monitor "
                "for the horizontal plane.") {
    setElement((new XMonitorRep("HMONITOR"))->makeAlignWrapper());
}


OpalHMonitor::OpalHMonitor(const std::string &name, OpalHMonitor *parent):
    OpalElement(name, parent) {
    setElement((new XMonitorRep(name))->makeAlignWrapper());
}


OpalHMonitor::~OpalHMonitor()
{}


OpalHMonitor *OpalHMonitor::clone(const std::string &name) {
    return new OpalHMonitor(name, this);
}


void OpalHMonitor::update() {
    OpalElement::update();

    XMonitorRep *mon =
        dynamic_cast<XMonitorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    mon->setElementLength(length);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mon);
}