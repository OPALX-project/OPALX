// ------------------------------------------------------------------------
// $RCSfile: OpalMonitor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalMonitor
//   The class of OPAL monitors for both planes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalMonitor.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MonitorRep.h"


// Class OpalMonitor
// ------------------------------------------------------------------------
 
OpalMonitor::OpalMonitor():
    OpalElement(SIZE, "MONITOR",
                "The \"MONITOR\" element defines a monitor for both planes.")
{
    itsAttr[OUTFN] = Attributes::makeString
        ("OUTFN", "Monitor output filename");

    registerStringAttribute("OUTFN");

    setElement((new MonitorRep("MONITOR"))->makeAlignWrapper());
}


OpalMonitor::OpalMonitor(const string &name, OpalMonitor *parent):
    OpalElement(name, parent)
{
    setElement((new MonitorRep(name))->makeAlignWrapper());
}


OpalMonitor::~OpalMonitor()
{}


OpalMonitor *OpalMonitor::clone(const string &name)
{
    return new OpalMonitor(name, this);
}


void OpalMonitor::update()
{
    MonitorRep *mon =
        dynamic_cast<MonitorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    mon->setElementLength(length);
    mon->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mon);
}
