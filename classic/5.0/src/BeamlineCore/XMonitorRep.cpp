// ------------------------------------------------------------------------
// $RCSfile: XMonitorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: XMonitorRep
//   Defines a concrete representation for a horizontal orbit monitor.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/XMonitorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(XMonitorRep::*get)() const;
        void (XMonitorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &XMonitorRep::getElementLength,
            &XMonitorRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class XMonitorRep
// ------------------------------------------------------------------------

XMonitorRep::XMonitorRep():
    MonitorRep()
{}


XMonitorRep::XMonitorRep(const XMonitorRep &rhs):
    MonitorRep(rhs)
{}


XMonitorRep::XMonitorRep(const string &name):
    MonitorRep(name)
{}


XMonitorRep::~XMonitorRep()
{}


ElementBase *XMonitorRep::clone() const {
    return new XMonitorRep(*this);
}


Channel *XMonitorRep::getChannel(const string &aKey) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<XMonitorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey);
}


ElementImage *XMonitorRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}


Monitor::Plane XMonitorRep::getPlane() const {
    return active ? X : OFF;
}


