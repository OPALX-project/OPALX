// ------------------------------------------------------------------------
// $RCSfile: YMonitorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: YMonitorRep
//   Defines a concrete representation for a vertical orbit monitor.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/YMonitorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(YMonitorRep::*get)() const;
        void (YMonitorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &YMonitorRep::getElementLength,
            &YMonitorRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class YMonitorRep
// ------------------------------------------------------------------------

YMonitorRep::YMonitorRep():
    MonitorRep()
{}


YMonitorRep::YMonitorRep(const YMonitorRep &right):
    MonitorRep(right)
{}


YMonitorRep::YMonitorRep(const string &name):
    MonitorRep(name)
{}


YMonitorRep::~YMonitorRep()
{}


ElementBase *YMonitorRep::clone() const {
    return new YMonitorRep(*this);
}


Channel *YMonitorRep::getChannel(const string &aKey) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<YMonitorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey);
}


ElementImage *YMonitorRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}


Monitor::Plane YMonitorRep::getPlane() const {
    return active ? Y : OFF;
}


