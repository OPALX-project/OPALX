#include "Utilities/OpalField.h"

extern Inform *gmsg;

OpalField::OpalField(Component *element, const double &start, const double &end):
    element_m(element),
    start_m(start),
    end_m(end),
    is_on_m(false)
{ }

OpalField::~OpalField() {
    element_m = NULL;
}

void OpalField::setOn() {
    if(!is_on_m) {
        element_m->goOnline();
        INFOMSG(element_m->getName() << " gone live" << endl);
        is_on_m = true;
    }
}

void OpalField::setOff() {
    if(is_on_m) {
        element_m->goOffline();
        INFOMSG(element_m->getName() << " gone off" << endl);
        is_on_m = false;
    }
}
