#include "Utilities/ClassicField.h"

extern Inform *gmsg;

ClassicField::ClassicField(std::shared_ptr<Component> element, const double &start, const double &end):
    element_m(element),
    start_m(start),
    end_m(end),
    is_on_m(false)
{ }

ClassicField::~ClassicField() {
    element_m = NULL;
}

void ClassicField::setOn() {
    if(!is_on_m) {
        element_m->goOnline();
        INFOMSG(element_m->getName() << " gone live" << endl);
        is_on_m = true;
    }
}

void ClassicField::setOff() {
    if(is_on_m) {
        element_m->goOffline();
        INFOMSG(element_m->getName() << " gone off" << endl);
        is_on_m = false;
    }
}
