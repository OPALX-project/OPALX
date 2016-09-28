// ------------------------------------------------------------------------
// $RCSfile: Element.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Element
//   The base class for all OPAL beamline elements.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"


// Class Element
// ------------------------------------------------------------------------

Element::~Element()
{}


bool Element::canReplaceBy(Object *object) {
    return (dynamic_cast<Element *>(object) != 0);
}


Element *Element::find(const std::string &name) {
    OpalData *opal = OpalData::getInstance();
    Element *element = dynamic_cast<Element *>(opal->find(name));
    if(element == 0) {
        throw OpalException("Element::find()",
                            "Element \"" + name + "\" not found.");
    }
    return element;
}


const std::string Element::getCategory() const {
    return "ELEMENT";
}


bool Element::shouldTrace() const {
    return false;

}

bool Element::shouldUpdate() const {
    return false;
}


double Element::getEntrance(ReferenceType ref) const {
    switch(ref) {

        case IS_CENTRE:
            return (- getLength() / 2.0);

        case IS_EXIT:
            return (- getLength());

        default:
            return 0.0;
    }
}


double Element::getExit(ReferenceType ref) const {
    switch(ref) {

        case IS_ENTRY:
            return getLength();

        case IS_CENTRE:
            return (getLength() / 2.0);

        default:
            return 0.0;
    }
}


void Element::setShared(bool flag) {
    Object::setShared(flag);
    if(flag) itsClassicElement->makeSharable();
}


Element::Element(int size, const char *name, const char *help):
    Object(size, name, help)
{}


Element::Element(const std::string &name, Element *parent):
    Object(name, parent)
{}
