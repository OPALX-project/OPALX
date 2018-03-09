// ------------------------------------------------------------------------
// $RCSfile: ElementFactory.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ElementFactory
//   This class is the main CLASSIC element factory.
//
// ------------------------------------------------------------------------
// Class category: Construction
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Construction/ElementFactory.h"
#include "AbsBeamline/AttributeSet.h"
#include "AbsBeamline/ElementBase.h"


// Class ElementFactory
// ------------------------------------------------------------------------

ElementFactory::ElementFactory():
    inventory()
{}


ElementFactory::~ElementFactory()
{}


bool ElementFactory::define(ElementBase *newElement) {
    std::string name = newElement->getName();
    MapType::value_type value(name, newElement);
    std::pair<MapType::iterator, bool> index = inventory.insert(value);

    if(index.second) {
        // Insertion took place (no duplicate).
        return true;
    } else {
        // Insertion rejected.
        delete newElement;
        return false;
    }
}


void ElementFactory::erase(const std::string &name) {
    inventory.erase(name);
}


ElementBase *ElementFactory::find(const std::string &name) const {
    MapType::const_iterator index = inventory.find(name);

    if(index == inventory.end()) {
        return 0;
    } else {
        return (*index).second;
    }
}


ElementBase *ElementFactory::makeElement(const std::string &type,
        const std::string &name,
        const AttributeSet &set) {
    ElementBase *model = find(type);
    ElementBase *copy  = 0;

    try {
        if(model != 0) {
            ElementBase *copy = model->clone();
            copy->setName(name);
            copy->update(set);
            storeElement(copy);
            return copy;
        }
    } catch(...) {
        delete copy;
    }

    return 0;
}


bool ElementFactory::storeElement(ElementBase *newElement) {
    std::string name = newElement->getName();
    MapType::value_type value(name, newElement);
    std::pair<MapType::iterator, bool> index = inventory.insert(value);

    if(index.second) {
        // Insertion took place (no duplicate).
        return true;
    } else {
        // Insertion was a replacement.
        (*index.first).second = newElement;
        return false;
    }
}
