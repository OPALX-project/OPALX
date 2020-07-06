//
// Class ElementFactory
//   Concrete factory class for CLASSIC elements.
//   Defines the beamline element creation via the factory pattern.
//   When the factory is constructed, empty elements are first created and
//   stored.
//   {p}
//   With the makeElement() method, these elements can be cloned.
//   The factory also implements an
//   element repository which can store beam lines.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Construction/ElementFactory.h"
#include "AbsBeamline/AttributeSet.h"
#include "AbsBeamline/ElementBase.h"


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
