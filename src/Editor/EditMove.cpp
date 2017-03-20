// ------------------------------------------------------------------------
// $RCSfile: EditMove.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditMove
//   The class for the OPAL sequence editor MOVE command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditMove.h"
#include "AbstractObjects/PlaceRep.h"
#include "Attributes/Attributes.h"
#include "Editor/Edit.h"
#include "Utilities/Options.h"
#include <iostream>


// Class EditMove
// ------------------------------------------------------------------------

// The attributes of class EditMove.

namespace {
    enum {
        SELECTED,  // If true, move selected elements.
        ELEMENT,   // The element to be moved.
        BY,        // The shift amount.
        TO,        // The new position.
        FROM,      // An origin name.
        SIZE
    };
}


EditMove::EditMove():
    Editor(SIZE, "MOVE",
           "The \"MOVE\" sub-command moves element(s) in the sequence "
           "being edited.") {
    itsAttr[SELECTED] = Attributes::makeBool
                        ("SELECTED", "If true, all selected elements are moved");
    itsAttr[ELEMENT] = Attributes::makePlace
                       ("ELEMENT", "Name of single element to be moved");
    itsAttr[BY] = Attributes::makeReal
                  ("BY", "Amount in m by which elements should be moved");
    itsAttr[TO] = Attributes::makeReal
                  ("TO", "New position in m relative to origin");
    itsAttr[FROM] = Attributes::makePlace
                    ("FROM", "Name of element defining the origin (default is start)");

    registerOwnership(AttributeHandler::SUB_COMMAND);
}


EditMove::EditMove(const std::string &name, EditMove *parent):
    Editor(name, parent)
{}


EditMove::~EditMove()
{}


EditMove *EditMove::clone(const std::string &name) {
    return new EditMove(name, this);
}


void EditMove::execute() {
    const PlaceRep cls  = Attributes::getPlace(itsAttr[ELEMENT]);
    int count = 0;

    if(Attributes::getBool(itsAttr[SELECTED])) {
        double by = Attributes::getReal(itsAttr[BY]);
        count = Edit::block->moveMultiple(by);
    } else if(itsAttr[BY]) {
        double by = Attributes::getReal(itsAttr[BY]);
        count = Edit::block->moveSingleRel(cls, cls, by);
    } else if(itsAttr[FROM]) {
        const PlaceRep from = Attributes::getPlace(itsAttr[FROM]);
        double to = Attributes::getReal(itsAttr[TO]);
        count = Edit::block->moveSingleRel(cls, from, to);
    } else {
        double to = Attributes::getReal(itsAttr[TO]);
        count = Edit::block->moveSingleAbs(cls, to);
    }

    if(Options::info) {
        if(count == 0) {
            std::cerr << "\nNo element";
        } else if(count == 1) {
            std::cerr << "\n1 element";
        } else {
            std::cerr << '\n' << count << " elements";
        }

        std::cerr << " moved.\n" << std::endl;
    }
}