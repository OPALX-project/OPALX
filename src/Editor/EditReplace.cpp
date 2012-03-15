// ------------------------------------------------------------------------
// $RCSfile: EditReplace.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditReplace
//   The class for the OPAL sequence editor REPLACE command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditReplace.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/PlaceRep.h"
#include "Attributes/Attributes.h"
#include "Editor/Edit.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <iostream>


// Class EditReplace
// ------------------------------------------------------------------------

// The attributes of class EditReplace.
namespace {
    enum {
        SELECTED,  // If true, replace selected elements.
        CLASS,     // The class to be replaced.
        BY,        // The replacement name.
        SIZE
    };
}


EditReplace::EditReplace():
    Editor(SIZE, "REPLACE",
           "The \"REPLACE\" sub-command replaces element(s) in the sequence "
           "being edited.") {
    itsAttr[SELECTED] = Attributes::makeBool
                        ("SELECTED", "If true, all selected elements are replaced");
    itsAttr[CLASS] = Attributes::makePlace
                     ("CLASS", "Name of element class to be replaced");
    itsAttr[BY] = Attributes::makeString("BY", "name of replacement class");
}


EditReplace::EditReplace(const string &name, EditReplace *parent):
    Editor(name, parent)
{}


EditReplace::~EditReplace()
{}


EditReplace *EditReplace::clone(const string &name) {
    return new EditReplace(name, this);
}


void EditReplace::execute() {
    int count = 0;
    string newName = Attributes::getString(itsAttr[BY]);

    // Check consistency.
    if(! bool(itsAttr[CLASS])  ||  newName.empty()) {
        throw OpalException("EditReplace::execute()",
                            "\"ELEMENT=SELECTED\" or \"ELEMENT=<name>\" and "
                            "\"BY=<name>\" are required for \"REPLACE\".");
    }

    // Look up new name.
    ElementBase *elem = Element::find(newName)->getElement();

    if(Attributes::getBool(itsAttr[SELECTED])) {
        count = Edit::block->replaceMultiple(elem);
    } else {
        const PlaceRep pos = Attributes::getPlace(itsAttr[CLASS]);
        count = Edit::block->replaceSingle(pos, elem);
    }

    if(Options::info) {
        if(count == 0) {
            std::cerr << "\nNo elements";
        } else if(count == 1) {
            std::cerr << "\n1 element";
        } else {
            std::cerr << '\n' << count << " elements";
        }

        std::cerr << " replaced.\n" << std::endl;
    }
}
