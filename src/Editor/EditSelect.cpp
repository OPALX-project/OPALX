// ------------------------------------------------------------------------
// $RCSfile: EditSelect.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditSelect
//   The class for the OPAL sequence editor SELECT command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditSelect.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/PlaceRep.h"
#include "AbstractObjects/RangeRep.h"
#include "Attributes/Attributes.h"
#include "Editor/Edit.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <iostream>


// Class EditSelect
// ------------------------------------------------------------------------

// The attributes of class EditSelect.
namespace {
    enum {
        RANGE,     // The range to be considered.
        CLASS,     // The class to be considered.
        TYPE,      // The element TYPE keyword.
        PATTERN,   // The regular expression to be matched.
        FULL,      // Flag for full selection.
        CLEAR,     // Flag for full de-selection.
        SIZE
    };
}


EditSelect::EditSelect():
    Editor(SIZE, "SELECT",
           "The \"SELECT\" sub-command marks element(s) in the sequence "
           "to be affected by subsequent edit sub-commands") {
    itsAttr[RANGE] = Attributes::makeRange
                     ("RANGE", "Range to be considered for selection");
    itsAttr[CLASS] = Attributes::makeString
                     ("CLASS", "Class of elements to be selected within range");
    itsAttr[TYPE] = Attributes::makeString
                    ("TYPE", "The element \"TYPE\" attribute.");
    itsAttr[PATTERN] = Attributes::makeString
                       ("PATTERN", "Regular expression for name matching");
    itsAttr[FULL] = Attributes::makeBool
                    ("FULL", "If true, all element are selected");
    itsAttr[CLEAR] = Attributes::makeBool
                     ("CLEAR", "If true, all selections are cleared");
}


EditSelect::EditSelect(const std::string &name, EditSelect *parent):
    Editor(name, parent)
{}


EditSelect::~EditSelect()
{}


EditSelect *EditSelect::clone(const std::string &name) {
    return new EditSelect(name, this);
}


void EditSelect::execute() {
    if(itsAttr[FULL] && Attributes::getBool(itsAttr[FULL])) {
        Edit::block->selectFull();

        if(Options::info) {
            std::cerr << "\nAll elements selected.\n";
        }
    } else if(itsAttr[CLEAR] && Attributes::getBool(itsAttr[CLEAR])) {
        Edit::block->selectClear();

        if(Options::info) {
            std::cerr << "\nAll elements de-selected.\n";
        }
    } else {
        // Now set select flags according to range, class and pattern.
        const RangeRep &range = Attributes::getRange(itsAttr[RANGE]);
        const std::string &clsName = Attributes::getString(itsAttr[CLASS]);
        const std::string &typName = Attributes::getString(itsAttr[TYPE]);
        const std::string &pattern = Attributes::getString(itsAttr[PATTERN]);
        int count = Edit::block->select(range, clsName, typName, pattern);

        // inform user about number of selections done
        if(Options::info) {
            if(count == 0) {
                std::cerr << "\nNo elements";
            } else if(count == 1) {
                std::cerr << "\n1 element";
            } else {
                std::cerr << '\n' << count << " elements";
            }

            std::cerr << " selected.\n";
        }
    }

    std::cerr << std::endl;
}
