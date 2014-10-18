// ------------------------------------------------------------------------
// $RCSfile: ErrorSelect.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorSelect
//   The class for the OPAL sequence error SELECT command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/ErrorSelect.h"
#include "Algorithms/Flagger.h"
#include "Attributes/Attributes.h"
#include "Errors/Error.h"
#include "Tables/Selector.h"
#include "Utilities/Options.h"
#include <iostream>


// Class ErrorSelect
// ------------------------------------------------------------------------

// The attributes of class ErrorSelect.
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


ErrorSelect::ErrorSelect():
    Action(SIZE, "SELECT",
           "The \"SELECT\" sub-command marks element(s) in the sequence "
           "to be affected by subsequent error sub-commands") {
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


ErrorSelect::ErrorSelect(const std::string &name, ErrorSelect *parent):
    Action(name, parent)
{}


ErrorSelect::~ErrorSelect()
{}


ErrorSelect *ErrorSelect::clone(const std::string &name) {
    return new ErrorSelect(name, this);
}


void ErrorSelect::execute() {
    if(Attributes::getBool(itsAttr[FULL])) {
        // Select all positions.
        Flagger flagger(*Error::block->itsLine, true);
        flagger.execute();

        if(Options::info) {
            std::cerr << "\nAll elements selected.\n";
        }
    } else if(Attributes::getBool(itsAttr[CLEAR])) {
        // Deselect all selections.
        Flagger flagger(*Error::block->itsLine, false);
        flagger.execute();

        if(Options::info) {
            std::cerr << "\nAll elements de-selected.\n";
        }
    } else {
        Selector sel(*Error::block->itsLine,
                     Attributes::getRange(itsAttr[RANGE]),
                     Attributes::getString(itsAttr[CLASS]),
                     Attributes::getString(itsAttr[TYPE]),
                     Attributes::getString(itsAttr[PATTERN]));
        sel.execute();

        if(Options::info) {
            int count = sel.getCount();

            if(count == 0) {
                std::cerr << "\nNo elements";
            } else if(count == 1) {
                std::cerr << "\n1 element";
            } else {
                std::cerr << '\n' << count << " elements";
            }

            std::cerr << "  selected.\n";
        }
    }

    std::cerr << std::endl;
}
