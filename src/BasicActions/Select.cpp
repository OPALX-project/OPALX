//
// Class Select
//   The class for OPAL SELECT command.
//
// Copyright (c) 2000 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "BasicActions/Select.h"
#include "Attributes/Attributes.h"
#include "Utility/Inform.h"

#include <iostream>

extern Inform* gmsg;

namespace {
    enum {
        LINE,     // The line to be affected.
        FULL,     // If true, all elements are selected.
        CLEAR,    // If true, all selections are cleared.
        RANGE,    // The range to be considered.
        CLASS,    // The class of elements to be selected.
        TYPE,     // The type name of elements to be selected.
        PATTERN,  // The regular expression for matching names.
        SIZE
    };
}

Select::Select()
    : Action(SIZE, "SELECT",
             "The \"SELECT\" sub-command selects the positions to be affected "
             "by subsequent error sub-commands.") {
    itsAttr[LINE] = Attributes::makeString(
            "LINE", "Name of the lattice to be affected by selections", "UNNAMED_USE");
    itsAttr[FULL]  = Attributes::makeBool("FULL", "If true, all elements are selected");
    itsAttr[CLEAR] = Attributes::makeBool("CLEAR", "If true, all selections are cleared");
    itsAttr[RANGE] = Attributes::makeRange(
            "RANGE", "Range to be considered for selection (default: full range)");
    itsAttr[CLASS] =
            Attributes::makeString("CLASS", "Name of class to be selected (default: all classes)");
    itsAttr[TYPE] = Attributes::makeString(
            "TYPE", "The type name of elements to be selected (default: all types)");
    itsAttr[PATTERN] = Attributes::makeString(
            "PATTERN", "Regular expression for matching names (default: all names)");

    registerOwnership(AttributeHandler::SUB_COMMAND);
}

Select::~Select() {}

void Select::execute() {
    // The selection flags SELECT set on line members were read by nothing.
    *gmsg << level1 << "* SELECT has no effect in OPALX." << endl;
}
