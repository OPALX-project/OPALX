// ------------------------------------------------------------------------
// $RCSfile: EditFlatten.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditFlatten
//   The class for the OPAL sequence editor FLATTEN command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditFlatten.h"
#include "Editor/Edit.h"
#include "Lines/Sequence.h"
#include "Utilities/Options.h"
#include <iostream>


// Class EditFlatten
// ------------------------------------------------------------------------


EditFlatten::EditFlatten():
    Editor(0, "FLATTEN",
           "The \"FLATTEN\" sub-command creates a flat sequence from the "
           "sequence being edited.")
{}


EditFlatten::EditFlatten(const std::string &name, EditFlatten *parent):
    Editor(name, parent)
{}


EditFlatten::~EditFlatten()
{}


EditFlatten *EditFlatten::clone(const std::string &name) {
    return new EditFlatten(name, this);
}


void EditFlatten::execute() {
    Edit::block->flatten();

    if(Options::info) {
        std::cerr << "\nSequence has been flattened.\n" << std::endl;
    }
}
