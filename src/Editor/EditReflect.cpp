// ------------------------------------------------------------------------
// $RCSfile: EditReflect.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditReflect
//   The class for the OPAL sequence editor REFLECT command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditReflect.h"
#include "Editor/Edit.h"
#include "Utilities/Options.h"
#include <iostream>


// Class EditReflect
// ------------------------------------------------------------------------

EditReflect::EditReflect():
    Editor(0, "REFLECT",
           "The \"REFLECT\" inverts the order of all sequence positions.")
{}


EditReflect::EditReflect(const string &name, EditReflect *parent):
    Editor(name, parent)
{}


EditReflect::~EditReflect()
{}


EditReflect *EditReflect::clone(const string &name) {
    return new EditReflect(name, this);
}


void EditReflect::execute() {
    Edit::block->reflect();

    if(Options::info) {
        std::cerr << "\nSequence has been reflected.\n" << std::endl;
    }
}
