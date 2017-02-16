// ------------------------------------------------------------------------
// $RCSfile: EditEnd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditEnd
//   The class for the OPAL sequence editor ENDEDIT command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditEnd.h"
#include "Attributes/Attributes.h"
#include "Editor/Edit.h"


// Class EditEnd
// ------------------------------------------------------------------------


EditEnd::EditEnd():
    Editor(1, "ENDEDIT",
           "The \"ENDEDIT\" sub-command terminates sequence editing mode.") {
    itsAttr[0] = Attributes::makeString
                 ("NAME", "New name for the edited sequence (default = old name)");

    registerOwnership(AttributeHandler::SUB_COMMAND);
}


EditEnd::EditEnd(const std::string &name, EditEnd *parent):
    Editor(name, parent)
{}


EditEnd::~EditEnd()
{}


EditEnd *EditEnd::clone(const std::string &name) {
    return new EditEnd(name, this);
}


void EditEnd::execute() {
    try {
        std::string newName = Edit::block->itsSequence->getOpalName();
        if(itsAttr[0]) newName = Attributes::getString(itsAttr[0]);
        Edit::block->finish(newName);
        Edit::block->parser.stop();
    } catch(...) {
        Edit::block->parser.stop();
        throw;
    }
}