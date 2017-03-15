// ------------------------------------------------------------------------
// $RCSfile: Title.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Title
//   The class for the OPAL Title command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "BasicActions/Title.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"


// Class Title
// ------------------------------------------------------------------------

Title::Title():
    Action(1, "TITLE",
           "The \"TITLE\" statement defines a new page title for subsequent "
           "output.") {
    itsAttr[0] = Attributes::makeString
                 ("STRING", "The title string");

    registerOwnership(AttributeHandler::STATEMENT);
}


Title::Title(const std::string &name, Title *parent):
    Action(name, parent)
{}


Title::~Title()
{}


Title *Title::clone(const std::string &name) {
    return new Title(name, this);
}


void Title::execute() {
    std::string title = Attributes::getString(itsAttr[0]);
    OpalData::getInstance()->storeTitle(title);
}


void Title::parse(Statement &statement) {
    parseShortcut(statement);
}