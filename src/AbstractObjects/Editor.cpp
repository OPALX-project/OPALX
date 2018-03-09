// ------------------------------------------------------------------------
// $RCSfile: Editor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Editor
//   The base class for all OPAL sequence editor commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:35 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Editor.h"


// Class Editor
// ------------------------------------------------------------------------

Editor::~Editor()
{}


const std::string Editor::getCategory() const {
    return "EDITOR";
}


bool Editor::shouldTrace() const {
    return true;

}

bool Editor::shouldUpdate() const {
    return false;
}


Editor::Editor(const std::string &name, Editor *parent):
    Object(name, parent)
{}


Editor::Editor(int size, const char *name, const char *help):
    Object(size, name, help)
{}
