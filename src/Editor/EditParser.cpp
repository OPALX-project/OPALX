// ------------------------------------------------------------------------
// $RCSfile: EditParser.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditParser
//   The parser class for the OPAL sequence editor.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditParser.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "Editor/EditCycle.h"
#include "Editor/EditEnd.h"
#include "Editor/EditFlatten.h"
#include "Editor/EditInstall.h"
#include "Editor/EditMove.h"
#include "Editor/EditReflect.h"
#include "Editor/EditRemove.h"
#include "Editor/EditReplace.h"
#include "Editor/EditSelect.h"
#include "MemoryManagement/Pointer.h"
#include "Parser/Statement.h"
#include "Utilities/ParseError.h"
#include <ctime>
#include <iostream>


// Class EditParser
// ------------------------------------------------------------------------

EditParser::EditParser()
{
  editDirectory.insert("CYCLE",   new EditCycle());
  editDirectory.insert("ENDEDIT", new EditEnd());
  editDirectory.insert("FLATTEN", new EditFlatten());
  editDirectory.insert("MOVE",    new EditMove());
  editDirectory.insert("REFLECT", new EditReflect());
  editDirectory.insert("REMOVE",  new EditRemove());
  editDirectory.insert("REPLACE", new EditReplace());
  editDirectory.insert("SELECT",  new EditSelect());
}


EditParser::~EditParser()
{}


Object *EditParser::find(const string &name) const
{
  return editDirectory.find(name);
}


void EditParser::parse(Statement &stat) const
{ 
  string name =
    Expressions::parseString(stat, "Identifier or keyword expected.");

  if (stat.delimiter('?')) {
    // "<class>?": give help for class.
    printHelp(name);
  } else if (stat.delimiter(':')) {
    // "<name>:...": INSTALL command.
    parseInstall(stat);
  } else if (stat.delimiter(',')) {
    // "<name>,...": executable or INSTALL command.
    if (find(name)) {
      parseAction(stat);
    } else {
      parseInstall(stat);
    }
  } else if (stat.delimiter(';') || stat.atEnd()) {  
    // "<name>;": Executable command.
    parseAction(stat);
  } else {
    throw ParseError("EditParser::parse()",
		     "Only sequence editor commands are accepted here.");
  }
}


void EditParser::parseInstall(Statement &statement) const
{
  // This command must be an "INSTALL" sub-command.
  statement.start();
  Object *copy = 0;
  try {
    copy = new EditInstall();
    copy->parse(statement);
    parseEnd(statement);
    execute(copy, "INSTALL");
  } catch (...) {
    delete copy;
    throw;
  }
}
