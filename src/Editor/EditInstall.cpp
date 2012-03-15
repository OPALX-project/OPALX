// ------------------------------------------------------------------------
// $RCSfile: EditInstall.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditInstall
//   The class for the OPAL sequence editor INSTALL command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditInstall.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/PlaceRep.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Editor/Edit.h"
#include "Parser/Statement.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"
#include "Utilities/Options.h"
#include <iostream>


// Class EditInstall
// ------------------------------------------------------------------------

// The attributes of class EditInstall.
namespace {
  enum {
    AT,        // The position.
    FROM,      // The place defining one or more positions.
    SIZE
  };
}


EditInstall::EditInstall():
  Editor(SIZE, "INSTALL",
	 "The \"INSTALL\" sub-command installs new elements in the sequence.")
{
  itsAttr[AT] = Attributes::makeReal
    ("AT",
     "Position in metres for installation relative to origin(s)");
  itsAttr[FROM] = Attributes::makePlace
    ("FROM", "Position of element defining origin (default is start)");
}


EditInstall::EditInstall(const string &name, EditInstall *parent):
  Editor(name, parent)
{}


EditInstall::~EditInstall()
{}


EditInstall *EditInstall::clone(const string &name)
{
  return new EditInstall(name, this);
}


void EditInstall::execute()
{
  if (! itsAttr[AT]) {
    throw OpalException("EditInstall::execute()",
		       "The \"AT\" attribute must be present.");
  }

  double at = Attributes::getReal(itsAttr[AT]);
  ElementBase *elm = newElement->getElement();
  const PlaceRep from = Attributes::getPlace(itsAttr[FROM]);
  int count;

  // Find where to install.
  if (from.isSelected()) {
    count = Edit::block->installMultiple(elm, at);
  } else {
    count = Edit::block->installSingle(from, elm, at);
  }
  
  if (Options::info) {
    const string &name = elm->getName();

    if (count == 0) {
      std::cerr << "\nNo \"" << name << '"';
    } else if (count == 1) {
      std::cerr << "\n1 \"" << name << '"';
    } else {
      std::cerr << '\n' << count << " \"" << name << "\"'s";
    }

    std::cerr << " installed.\n" << std::endl;
  }
}


void EditInstall::parse(Statement &stat)
{
  // Read object identifier.
  string className;
  string objectName = Expressions::parseString(stat, "Object name expected.");

  // Read class identifier.
  if (stat.delimiter(':')) {
    className = Expressions::parseString(stat, "Class name expected.");
  } else {
    className = objectName;
    objectName = "";
  }

  // Find exemplar object.
  if (Object *object = OPAL.find(className)) {
    // Instantiate, if necessary.
    int defined = 0;
    Pointer<Object> copy;

    if (stat.delimiter('(')) {
      if (objectName.empty()) {
	throw ParseError("EditInstall::parse()",
			 "An instantiation of \"" + className +
			 "\" in a sequence shoud have name.");
      }
      copy = object->makeInstance(objectName, stat, 0);
      defined = 1;
    } else if (! objectName.empty()) {
      copy = object->clone(objectName);
      defined = 2;
    } else {
      copy = object;
    }

    // Check that there is no invalid redefinition.
    if (defined) {
      if (! objectName.empty()  &&  OPAL.find(objectName)) {
	throw ParseError("EditInstall::parse()",
			 "You cannot redefine \"" + objectName +
			 "\" within the sequence editor.");
      }
    }

    if (Element *element = dynamic_cast<Element *>(&*copy)) {
      newElement = element;
      
      while (stat.delimiter(',')) {
	string attrName =
	  Expressions::parseString(stat, "Attribute name expected.");
	Attribute *attr = 0;
	  
	if ((attr = findAttribute(attrName)) != 0) {
	  // An attribute of the "INSTALL" command.
	  Expressions::parseDelimiter(stat, '=');
	  attr->parse(stat, true);
	} else if (defined) {
	  // An attribute of the new element.
	  if ((attr = newElement->findAttribute(attrName)) != 0) {
	    if (stat.delimiter('=')) {
	      attr->parse(stat, true);
	    } else if (stat.delimiter("=")) {
	      attr->parse(stat, false);
	    } else {
	      throw ParseError("EditInstall::parsePosition()",
			       "Delimiter \"=\" or \":=\" expected.");
	    }
	  } else {
	    throw ParseError("EditInstall::parsePosition()",
			     "Element \"" + newElement->getOpalName() +
			     "\" has no attribute \"" + attrName + "\".");
	  }
	} else {
	  throw ParseError("EditInstall::parsePosition()",
			   "Overriding attributes not permitted here.");
	}
      }

      if (defined) OPAL.define(&*copy);
    } else {
      throw ParseError("EditInstall::parse()", "Object \"" +
		       objectName + "\" is not an element.");
    }
  } else {
    throw ParseError("EditInstall::parse()",
		     "Element class name \"" + className + "\" is unknown.");
  }
}
