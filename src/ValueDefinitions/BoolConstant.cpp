// ------------------------------------------------------------------------
// $RCSfile: BoolConstant.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BoolConstant
//   Implements a BOOL_CONSTANT definition.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:26:42 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "ValueDefinitions/BoolConstant.h"
#include "Attributes/Attributes.h"
#include "Utilities/Options.h"
#include <iostream>


// Class BoolConstant
// ------------------------------------------------------------------------

BoolConstant::BoolConstant():
    ValueDefinition(1, "BOOL_CONSTANT",
                    "The \"BOOL CONSTANT\" statement defines a global "
                    "logical constant:\n"
                    "\tBOOL CONSTANT <name> = <Bool-expression>;\n") {
    itsAttr[0] = Attributes::makeBool("VALUE", "The constant value");

    registerOwnership(AttributeHandler::STATEMENT);
}


BoolConstant::BoolConstant(const std::string &name, BoolConstant *parent):
    ValueDefinition(name, parent)
{}


BoolConstant::~BoolConstant()
{}


bool BoolConstant::canReplaceBy(Object *) {
    return false;
}


BoolConstant *BoolConstant::clone(const std::string &name) {
    return new BoolConstant(name, this);
}


bool BoolConstant::getBool() const {
    return Attributes::getBool(itsAttr[0]);
}


void BoolConstant::print(std::ostream &os) const {
    os << "BOOL CONST " << getOpalName() << '=' << itsAttr[0] << ';';
    os << std::endl;
}

void BoolConstant::printValue(std::ostream &os) const {
    os << itsAttr[0];
}