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
#include "AbstractObjects/DoomWriter.h"
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
}


BoolConstant::BoolConstant(const string &name, BoolConstant *parent):
    ValueDefinition(name, parent)
{}


BoolConstant::~BoolConstant()
{}


bool BoolConstant::canReplaceBy(Object *) {
    return false;
}


BoolConstant *BoolConstant::clone(const string &name) {
    return new BoolConstant(name, this);
}


void BoolConstant::doomGet(const DoomReader &reader) {
    itsAttr[0].doomGet(reader, 0);
}


void BoolConstant::doomPut(DoomWriter &writer) const {
    itsAttr[0].doomPut(writer, 0);
}


bool BoolConstant::getBool() const {
    return Attributes::getBool(itsAttr[0]);
}


//JMJ 30/11/2000 added harmless semi-colon in OPAL8 output
void BoolConstant::print(std::ostream &os) const {
    if(Options::opal8) {
        os << getOpalName() << '=' << (Attributes::getBool(itsAttr[0]) ? 1 : 0) << ';';
    } else {
        os << "BOOL CONST " << getOpalName() << '=' << itsAttr[0] << ';';
    }
    os << std::endl;
}
