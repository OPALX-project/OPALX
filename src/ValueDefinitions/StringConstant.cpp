// ------------------------------------------------------------------------
// $RCSfile: StringConstant.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: StringConstant
//   Implements a OPAL STRING_CONSTANT definition.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:26:42 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "ValueDefinitions/StringConstant.h"
#include "AbstractObjects/DoomWriter.h"
#include "Attributes/Attributes.h"
#include "Utilities/Options.h"
#include <iostream>


// Class StringConstant
// ------------------------------------------------------------------------

StringConstant::StringConstant():
    ValueDefinition(1, "STRING_CONSTANT",
                    "The \"STRING CONSTANT\" statement defines a global "
                    "string constant:\n"
                    "\tSTRING CONSTANT <name> = <String-expression>;\n") {
    itsAttr[0] = Attributes::makeString("VALUE", "The constant value");
}


StringConstant::StringConstant(const string &name, StringConstant *parent):
    ValueDefinition(name, parent)
{}


StringConstant::~StringConstant()
{}


bool StringConstant::canReplaceBy(Object *) {
    return false;
}


StringConstant *StringConstant::clone(const string &name) {
    return new StringConstant(name, this);
}


void StringConstant::doomGet(const DoomReader &reader) {
    itsAttr[0].doomGet(reader, 0);
}


void StringConstant::doomPut(DoomWriter &writer) const {
    itsAttr[0].doomPut(writer, 0);
}


string StringConstant::getString() const {
    return Attributes::getString(itsAttr[0]);
}


//JMJ 30/11/2000 added harmless semi-colon in OPAL8 output
void StringConstant::print(std::ostream &os) const {
    if(Options::opal8) {
        os << getOpalName() << ":STRING,\"" << itsAttr[0] << '"' << ';';
    } else {
        os << "STRING " << getOpalName() << '=' << itsAttr[0] << ';';
    }
    os << std::endl;
}
