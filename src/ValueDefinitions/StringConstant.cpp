//
// Class StringConstant
//   The STRING CONSTANT definition.
//
// Copyright (c) 2000 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "ValueDefinitions/StringConstant.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/Util.h"

#include <iostream>


StringConstant::StringConstant():
    ValueDefinition(1, "STRING_CONSTANT",
                    "The \"STRING CONSTANT\" statement defines a global "
                    "string constant:\n"
                    "\tSTRING CONSTANT <name> = <String-expression>;\n") {
    itsAttr[0] = Attributes::makeString("VALUE", "The constant value");

    registerOwnership(AttributeHandler::STATEMENT);

    OpalData *opal = OpalData::getInstance();
    opal->create(new StringConstant("GITREVISION", this, Util::getGitRevision()));
}


StringConstant::StringConstant(const std::string &name, StringConstant *parent):
    ValueDefinition(name, parent)
{}


StringConstant::StringConstant(const std::string &name, StringConstant *parent, const std::string &value):
    ValueDefinition(name, parent)
{
    Attributes::setString(itsAttr[0], value);
    itsAttr[0].setReadOnly(true);
    builtin = true;
}


StringConstant::~StringConstant()
{}


bool StringConstant::canReplaceBy(Object *) {
    return false;
}


StringConstant *StringConstant::clone(const std::string &name) {
    return new StringConstant(name, this);
}


std::string StringConstant::getString() const {
    return Attributes::getString(itsAttr[0]);
}


void StringConstant::print(std::ostream &os) const {
    os << "STRING " << getOpalName() << '=' << itsAttr[0] << ';';
    os << std::endl;
}

void StringConstant::printValue(std::ostream &os) const {
    os << itsAttr[0];
}
