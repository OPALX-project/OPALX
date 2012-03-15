// ------------------------------------------------------------------------
// $RCSfile: RealConstant.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RealConstant
//   Implements a REAL_CONSTANT definition.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:26:42 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "ValueDefinitions/RealConstant.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/Options.h"
#include <cmath>
#include <iostream>

static const double MATH_E  = exp(1.0);
static const double MATH_PI = 4.0 * atan(1.0);


// Class RealConstant
// ------------------------------------------------------------------------

RealConstant::RealConstant():
    ValueDefinition(1, "REAL_CONSTANT",
                    "The \"REAL CONSTANT\" statement defines a global "
                    "real constant:\n"
                    "\tREAL CONSTANT <name> = <real-expression>;\n") {
    itsAttr[0] = Attributes::makeReal("VALUE", "The constant value", 0.0);

    // Define the standard constants.
    OPAL.create(new RealConstant("PI",     this, MATH_PI));
    OPAL.create(new RealConstant("TWOPI",  this, 2.0 * MATH_PI));
    OPAL.create(new RealConstant("RADDEG", this, MATH_PI / 180.0));
    OPAL.create(new RealConstant("DEGRAD", this, 180.0 / MATH_PI));
    OPAL.create(new RealConstant("E",      this, MATH_E));
    OPAL.create(new RealConstant("EMASS",  this, 0.51099906e-3));
    OPAL.create(new RealConstant("PMASS",  this, 0.93827231));
    OPAL.create(new RealConstant("CLIGHT", this, 299792458.));
}


RealConstant::RealConstant(const string &name, RealConstant *parent):
    ValueDefinition(name, parent)
{}


RealConstant::RealConstant(const string &name, RealConstant *parent,
                           double value):
    ValueDefinition(name, parent) {
    Attributes::setReal(itsAttr[0], value);
    itsAttr[0].setReadOnly(true);
    builtin = true;
}


RealConstant::~RealConstant()
{}


bool RealConstant::canReplaceBy(Object *) {
    return false;
}


RealConstant *RealConstant::clone(const string &name) {
    return new RealConstant(name, this);
}


void RealConstant::doomGet(const DoomReader &reader) {
    itsAttr[0].doomGet(reader, 0);
}


void RealConstant::doomPut(DoomWriter &writer) const {
    itsAttr[0].doomPut(writer, 0);
}


double RealConstant::getReal() const {
    return Attributes::getReal(itsAttr[0]);
}


//JMJ 30/11/2000 added harmless semi-colon in OPAL8 output
void RealConstant::print(std::ostream &os) const {
    if(Options::opal8) {
        os << getOpalName() << ":CONSTANT=" << itsAttr[0] << ';';
    } else {
        os << "REAL CONST " << getOpalName() << '=' << itsAttr[0] << ';';
    }
    os << std::endl;
}
