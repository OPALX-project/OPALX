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
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/Options.h"
#include "Physics/Physics.h"

#include <cmath>
#include <iostream>

// Class RealConstant
// ------------------------------------------------------------------------

RealConstant::RealConstant():
    ValueDefinition(1, "REAL_CONSTANT",
                    "The \"REAL CONSTANT\" statement defines a global "
                    "real constant:\n"
                    "\tREAL CONSTANT <name> = <real-expression>;\n") {
    itsAttr[0] = Attributes::makeReal("VALUE", "The constant value", 0.0);

    // Define the standard constants.
    OpalData *OPAL = OpalData::getInstance();
    OPAL->create(new RealConstant("PI",     this, Physics::pi));
    OPAL->create(new RealConstant("TWOPI",  this, Physics::two_pi));
    OPAL->create(new RealConstant("RADDEG", this, 180.0 / Physics::pi / 180.0));
    OPAL->create(new RealConstant("DEGRAD", this, Physics::pi / 180.0));
    OPAL->create(new RealConstant("E",      this, Physics::e));

    OPAL->create(new RealConstant("EMASS",  this, Physics::m_e));
    OPAL->create(new RealConstant("PMASS",  this, Physics::m_p));
    OPAL->create(new RealConstant("HMMASS", this, Physics::m_hm));
    OPAL->create(new RealConstant("UMASS", this, Physics::m_u));
    OPAL->create(new RealConstant("CMASS", this, Physics::m_c));
    OPAL->create(new RealConstant("MMASS", this, Physics::m_mu));
    OPAL->create(new RealConstant("DMASS", this, Physics::m_d));
    OPAL->create(new RealConstant("XEMASS", this, Physics::m_xe));

    OPAL->create(new RealConstant("CLIGHT", this, Physics::c));

    OPAL->create(new RealConstant("OPALVERSION", this, 120));
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


double RealConstant::getReal() const {
    return Attributes::getReal(itsAttr[0]);
}


void RealConstant::print(std::ostream &os) const {
    os << "REAL CONST " << getOpalName() << '=' << itsAttr[0] << ';';
    os << std::endl;
}
