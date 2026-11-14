//
// Class RealConstant
//   The REAL CONSTANT definition.
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
#include "ValueDefinitions/RealConstant.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BuildInfo.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/Options.h"
#include "Utility/IpplInfo.h"

#include <cmath>
#include <iostream>

RealConstant::RealConstant()
    : ValueDefinition(
              1, "REAL_CONSTANT",
              "The \"REAL CONSTANT\" statement defines a global "
              "real constant:\n"
              "\tREAL CONSTANT <name> = <real-expression>;\n") {
    itsAttr[0] = Attributes::makeReal("VALUE", "The constant value", 0.0);

    registerOwnership(AttributeHandler::STATEMENT);
}

void RealConstant::defineBuiltins() {
    OpalData* opal = OpalData::getInstance();
    opal->create(new RealConstant("PI", Physics::pi));
    opal->create(new RealConstant("TWOPI", Physics::two_pi));
    opal->create(new RealConstant("RADDEG", Units::rad2deg));
    opal->create(new RealConstant("DEGRAD", Units::deg2rad));
    opal->create(new RealConstant("E", Physics::e));

    opal->create(new RealConstant("MUMASS", Physics::m_mu));
    opal->create(new RealConstant("EMASS", Physics::m_e));
    opal->create(new RealConstant("PMASS", Physics::m_p));
    opal->create(new RealConstant("HMMASS", Physics::m_hm));
    opal->create(new RealConstant("H2PMASS", Physics::m_h2p));
    opal->create(new RealConstant("DMASS", Physics::m_d));
    opal->create(new RealConstant("ALPHAMASS", Physics::m_alpha));
    opal->create(new RealConstant("CMASS", Physics::m_c));
    opal->create(new RealConstant("XEMASS", Physics::m_xe));
    opal->create(new RealConstant("UMASS", Physics::m_u));

    opal->create(new RealConstant("CLIGHT", Physics::c));
    opal->create(new RealConstant("AMU", Physics::amu));

    opal->create(new RealConstant(
            "OPALVERSION", buildinfo::version_major * 10000 + buildinfo::version_minor * 100
                                   + buildinfo::version_patch));
    opal->create(new RealConstant("RANK", ippl::Comm->rank()));
}

RealConstant::RealConstant(const std::string& name, RealConstant* parent)
    : ValueDefinition(name, parent) {}

RealConstant::RealConstant(const std::string& name, double value) : RealConstant() {
    setOpalName(name);
    Attributes::setReal(itsAttr[0], value);
    itsAttr[0].setReadOnly(true);
    builtin = true;
}

RealConstant::~RealConstant() {}

bool RealConstant::canReplaceBy(Object*) { return false; }

RealConstant* RealConstant::clone(const std::string& name) { return new RealConstant(name, this); }

double RealConstant::getReal() const { return Attributes::getReal(itsAttr[0]); }

void RealConstant::print(std::ostream& os) const {
    os << "REAL CONST " << getOpalName() << '=' << itsAttr[0] << ';';
    os << std::endl;
}

void RealConstant::printValue(std::ostream& os) const { os << itsAttr[0]; }
