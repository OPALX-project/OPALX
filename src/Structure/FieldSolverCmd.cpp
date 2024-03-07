//
// Class FieldSolverCmd
//   The class for the OPAL FIELDSOLVER command.
//   A FieldSolverCmd definition is used by most physics commands to define the
//   particle charge and the reference momentum, together with some other data.
//
// Copyright (c) 200x - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
//
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
#include "Structure/FieldSolverCmd.h"

#include <map>
#include "AbstractObjects/Element.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

using namespace Expressions;

// TODO: o add a FIELD for DISCRETIZATION, MAXITERS, TOL...

// The attributes of class FieldSolverCmd.
namespace {
    enum {
        TYPE,  // The field solver name
        NX,    // mesh sixe in x
        NY,    // mesh sixe in y
        NZ,    // mesh sixe in z
        SIZE
    };
}

FieldSolverCmd::FieldSolverCmd()
    : Definition(
        SIZE, "FIELDSOLVER", "The \"FIELDSOLVER\" statement defines data for a the field solver") {
    itsAttr[TYPE] =
        Attributes::makePredefinedString("TYPE", "Name of the attached field solver.", {"NONE"});

    itsAttr[NX] = Attributes::makeReal("NX", "Meshsize in x");
    itsAttr[NY] = Attributes::makeReal("NY", "Meshsize in y");
    itsAttr[NZ] = Attributes::makeReal("NZ", "Meshsize in z");

    // \todo does not work   registerOwnership(AttributeHandler::STATEMENT);
}

FieldSolverCmd::FieldSolverCmd(const std::string& name, FieldSolverCmd* parent)
    : Definition(name, parent) {
}

FieldSolverCmd::~FieldSolverCmd() {
}

FieldSolverCmd* FieldSolverCmd::clone(const std::string& name) {
    return new FieldSolverCmd(name, this);
}

void FieldSolverCmd::execute() {
    setFieldSolverCmdType();
    update();
}

FieldSolverCmd* FieldSolverCmd::find(const std::string& name) {
    FieldSolverCmd* fs = dynamic_cast<FieldSolverCmd*>(OpalData::getInstance()->find(name));

    if (fs == 0) {
        throw OpalException("FieldSolverCmd::find()", "FieldSolverCmd \"" + name + "\" not found.");
    }
    return fs;
}

std::string FieldSolverCmd::getType() {
    return Attributes::getString(itsAttr[TYPE]);
}

double FieldSolverCmd::getNX() const {
    return Attributes::getReal(itsAttr[NX]);
}

double FieldSolverCmd::getNY() const {
    return Attributes::getReal(itsAttr[NY]);
}

double FieldSolverCmd::getNZ() const {
    return Attributes::getReal(itsAttr[NZ]);
}

void FieldSolverCmd::setNX(double value) {
    Attributes::setReal(itsAttr[NX], value);
}

void FieldSolverCmd::setNY(double value) {
    Attributes::setReal(itsAttr[NY], value);
}

void FieldSolverCmd::setNZ(double value) {
    Attributes::setReal(itsAttr[NZ], value);
}

void FieldSolverCmd::update() {
    if (itsAttr[TYPE]) {
        fsName_m = getType();
    }
}

void FieldSolverCmd::setFieldSolverCmdType() {
    static const std::map<std::string, FieldSolverCmdType> stringType_s = {
        {"NONE", FieldSolverCmdType::NONE},
        {"FFT", FieldSolverCmdType::FFT},
    };

    fsName_m = getType();

    if (fsName_m.empty()) {
        throw OpalException(
            "FieldSolverCmd::setFieldSolverCmdType",
            "The attribute \"TYPE\" isn't set for \"FIELDSOLVER\"!");
    } else {
        fsType_m = stringType_s.at(fsName_m);
    }
}

bool FieldSolverCmd::hasValidSolver() {
    return false;
}

Inform& FieldSolverCmd::printInfo(Inform& os) const {
    os << "* ************* F I E L D S O L V E R ********************************************** "
       << endl;
    os << "* FIELDSOLVER  " << getOpalName() << '\n'
       << "* TYPE         " << fsName_m << '\n'
       << "* N-PROCESSORS " << ippl::Comm->size() << '\n'
       << "* NX           " << Attributes::getReal(itsAttr[NX]) << '\n'
       << "* NY           " << Attributes::getReal(itsAttr[NY]) << '\n'
       << "* NZ           " << Attributes::getReal(itsAttr[NZ]) << '\n';
    os << "* ********************************************************************************** "
       << endl;
    return os;
}
