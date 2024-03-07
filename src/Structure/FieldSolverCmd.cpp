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
        MX,    // mesh sixe in x
        MY,    // mesh sixe in y
        MZ,    // mesh sixe in z
        SIZE
    };
}

FieldSolverCmd::FieldSolverCmd()
    : Definition(
        SIZE, "FIELDSOLVER", "The \"FIELDSOLVER\" statement defines data for a the field solver") {
    itsAttr[TYPE] =
        Attributes::makePredefinedString("TYPE", "Name of the attached field solver.", {"NONE"});

    itsAttr[MX] = Attributes::makeReal("MX", "Meshsize in x");
    itsAttr[MY] = Attributes::makeReal("MY", "Meshsize in y");
    itsAttr[MZ] = Attributes::makeReal("MT", "Meshsize in z(t)");

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
    //    setFieldSolverCmdType();
    // update();
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

double FieldSolverCmd::getMX() const {
    return Attributes::getReal(itsAttr[MX]);
}

double FieldSolverCmd::getMY() const {
    return Attributes::getReal(itsAttr[MY]);
}

double FieldSolverCmd::getMZ() const {
    return Attributes::getReal(itsAttr[MZ]);
}

void FieldSolverCmd::setMX(double value) {
    Attributes::setReal(itsAttr[MX], value);
}

void FieldSolverCmd::setMY(double value) {
    Attributes::setReal(itsAttr[MY], value);
}

void FieldSolverCmd::setMZ(double value) {
    Attributes::setReal(itsAttr[MZ], value);
}

void FieldSolverCmd::update() {
}

void FieldSolverCmd::setFieldSolverCmdType() {
    static const std::map<std::string, FieldSolverCmdType> stringType_s = {
        {"NONE", FieldSolverCmdType::NONE}};

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
       << "* MX           " << Attributes::getReal(itsAttr[MX]) << '\n'
       << "* MY           " << Attributes::getReal(itsAttr[MY]) << '\n'
       << "* MZ           " << Attributes::getReal(itsAttr[MZ]) << '\n';
    os << "* ********************************************************************************** "
       << endl;
    return os;
}
