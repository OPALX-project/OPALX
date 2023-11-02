//
// Class FieldSolver
//   The class for the OPAL FIELDSOLVER command.
//   A FieldSolver definition is used by most physics commands to define the
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
#include "Structure/FieldSolver.h"

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

// The attributes of class FieldSolver.
namespace {
    enum {
        TYPE,  // The field solver name
        MX,    // mesh sixe in x
        MY,    // mesh sixe in y
        MZ,    // mesh sixe in z
        SIZE
    };
}

FieldSolver::FieldSolver()
    : Definition(
        SIZE, "FIELDSOLVER", "The \"FIELDSOLVER\" statement defines data for a the field solver") {
    itsAttr[TYPE] =
        Attributes::makePredefinedString("TYPE", "Name of the attached field solver.", {"NONE"});

    itsAttr[MX] = Attributes::makeReal("MX", "Meshsize in x");
    itsAttr[MY] = Attributes::makeReal("MY", "Meshsize in y");
    itsAttr[MZ] = Attributes::makeReal("MT", "Meshsize in z(t)");

    // \todo does not work   registerOwnership(AttributeHandler::STATEMENT);
}

FieldSolver::FieldSolver(const std::string& name, FieldSolver* parent) : Definition(name, parent) {
}

FieldSolver::~FieldSolver() {
}

FieldSolver* FieldSolver::clone(const std::string& name) {
    return new FieldSolver(name, this);
}

void FieldSolver::execute() {
    //    setFieldSolverType();
    // update();
}

FieldSolver* FieldSolver::find(const std::string& name) {
    FieldSolver* fs = dynamic_cast<FieldSolver*>(OpalData::getInstance()->find(name));

    if (fs == 0) {
        throw OpalException("FieldSolver::find()", "FieldSolver \"" + name + "\" not found.");
    }
    return fs;
}

std::string FieldSolver::getType() {
    return Attributes::getString(itsAttr[TYPE]);
}

double FieldSolver::getMX() const {
    return Attributes::getReal(itsAttr[MX]);
}

double FieldSolver::getMY() const {
    return Attributes::getReal(itsAttr[MY]);
}

double FieldSolver::getMZ() const {
    return Attributes::getReal(itsAttr[MZ]);
}

void FieldSolver::setMX(double value) {
    Attributes::setReal(itsAttr[MX], value);
}

void FieldSolver::setMY(double value) {
    Attributes::setReal(itsAttr[MY], value);
}

void FieldSolver::setMZ(double value) {
    Attributes::setReal(itsAttr[MZ], value);
}

void FieldSolver::update() {
}

void FieldSolver::setFieldSolverType() {
    static const std::map<std::string, FieldSolverType> stringType_s = {
        {"NONE", FieldSolverType::NONE}};

    fsName_m = getType();
    if (fsName_m.empty()) {
        throw OpalException(
            "FieldSolver::setFieldSolverType",
            "The attribute \"TYPE\" isn't set for \"FIELDSOLVER\"!");
    } else {
        fsType_m = stringType_s.at(fsName_m);
    }
}

bool FieldSolver::hasValidSolver() {
    return false;
}

Inform& FieldSolver::printInfo(Inform& os) const {
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
