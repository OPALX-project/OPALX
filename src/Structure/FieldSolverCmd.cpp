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

FieldSolverCmd::FieldSolverCmd()
    : Definition(
        FIELDSOLVER::SIZE, "FIELDSOLVER",
        "The \"FIELDSOLVER\" statement defines data for a the field solver") {
    itsAttr[FIELDSOLVER::TYPE] = Attributes::makePredefinedString(
        "TYPE", "Name of the attached field solver.", {"NONE", "FFT", "FFTPERIODIC", "CG ", "OPEN", "P3M"});

    itsAttr[FIELDSOLVER::NX] = Attributes::makeReal("NX", "Meshsize in x");
    itsAttr[FIELDSOLVER::NY] = Attributes::makeReal("NY", "Meshsize in y");
    itsAttr[FIELDSOLVER::NZ] = Attributes::makeReal("NZ", "Meshsize in z");

    itsAttr[FIELDSOLVER::PARFFTX] =
        Attributes::makeBool("PARFFTX", "True, dimension 0 i.e x is parallelized", false);

    itsAttr[FIELDSOLVER::PARFFTY] =
        Attributes::makeBool("PARFFTY", "True, dimension 1 i.e y is parallelized", false);

    itsAttr[FIELDSOLVER::PARFFTZ] =
        Attributes::makeBool("PARFFTZ", "True, dimension 2 i.e z is parallelized", true);

    itsAttr[FIELDSOLVER::BCFFTX] = Attributes::makePredefinedString(
        "BCFFTX", "Boundary conditions in x.", {"OPEN", "DIRICHLET", "PERIODIC"}, "OPEN");

    itsAttr[FIELDSOLVER::BCFFTY] = Attributes::makePredefinedString(
        "BCFFTY", "Boundary conditions in y.", {"OPEN", "DIRICHLET", "PERIODIC"}, "OPEN");

    itsAttr[FIELDSOLVER::BCFFTZ] = Attributes::makePredefinedString(
        "BCFFTZ", "Boundary conditions in z.", {"OPEN", "DIRICHLET", "PERIODIC"}, "OPEN");

    itsAttr[FIELDSOLVER::GREENSF] = Attributes::makePredefinedString(
        "GREENSF", "Which Greensfunction to be used.", {"STANDARD", "INTEGRATED"}, "INTEGRATED");

    itsAttr[FIELDSOLVER::BBOXINCR] =
        Attributes::makeReal("BBOXINCR", "Increase of bounding box in % ", 2.0);

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
    return Attributes::getString(itsAttr[FIELDSOLVER::TYPE]);
}

double FieldSolverCmd::getNX() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::NX]);
}

double FieldSolverCmd::getNY() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::NY]);
}

double FieldSolverCmd::getNZ() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::NZ]);
}

void FieldSolverCmd::setNX(double value) {
    Attributes::setReal(itsAttr[FIELDSOLVER::NX], value);
}

void FieldSolverCmd::setNY(double value) {
    Attributes::setReal(itsAttr[FIELDSOLVER::NY], value);
}

void FieldSolverCmd::setNZ(double value) {
    Attributes::setReal(itsAttr[FIELDSOLVER::NZ], value);
}

double FieldSolverCmd::getBoxIncr() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::BBOXINCR]);
}


void FieldSolverCmd::update() {
    if (itsAttr[FIELDSOLVER::TYPE]) {
        fsName_m = getType();
    }
}

void FieldSolverCmd::setFieldSolverCmdType() {
    static const std::map<std::string, FieldSolverCmdType> stringType_s = {
        {"NONE", FieldSolverCmdType::NONE},
        {"FFT", FieldSolverCmdType::FFT},
        {"FFTPERIODIC", FieldSolverCmdType::FFTPERIODIC},
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
       << "* RANKS       " << ippl::Comm->size() << '\n'
       << "* NX           " << Attributes::getReal(itsAttr[FIELDSOLVER::NX]) << '\n'
       << "* NY           " << Attributes::getReal(itsAttr[FIELDSOLVER::NY]) << '\n'
       << "* NZ           " << Attributes::getReal(itsAttr[FIELDSOLVER::NZ]) << '\n'
       << "* BBOXINCR     " << Attributes::getReal(itsAttr[FIELDSOLVER::BBOXINCR]) << '\n'
       << "* GREENSF      " << Attributes::getString(itsAttr[FIELDSOLVER::GREENSF]) << endl;

    if (Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTX])) {
        os << "* XDIM         parallel  " << endl;
    } else {
        os << "* XDIM         serial  " << endl;
    }

    if (Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTY])) {
        os << "* YDIM         parallel  " << endl;
    } else {
        os << "* YDIM         serial  " << endl;
    }

    if (Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTZ])) {
        os << "* ZDIM      parallel  " << endl;
    } else {
        os << "* ZDIM      serial  " << endl;
    }

    os << "* ********************************************************************************** "
       << endl;
    return os;
}
