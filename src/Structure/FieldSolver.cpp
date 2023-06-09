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

#include "AbstractObjects/Element.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Solvers/FFTBoxPoissonSolver.h"
#include "Solvers/FFTPoissonSolver.h"
#include "Structure/BoundaryGeometry.h"
#include "Utilities/OpalException.h"

#include <map>

using namespace Expressions;

// The attributes of class FieldSolver.
namespace FIELDSOLVER {
    enum {
        TYPE,   // The field solver name
        MX,         // mesh sixe in x
        MY,         // mesh sixe in y
        MZ,         // mesh sixe in z
        PARFFTX,    // parallelized grid in x
        PARFFTY,    // parallelized grid in y
        PARFFTZ,    // parallelized grid in z
        BCFFTX,     // boundary condition in x [FFT]
        BCFFTY,     // boundary condition in y [FFT]
        BCFFTZ,     // boundary condition in z [FFT]
        GREENSF,    // holds greensfunction to be used [FFT only]
        BBOXINCR,   // how much the boundingbox is increased
        SIZE
    };
}


FieldSolver::FieldSolver():
    Definition(FIELDSOLVER::SIZE, "FIELDSOLVER",
               "The \"FIELDSOLVER\" statement defines data for a the field solver") {

    itsAttr[FIELDSOLVER::TYPE] = Attributes::makePredefinedString("TYPE", "Name of the attached field solver.",
                                                       {"FFT", "FFTPERIODIC", "NONE"});

    itsAttr[FIELDSOLVER::MX] = Attributes::makeReal("MX", "Meshsize in x");
    itsAttr[FIELDSOLVER::MY] = Attributes::makeReal("MY", "Meshsize in y");
    itsAttr[FIELDSOLVER::MZ] = Attributes::makeReal("MZ", "Meshsize in z");

    itsAttr[FIELDSOLVER::PARFFTX] = Attributes::makeBool("PARFFTX",
                                            "True, dimension 0 i.e x is parallelized",
                                            false);

    itsAttr[FIELDSOLVER::PARFFTY] = Attributes::makeBool("PARFFTY",
                                            "True, dimension 1 i.e y is parallelized",
                                            false);

    itsAttr[FIELDSOLVER::PARFFTZ] = Attributes::makeBool("PARFFTZ",
                                            "True, dimension 2 i.e z is parallelized",
                                            true);

    //FFT ONLY:
    itsAttr[FIELDSOLVER::BCFFTX] = Attributes::makePredefinedString("BCFFTX",
                                                       "Boundary conditions in x.",
                                                       {"OPEN", "DIRICHLET", "PERIODIC"},
                                                       "OPEN");

    itsAttr[FIELDSOLVER::BCFFTY] = Attributes::makePredefinedString("BCFFTY",
                                                       "Boundary conditions in y.",
                                                       {"OPEN", "DIRICHLET", "PERIODIC"},
                                                       "OPEN");

    itsAttr[FIELDSOLVER::BCFFTZ] = Attributes::makePredefinedString("BCFFTZ",
                                                       "Boundary conditions in z.",
                                                       {"OPEN", "DIRICHLET", "PERIODIC"},
                                                       "OPEN");

    itsAttr[FIELDSOLVER::GREENSF]  = Attributes::makePredefinedString("GREENSF",
                                                         "Which Greensfunction to be used.",
                                                         {"STANDARD", "INTEGRATED"},
                                                         "INTEGRATED");

    itsAttr[FIELDSOLVER::BBOXINCR] = Attributes::makeReal("BBOXINCR",
                                             "Increase of bounding box in % ",
                                             2.0);

    mesh_m = 0;
    FL_m = 0;
    PL_m.reset(nullptr);
    solver_m = 0;

    registerOwnership(AttributeHandler::STATEMENT);
}


FieldSolver::FieldSolver(const std::string& name, FieldSolver* parent):
    Definition(name, parent)
{
    mesh_m = 0;
    FL_m = 0;
    PL_m.reset(nullptr);
    solver_m = 0;
}


FieldSolver::~FieldSolver() {
    if (mesh_m) {
        delete mesh_m;
        mesh_m = 0;
    }
    if (FL_m) {
        delete FL_m;
        FL_m = 0;
    }
    if (solver_m) {
       delete solver_m;
       solver_m = 0;
    }
}

FieldSolver* FieldSolver::clone(const std::string& name) {
    return new FieldSolver(name, this);
}

void FieldSolver::execute() {
    setFieldSolverType();
    update();
}

FieldSolver* FieldSolver::find(const std::string& name) {
    FieldSolver* fs = dynamic_cast<FieldSolver*>(OpalData::getInstance()->find(name));

    if (fs == 0) {
        throw OpalException("FieldSolver::find()", "FieldSolver \"" + name + "\" not found.");
    }
    return fs;
}

std::string FieldSolver::getType() {
    return Attributes::getString(itsAttr[FIELDSOLVER::TYPE]);
}

double FieldSolver::getMX() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::MX]);
}

double FieldSolver::getMY() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::MY]);
}

double FieldSolver::getMZ() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::MZ]);
}

void FieldSolver::setMX(double value) {
    Attributes::setReal(itsAttr[FIELDSOLVER::MX], value);
}

void FieldSolver::setMY(double value) {
    Attributes::setReal(itsAttr[FIELDSOLVER::MY], value);
}

void FieldSolver::setMZ(double value) {
    Attributes::setReal(itsAttr[FIELDSOLVER::MZ], value);
}

void FieldSolver::update() {

}

void FieldSolver::initCartesianFields() {

    e_dim_tag decomp[3] = {SERIAL, SERIAL, SERIAL};

    NDIndex<3> domain;
    domain[0] = Index((int)getMX() + 1);
    domain[1] = Index((int)getMY() + 1);
    domain[2] = Index((int)getMZ() + 1);

    if (Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTX]))
        decomp[0] = PARALLEL;
    if (Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTY]))
        decomp[1] = PARALLEL;
    if (Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTZ]))
        decomp[2] = PARALLEL;

    if (Attributes::getString(itsAttr[FIELDSOLVER::TYPE]) == "FFTPERIODIC") {
        decomp[0] = decomp[1] = SERIAL;
        decomp[2] = PARALLEL;
    }
    // create prototype mesh and layout objects for this problem domain

    mesh_m   = new Mesh_t(domain);
    FL_m     = new FieldLayout_t(*mesh_m, decomp);
    PL_m.reset(new Layout_t(*FL_m, *mesh_m));
}

bool FieldSolver::hasPeriodicZ() {

    return (Attributes::getString(itsAttr[FIELDSOLVER::BCFFTZ]) == "PERIODIC");
}

void FieldSolver::setFieldSolverType() {
    static const std::map<std::string, FieldSolverType> stringFSType_s = {
        {"NONE",   FieldSolverType::NONE},
        {"FFT",    FieldSolverType::FFT},
        {"FFTBOX", FieldSolverType::FFTBOX}
    };

    fsName_m = getType();
    if (fsName_m.empty()) {
        throw OpalException("FieldSolver::setFieldSolverType",
                            "The attribute \"TYPE\" isn't set for \"FIELDSOLVER\"!");
    } else {
        fsType_m = stringFSType_s.at(fsName_m);
    }
}

void FieldSolver::initSolver(PartBunchBase<double, 3>* b) {
    itsBunch_m = b;

    std::string greens = Attributes::getString(itsAttr[FIELDSOLVER::GREENSF]);
    std::string bcx = Attributes::getString(itsAttr[FIELDSOLVER::BCFFTX]);
    std::string bcy = Attributes::getString(itsAttr[FIELDSOLVER::BCFFTY]);
    std::string bcz = Attributes::getString(itsAttr[FIELDSOLVER::BCFFTZ]);
}

bool FieldSolver::hasValidSolver() {
    return (solver_m != 0);
}

Inform& FieldSolver::printInfo(Inform& os) const {
    os << "* ************* F I E L D S O L V E R ********************************************** " << endl;
    os << "* FIELDSOLVER  " << getOpalName() << '\n'
       << "* TYPE         " << fsName_m << '\n'
       << "* N-PROCESSORS " << Ippl::getNodes() << '\n'
       << "* MX           " << Attributes::getReal(itsAttr[FIELDSOLVER::MX])   << '\n'
       << "* MY           " << Attributes::getReal(itsAttr[FIELDSOLVER::MY])   << '\n'
       << "* MZ           " << Attributes::getReal(itsAttr[FIELDSOLVER::MZ])   << '\n'
       << "* BBOXINCR     " << Attributes::getReal(itsAttr[FIELDSOLVER::BBOXINCR]) << '\n'
       << "* GRRENSF      " << Attributes::getString(itsAttr[FIELDSOLVER::GREENSF]) << endl;
    

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
    }else {
        os << "* ZDIM      serial  " << endl;
    }

    if (solver_m)
        os << *solver_m << endl;
    os << "* ********************************************************************************** " << endl;
    return os;
}

