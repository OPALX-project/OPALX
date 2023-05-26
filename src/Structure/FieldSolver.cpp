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

//TODO: o add a FIELD for DISCRETIZATION, MAXITERS, TOL...

// The attributes of class FieldSolver.
namespace {
    namespace deprecated {
        enum {
            BCFFTT,
            SIZE
        };
    }
    enum {
        FSTYPE = deprecated::SIZE,   // The field solver name
        // FOR FFT BASED SOLVER
        MX,         // mesh sixe in x
        MY,         // mesh sixe in y
        MT,         // mesh sixe in z
        PARFFTX,    // parallelized grid in x
        PARFFTY,    // parallelized grid in y
        PARFFTT,    // parallelized grid in z
        BCFFTX,     // boundary condition in x [FFT]
        BCFFTY,     // boundary condition in y [FFT]
        BCFFTZ,     // boundary condition in z [FFT]
        GREENSF,    // holds greensfunction to be used [FFT only]
        BBOXINCR,   // how much the boundingbox is increased
        RC,         // cutoff radius for PP interactions
        ALPHA,      // Green’s function splitting parameter
        EPSILON,    // regularization for PP interaction
        SIZE
    };
}


FieldSolver::FieldSolver():
    Definition(SIZE, "FIELDSOLVER",
               "The \"FIELDSOLVER\" statement defines data for a the field solver") {

    itsAttr[FSTYPE] = Attributes::makePredefinedString("FSTYPE", "Name of the attached field solver.",
                                                       {"FFT", "FFTPERIODIC", "NONE"});

    itsAttr[MX] = Attributes::makeReal("MX", "Meshsize in x");
    itsAttr[MY] = Attributes::makeReal("MY", "Meshsize in y");
    itsAttr[MT] = Attributes::makeReal("MT", "Meshsize in z(t)");

    itsAttr[PARFFTX] = Attributes::makeBool("PARFFTX",
                                            "True, dimension 0 i.e x is parallelized",
                                            false);

    itsAttr[PARFFTY] = Attributes::makeBool("PARFFTY",
                                            "True, dimension 1 i.e y is parallelized",
                                            false);

    itsAttr[PARFFTT] = Attributes::makeBool("PARFFTT",
                                            "True, dimension 2 i.e z(t) is parallelized",
                                            true);

    //FFT ONLY:
    itsAttr[BCFFTX] = Attributes::makePredefinedString("BCFFTX",
                                                       "Boundary conditions in x.",
                                                       {"OPEN", "DIRICHLET", "PERIODIC"},
                                                       "OPEN");

    itsAttr[BCFFTY] = Attributes::makePredefinedString("BCFFTY",
                                                       "Boundary conditions in y.",
                                                       {"OPEN", "DIRICHLET", "PERIODIC"},
                                                       "OPEN");

    itsAttr[BCFFTZ] = Attributes::makePredefinedString("BCFFTZ",
                                                       "Boundary conditions in z(t).",
                                                       {"OPEN", "DIRICHLET", "PERIODIC"},
                                                       "OPEN");

    itsAttr[deprecated::BCFFTT] = Attributes::makePredefinedString("BCFFTT",
                                                                  "Boundary conditions in z(t).",
                                                                  {"OPEN", "DIRICHLET", "PERIODIC"},
                                                                  "OPEN");

    itsAttr[GREENSF]  = Attributes::makePredefinedString("GREENSF",
                                                         "Which Greensfunction to be used.",
                                                         {"STANDARD", "INTEGRATED"},
                                                         "INTEGRATED");

    itsAttr[BBOXINCR] = Attributes::makeReal("BBOXINCR",
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
    return Attributes::getString(itsAttr[FSTYPE]);
}

double FieldSolver::getMX() const {
    return Attributes::getReal(itsAttr[MX]);
}

double FieldSolver::getMY() const {
    return Attributes::getReal(itsAttr[MY]);
}

double FieldSolver::getMT() const {
    return Attributes::getReal(itsAttr[MT]);
}

void FieldSolver::setMX(double value) {
    Attributes::setReal(itsAttr[MX], value);
}

void FieldSolver::setMY(double value) {
    Attributes::setReal(itsAttr[MY], value);
}

void FieldSolver::setMT(double value) {
    Attributes::setReal(itsAttr[MT], value);
}

void FieldSolver::update() {

}

void FieldSolver::initCartesianFields() {

    e_dim_tag decomp[3] = {SERIAL, SERIAL, SERIAL};

    NDIndex<3> domain;
    domain[0] = Index((int)getMX() + 1);
    domain[1] = Index((int)getMY() + 1);
    domain[2] = Index((int)getMT() + 1);

    if (Attributes::getBool(itsAttr[PARFFTX]))
        decomp[0] = PARALLEL;
    if (Attributes::getBool(itsAttr[PARFFTY]))
        decomp[1] = PARALLEL;
    if (Attributes::getBool(itsAttr[PARFFTT]))
        decomp[2] = PARALLEL;

    if (Attributes::getString(itsAttr[FSTYPE]) == "FFTPERIODIC") {
        decomp[0] = decomp[1] = SERIAL;
        decomp[2] = PARALLEL;
    }
    // create prototype mesh and layout objects for this problem domain

    mesh_m   = new Mesh_t(domain);
    FL_m     = new FieldLayout_t(*mesh_m, decomp);
    PL_m.reset(new Layout_t(*FL_m, *mesh_m));
}

bool FieldSolver::hasPeriodicZ() {
    if (itsAttr[deprecated::BCFFTT])
        return (Attributes::getString(itsAttr[deprecated::BCFFTT]) == "PERIODIC");

    return (Attributes::getString(itsAttr[BCFFTZ]) == "PERIODIC");
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
                            "The attribute \"FSTYPE\" isn't set for \"FIELDSOLVER\"!");
    } else {
        fsType_m = stringFSType_s.at(fsName_m);
    }
}

void FieldSolver::initSolver(PartBunchBase<double, 3>* b) {
    itsBunch_m = b;

    std::string greens = Attributes::getString(itsAttr[GREENSF]);
    std::string bcx = Attributes::getString(itsAttr[BCFFTX]);
    std::string bcy = Attributes::getString(itsAttr[BCFFTY]);
    std::string bcz = Attributes::getString(itsAttr[deprecated::BCFFTT]);
    if (bcz.empty()) {
        bcz = Attributes::getString(itsAttr[BCFFTZ]);
    }

    /* Boundary geometry must be set */

}

bool FieldSolver::hasValidSolver() {
    return (solver_m != 0);
}

Inform& FieldSolver::printInfo(Inform& os) const {
    os << "* ************* F I E L D S O L V E R ********************************************** " << endl;
    os << "* FIELDSOLVER  " << getOpalName() << '\n'
       << "* TYPE         " << fsName_m << '\n'
       << "* N-PROCESSORS " << Ippl::getNodes() << '\n'
       << "* MX           " << Attributes::getReal(itsAttr[MX])   << '\n'
       << "* MY           " << Attributes::getReal(itsAttr[MY])   << '\n'
       << "* MT           " << Attributes::getReal(itsAttr[MT])   << '\n'
       << "* BBOXINCR     " << Attributes::getReal(itsAttr[BBOXINCR]) << '\n'
       << "* GRRENSF      " << Attributes::getString(itsAttr[GREENSF]) << endl;
    

    if (Attributes::getBool(itsAttr[PARFFTX])) {
        os << "* XDIM         parallel  " << endl;
    } else {
        os << "* XDIM         serial  " << endl;
    }

    if (Attributes::getBool(itsAttr[PARFFTY])) {
        os << "* YDIM         parallel  " << endl;
    } else {
        os << "* YDIM         serial  " << endl;
    }

    if (Attributes::getBool(itsAttr[PARFFTT])) {
        os << "* Z(T)DIM      parallel  " << endl;
    }else {
        os << "* Z(T)DIM      serial  " << endl;
    }

    if (solver_m)
        os << *solver_m << endl;
    os << "* ********************************************************************************** " << endl;
    return os;
}

