//
// Class FieldSolverCmd
//   The class for the OPAL FIELDSOLVER command.
//   Stores parsed mesh, boundary, and algorithm attributes. SpaceCharge owns
//   their conversion and solver compatibility validation.
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
#include "Structure/BinningCmd.h"
#include "Utilities/OpalException.h"

using namespace Expressions;

// TODO: o add a FIELD for DISCRETIZATION, MAXITERS, TOL...

FieldSolverCmd::FieldSolverCmd()
    : Definition(
              FIELDSOLVER::SIZE, "FIELDSOLVER",
              "The \"FIELDSOLVER\" statement defines data for a the field solver") {
    itsAttr[FIELDSOLVER::TYPE] = Attributes::makePredefinedString(
            "TYPE", "Name of the attached field solver.",
            {"NONE", "FFT", "P3M", "OPEN", "CG", "FFT2D5"});

    itsAttr[FIELDSOLVER::BINS] = Attributes::makeString(
            "BINS", "Name of BINNING definition to be used, or NONE for no binning.", "NONE");

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
            "BCFFTX", "Poisson-domain boundary in x; Dirichlet planes are configured separately.",
            {"OPEN", "DIRICHLET", "PERIODIC"}, "OPEN");

    itsAttr[FIELDSOLVER::BCFFTY] = Attributes::makePredefinedString(
            "BCFFTY", "Poisson-domain boundary in y.", {"OPEN", "DIRICHLET", "PERIODIC"}, "OPEN");

    itsAttr[FIELDSOLVER::BCFFTZ] = Attributes::makePredefinedString(
            "BCFFTZ", "Poisson-domain boundary in z; distinct from a Dirichlet plane.",
            {"OPEN", "DIRICHLET", "PERIODIC"}, "OPEN");

    itsAttr[FIELDSOLVER::GREENSF] = Attributes::makePredefinedString(
            "GREENSF", "Green function for TYPE=OPEN; TYPE=P3M selects its kernel internally.",
            {"STANDARD", "INTEGRATED"}, "INTEGRATED");

    itsAttr[FIELDSOLVER::P3MRCUT] = Attributes::makeReal(
            "RCUT",
            "Fixed P3M particle-particle cutoff radius in m (ALPHA = 2/RCUT); excludes "
            "RCUT_CELLS.",
            0.0);

    itsAttr[FIELDSOLVER::P3MRCUTCELLS] = Attributes::makeReal(
            "RCUT_CELLS",
            "P3M cutoff in units of the largest mesh spacing, updated with the mesh (ALPHA = "
            "2/cutoff); used unless RCUT is given.",
            2.0);

    itsAttr[FIELDSOLVER::PERIODICBOX] = Attributes::makeReal(
            "PERIODICBOX",
            "Edge length in m of a fixed all-periodic cube centered on the solve-frame origin; "
            "zero derives the box from the particle extents at every solve.",
            0.0);

    itsAttr[FIELDSOLVER::BBOXINCR] =
            Attributes::makeReal("BBOXINCR", "Increase of bounding box in % ", 2.0);

    // Attributes for FFT2D5 mode
    itsAttr[FIELDSOLVER::FFT2D5MODEL] = Attributes::makePredefinedString(
            "FFT2D5MODEL",
            "FFT2D5 model identity. ORIGINAL retains independent transverse 2D slices "
            "and the PIPEMODE longitudinal closure; numerical behavior is unchanged.",
            {"ORIGINAL"}, "ORIGINAL");
    itsAttr[FIELDSOLVER::PIPEMODE] = Attributes::makePredefinedString(
            "PIPEMODE",
            "Longitudinal-field geometry model; transverse Poisson slices remain open [FFT2D5].",
            {"OPEN", "CIRCULAR", "PLATES", "NONE"}, "OPEN");
    itsAttr[FIELDSOLVER::BEAMR] =
            Attributes::makeReal("BEAMR", "Beam radius in metres [FFT2D5 only]", 1.0);
    itsAttr[FIELDSOLVER::CLOSEDRING] =
            Attributes::makeBool("CLOSEDRING", "TRUE if the ring is closed [FFT2D5 only]", false);
    itsAttr[FIELDSOLVER::SCATTERLONGITUDINALLY] = Attributes::makeBool(
            "SCATTERLONGITUDINALLY",
            "TRUE to scatter charge longitudinally across slices [FFT2D5 only]", true);
    itsAttr[FIELDSOLVER::PIPESIZEX] = Attributes::makeReal(
            "PIPESIZEX", "Beam pipe horizontal size in metres [FFT2D5 only]", 1.0);
    itsAttr[FIELDSOLVER::PIPESIZEY] = Attributes::makeReal(
            "PIPESIZEY", "Beam pipe vertical size in metres [FFT2D5 only]", 1.0);
    itsAttr[FIELDSOLVER::REFPATHFNAME] =
            Attributes::makeString("REFPATHFNAME", "Reference path file name [FFT2D5 only]", "");

    // \todo does not work   registerOwnership(AttributeHandler::STATEMENT);
}

FieldSolverCmd::~FieldSolverCmd() {}

void FieldSolverCmd::execute() {}

FieldSolverCmd* FieldSolverCmd::find(const std::string& name) {
    FieldSolverCmd* fs = dynamic_cast<FieldSolverCmd*>(OpalData::getInstance()->find(name));

    if (fs == 0) {
        throw OpalException("FieldSolverCmd::find()", "FieldSolverCmd \"" + name + "\" not found.");
    }
    return fs;
}

std::string FieldSolverCmd::getType(const ObjectValues& values) { return values.getString("TYPE"); }

std::string FieldSolverCmd::getBinsName(const ObjectValues& values) {
    return values.getString("BINS");
}

std::string FieldSolverCmd::getGreensFunction(const ObjectValues& values) {
    return values.getString("GREENSF");
}

double FieldSolverCmd::getP3MCutoff(const ObjectValues& values) { return values.getReal("RCUT"); }

double FieldSolverCmd::getP3MCutoffCells(const ObjectValues& values) {
    return values.getReal("RCUT_CELLS");
}

double FieldSolverCmd::getPeriodicBox(const ObjectValues& values) {
    return values.getReal("PERIODICBOX");
}

bool FieldSolverCmd::hasP3MCutoff(const ObjectValues& values) {
    return !values.defaultUsed("RCUT");
}

bool FieldSolverCmd::hasP3MCutoffCells(const ObjectValues& values) {
    return !values.defaultUsed("RCUT_CELLS");
}

BCHandler<3> FieldSolverCmd::constructBCHandler(const ObjectValues& values) {
    using BCH_t = BCHandler<3>;

    BCH_t boundary_conditions(
            BCH_t::strToBCType(values.getString("BCFFTX")),
            BCH_t::strToBCType(values.getString("BCFFTY")),
            BCH_t::strToBCType(values.getString("BCFFTZ")));

    return boundary_conditions;
}

double FieldSolverCmd::getNX(const ObjectValues& values) { return values.getReal("NX"); }

double FieldSolverCmd::getNY(const ObjectValues& values) { return values.getReal("NY"); }

double FieldSolverCmd::getNZ(const ObjectValues& values) { return values.getReal("NZ"); }

void FieldSolverCmd::setNX(double value) { Attributes::setReal(itsAttr[FIELDSOLVER::NX], value); }

void FieldSolverCmd::setNY(double value) { Attributes::setReal(itsAttr[FIELDSOLVER::NY], value); }

void FieldSolverCmd::setNZ(double value) { Attributes::setReal(itsAttr[FIELDSOLVER::NZ], value); }

double FieldSolverCmd::getBoxIncr(const ObjectValues& values) { return values.getReal("BBOXINCR"); }

std::string FieldSolverCmd::getFFT2D5Model(const ObjectValues& values) {
    return values.getString("FFT2D5MODEL");
}
bool FieldSolverCmd::hasFFT2D5Model(const ObjectValues& values) {
    return !values.defaultUsed("FFT2D5MODEL");
}
void FieldSolverCmd::setFFT2D5Model(const std::string& model) {
    Attributes::setPredefinedString(itsAttr[FIELDSOLVER::FFT2D5MODEL], model);
}
std::string FieldSolverCmd::getPipeMode(const ObjectValues& values) {
    return values.getString("PIPEMODE");
}
double FieldSolverCmd::getBeamRadius(const ObjectValues& values) { return values.getReal("BEAMR"); }
bool FieldSolverCmd::getClosedRing(const ObjectValues& values) {
    return values.getBool("CLOSEDRING");
}
bool FieldSolverCmd::getScatterLongitudinally(const ObjectValues& values) {
    return values.getBool("SCATTERLONGITUDINALLY");
}
double FieldSolverCmd::getPipeSizeX(const ObjectValues& values) {
    return values.getReal("PIPESIZEX");
}
double FieldSolverCmd::getPipeSizeY(const ObjectValues& values) {
    return values.getReal("PIPESIZEY");
}
std::string FieldSolverCmd::getRefPathFileName(const ObjectValues& values) {
    return values.getString("REFPATHFNAME");
}
void FieldSolverCmd::setPipeMode(const std::string& pipeMode) {
    Attributes::setPredefinedString(itsAttr[FIELDSOLVER::PIPEMODE], pipeMode);
}
void FieldSolverCmd::setBeamRadius(const double beamRadius) {
    Attributes::setReal(itsAttr[FIELDSOLVER::BEAMR], beamRadius);
}
void FieldSolverCmd::setClosedRing(const bool closedRing) {
    Attributes::setBool(itsAttr[FIELDSOLVER::CLOSEDRING], closedRing);
}
void FieldSolverCmd::setScatterLongitudinally(const bool val) {
    Attributes::setBool(itsAttr[FIELDSOLVER::SCATTERLONGITUDINALLY], val);
}
void FieldSolverCmd::setPipeSizeX(const double pipeSizeX) {
    Attributes::setReal(itsAttr[FIELDSOLVER::PIPESIZEX], pipeSizeX);
}
void FieldSolverCmd::setPipeSizeY(const double pipeSizeY) {
    Attributes::setReal(itsAttr[FIELDSOLVER::PIPESIZEY], pipeSizeY);
}
void FieldSolverCmd::setRefPathFileName(const std::string& refPathFileName) {
    Attributes::setString(itsAttr[FIELDSOLVER::REFPATHFNAME], refPathFileName);
}

FieldSolverCmdType FieldSolverCmd::getFieldSolverCmdType(const ObjectValues& values) {
    static const std::map<std::string, FieldSolverCmdType> types = {
            {"NONE", FieldSolverCmdType::NONE}, {"FFT", FieldSolverCmdType::FFT},
            {"P3M", FieldSolverCmdType::P3M},   {"OPEN", FieldSolverCmdType::OPEN},
            {"CG", FieldSolverCmdType::CG},     {"FFT2D5", FieldSolverCmdType::FFT2D5}};
    const auto found = types.find(getType(values));
    if (found == types.end()) {
        throw OpalException("FieldSolverCmd::getFieldSolverCmdType", "Unknown or missing TYPE.");
    }
    return found->second;
}

ippl::Vector<bool, 3> FieldSolverCmd::getDomainDecomposition(const ObjectValues& values) {
    return {values.getBool("PARFFTX"), values.getBool("PARFFTY"), values.getBool("PARFFTZ")};
}

BinningCmd* FieldSolverCmd::getBinningCmd() const {
    const std::string binsName = Attributes::getString(itsAttr[FIELDSOLVER::BINS]);
    if (binsName == "NONE" || binsName.empty()) {
        return nullptr;
    }
    return BinningCmd::find(binsName);
}

Inform& FieldSolverCmd::printInfo(Inform& os, const ObjectValues& values) {
    os << "* ************* F I E L D S O L V E R ********************************************** "
       << endl;
    os << "* FIELDSOLVER  " << values.getName() << '\n'
       << "* BINS         " << getBinsName(values) << '\n'
       << "* TYPE         " << getType(values) << '\n'
       << "* RANKS        " << ippl::Comm->size() << '\n'
       << "* NX           " << values.getReal("NX") << '\n'
       << "* NY           " << values.getReal("NY") << '\n'
       << "* NZ           " << values.getReal("NZ") << '\n'
       << "* BBOXINCR     " << values.getReal("BBOXINCR") << '\n'
       << "* GREENSF      " << values.getString("GREENSF") << endl;

    if (getType(values) == "FFT2D5") {
        os << "* FFT2D5MODEL  " << getFFT2D5Model(values) << endl;
    }

    if (getType(values) == "P3M") {
        if (hasP3MCutoff(values)) {
            const double cutoff = getP3MCutoff(values);
            os << "* RCUT         " << cutoff << " [m]" << '\n'
               << "* ALPHA        " << 2.0 / cutoff << " [1/m]" << endl;
        } else {
            os << "* RCUT_CELLS   " << getP3MCutoffCells(values) << " [largest mesh spacing]"
               << '\n'
               << "* ALPHA        2/RCUT, updated with the mesh" << endl;
        }
    }

    if (values.getBool("PARFFTX")) {
        os << "* XDIM         parallel  " << endl;
    } else {
        os << "* XDIM         serial  " << endl;
    }

    if (values.getBool("PARFFTY")) {
        os << "* YDIM         parallel  " << endl;
    } else {
        os << "* YDIM         serial  " << endl;
    }

    if (values.getBool("PARFFTZ")) {
        os << "* ZDIM         parallel  " << endl;
    } else {
        os << "* ZDIM         serial  " << endl;
    }

    os << "* ********************************************************************************** "
       << endl;
    return os;
}
