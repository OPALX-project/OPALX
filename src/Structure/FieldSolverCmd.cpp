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
            "BCFFTX",
            "Poisson-domain boundary in x; source-plane corrections are configured separately.",
            {"OPEN", "DIRICHLET", "PERIODIC"}, "OPEN");

    itsAttr[FIELDSOLVER::BCFFTY] = Attributes::makePredefinedString(
            "BCFFTY", "Poisson-domain boundary in y.", {"OPEN", "DIRICHLET", "PERIODIC"}, "OPEN");

    itsAttr[FIELDSOLVER::BCFFTZ] = Attributes::makePredefinedString(
            "BCFFTZ", "Poisson-domain boundary in z; distinct from a source-plane correction.",
            {"OPEN", "DIRICHLET", "PERIODIC"}, "OPEN");

    itsAttr[FIELDSOLVER::GREENSF] = Attributes::makePredefinedString(
            "GREENSF", "Green function for TYPE=OPEN; TYPE=P3M selects its kernel internally.",
            {"STANDARD", "INTEGRATED"}, "INTEGRATED");

    itsAttr[FIELDSOLVER::P3MRCUT] = Attributes::makeReal(
            "RCUT", "P3M particle-particle cutoff radius in m (ALPHA is derived as 2/RCUT).", 0.0);

    itsAttr[FIELDSOLVER::BBOXINCR] =
            Attributes::makeReal("BBOXINCR", "Increase of bounding box in % ", 2.0);

    // Attributes for FFT2D5 mode
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

FieldSolverCmd::FieldSolverCmd(const std::string& name, FieldSolverCmd* parent)
    : Definition(name, parent) {}

FieldSolverCmd::~FieldSolverCmd() {}

FieldSolverCmd* FieldSolverCmd::clone(const std::string& name) {
    return new FieldSolverCmd(name, this);
}

void FieldSolverCmd::execute() { update(); }

FieldSolverCmd* FieldSolverCmd::find(const std::string& name) {
    FieldSolverCmd* fs = dynamic_cast<FieldSolverCmd*>(OpalData::getInstance()->find(name));

    if (fs == 0) {
        throw OpalException("FieldSolverCmd::find()", "FieldSolverCmd \"" + name + "\" not found.");
    }
    return fs;
}

std::string FieldSolverCmd::getType() const {
    return Attributes::getString(itsAttr[FIELDSOLVER::TYPE]);
}

std::string FieldSolverCmd::getBinsName() const {
    return Attributes::getString(itsAttr[FIELDSOLVER::BINS]);
}

std::string FieldSolverCmd::getGreensFunction() const {
    return Attributes::getString(itsAttr[FIELDSOLVER::GREENSF]);
}

double FieldSolverCmd::getP3MCutoff() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::P3MRCUT]);
}

BCHandler<3> FieldSolverCmd::constructBCHandler() const {
    using BCH_t = BCHandler<3>;

    BCH_t boundary_conditions(
            BCH_t::strToBCType(Attributes::getString(itsAttr[FIELDSOLVER::BCFFTX])),
            BCH_t::strToBCType(Attributes::getString(itsAttr[FIELDSOLVER::BCFFTY])),
            BCH_t::strToBCType(Attributes::getString(itsAttr[FIELDSOLVER::BCFFTZ])));

    return boundary_conditions;
}

double FieldSolverCmd::getNX() const { return Attributes::getReal(itsAttr[FIELDSOLVER::NX]); }

double FieldSolverCmd::getNY() const { return Attributes::getReal(itsAttr[FIELDSOLVER::NY]); }

double FieldSolverCmd::getNZ() const { return Attributes::getReal(itsAttr[FIELDSOLVER::NZ]); }

void FieldSolverCmd::setNX(double value) { Attributes::setReal(itsAttr[FIELDSOLVER::NX], value); }

void FieldSolverCmd::setNY(double value) { Attributes::setReal(itsAttr[FIELDSOLVER::NY], value); }

void FieldSolverCmd::setNZ(double value) { Attributes::setReal(itsAttr[FIELDSOLVER::NZ], value); }

double FieldSolverCmd::getBoxIncr() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::BBOXINCR]);
}

void FieldSolverCmd::update() {}

std::string FieldSolverCmd::getPipeMode() const {
    return Attributes::getString(itsAttr[FIELDSOLVER::PIPEMODE]);
}
double FieldSolverCmd::getBeamRadius() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::BEAMR]);
}
bool FieldSolverCmd::getClosedRing() const {
    return Attributes::getBool(itsAttr[FIELDSOLVER::CLOSEDRING]);
}
bool FieldSolverCmd::getScatterLongitudinally() const {
    return Attributes::getBool(itsAttr[FIELDSOLVER::SCATTERLONGITUDINALLY]);
}
double FieldSolverCmd::getPipeSizeX() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::PIPESIZEX]);
}
double FieldSolverCmd::getPipeSizeY() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::PIPESIZEY]);
}
std::string FieldSolverCmd::getRefPathFileName() const {
    return Attributes::getString(itsAttr[FIELDSOLVER::REFPATHFNAME]);
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

FieldSolverCmdType FieldSolverCmd::getFieldSolverCmdType() const {
    static const std::map<std::string, FieldSolverCmdType> types = {
            {"NONE", FieldSolverCmdType::NONE}, {"FFT", FieldSolverCmdType::FFT},
            {"P3M", FieldSolverCmdType::P3M},   {"OPEN", FieldSolverCmdType::OPEN},
            {"CG", FieldSolverCmdType::CG},     {"FFT2D5", FieldSolverCmdType::FFT2D5}};
    const auto found = types.find(getType());
    if (found == types.end()) {
        throw OpalException("FieldSolverCmd::getFieldSolverCmdType", "Unknown or missing TYPE.");
    }
    return found->second;
}

ippl::Vector<bool, 3> FieldSolverCmd::getDomainDecomposition() const {
    return {Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTX]),
            Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTY]),
            Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTZ])};
}

BinningCmd* FieldSolverCmd::getBinningCmd() const {
    const std::string binsName = getBinsName();
    if (binsName == "NONE" || binsName.empty()) {
        return nullptr;
    }
    return BinningCmd::find(binsName);
}

Inform& FieldSolverCmd::printInfo(Inform& os) const {
    os << "* ************* F I E L D S O L V E R ********************************************** "
       << endl;
    os << "* FIELDSOLVER  " << getOpalName() << '\n'
       << "* BINS         " << getBinsName() << '\n'
       << "* TYPE         " << getType() << '\n'
       << "* RANKS        " << ippl::Comm->size() << '\n'
       << "* NX           " << Attributes::getReal(itsAttr[FIELDSOLVER::NX]) << '\n'
       << "* NY           " << Attributes::getReal(itsAttr[FIELDSOLVER::NY]) << '\n'
       << "* NZ           " << Attributes::getReal(itsAttr[FIELDSOLVER::NZ]) << '\n'
       << "* BBOXINCR     " << Attributes::getReal(itsAttr[FIELDSOLVER::BBOXINCR]) << '\n'
       << "* GREENSF      " << Attributes::getString(itsAttr[FIELDSOLVER::GREENSF]) << endl;

    if (getType() == "P3M") {
        const double cutoff = getP3MCutoff();
        os << "* RCUT         " << cutoff << " [m]" << '\n'
           << "* ALPHA        " << 2.0 / cutoff << " [1/m]" << endl;
    }

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
        os << "* ZDIM         parallel  " << endl;
    } else {
        os << "* ZDIM         serial  " << endl;
    }

    os << "* ********************************************************************************** "
       << endl;
    return os;
}
