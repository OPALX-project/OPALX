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
            {"NONE", "FFT", "P3M", "OPEN", "CG", "FFT2D5", "BH"});

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

    // Attributes for the gridless Barnes-Hut solver
    itsAttr[FIELDSOLVER::BHTHETA] = Attributes::makeReal(
            "BHTHETA", "Barnes-Hut multipole acceptance angle in (0, 1] [BH only]", 0.5);
    itsAttr[FIELDSOLVER::BHSOFTENING] = Attributes::makeReal(
            "BHSOFTENING",
            "Barnes-Hut softening length h in m; pairs closer than 2h interact as if 2h apart. "
            "Required unless BHLEAFH=TRUE [BH only]",
            0.0);
    itsAttr[FIELDSOLVER::BHLEAFH] = Attributes::makeBool(
            "BHLEAFH",
            "TRUE to use each particle's octree leaf edge length as its softening [BH only]",
            false);
    itsAttr[FIELDSOLVER::BHDIRECT] = Attributes::makeBool(
            "BHDIRECT",
            "TRUE to evaluate the field by direct O(N^2) summation with the same softened kernel "
            "(and image charges) instead of the octree; for reference runs [BH only]",
            false);
    itsAttr[FIELDSOLVER::BHSOFTENINGREL] = Attributes::makeReal(
            "BHSOFTENINGREL",
            "Softening h as a multiple of the mean particle spacing in the solve (rest) frame, "
            "recomputed every solve from the rms sizes; replaces BHSOFTENING when > 0 [BH only]",
            0.0);
    itsAttr[FIELDSOLVER::BHBOXMARGIN] = Attributes::makeReal(
            "BHBOXMARGIN",
            "Margin per side of a newly built Barnes-Hut domain box, as a fraction of the bunch "
            "extent on that axis; the bunch can grow by 1 + 2 BHBOXMARGIN before the domain is "
            "rebuilt [BH only]",
            0.05);
    itsAttr[FIELDSOLVER::BHBOXSLACK] = Attributes::makeReal(
            "BHBOXSLACK",
            "The Barnes-Hut domain is rebuilt with a smaller box once its box exceeds a refit by "
            "more than this fraction on some axis [BH only]",
            0.25);
    itsAttr[FIELDSOLVER::BHSTATS] = Attributes::makeBool(
            "BHSTATS",
            "TRUE to print the Barnes-Hut interaction counts (P2P, M2P) of every solve and to time "
            "the traversal stages (bh.groups, bh.upsweep, bh.compute); costs a few collectives per "
            "solve [BH only]",
            false);
    itsAttr[FIELDSOLVER::BHPRECISION] = Attributes::makePredefinedString(
            "BHPRECISION",
            "Barnes-Hut precision: DOUBLE everywhere, or MIXED (double coordinates; float charges, "
            "softening, multipoles and accelerations) [BH only]",
            {"DOUBLE", "MIXED"}, "DOUBLE");
    itsAttr[FIELDSOLVER::BHRELATIVISTIC] = Attributes::makePredefinedString(
            "BHRELATIVISTIC",
            "Barnes-Hut source velocities: FRAME boosts one rest-frame solve with the mean "
            "momentum (like one energy bin); AUTO follows BINS: none is electrostatic, MAXBINS=1 "
            "is FRAME [BH only]",
            {"AUTO", "FRAME"}, "AUTO");

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

double FieldSolverCmd::getP3MCutoffCells() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::P3MRCUTCELLS]);
}

bool FieldSolverCmd::hasP3MCutoff() const { return !itsAttr[FIELDSOLVER::P3MRCUT].defaultUsed(); }

bool FieldSolverCmd::hasP3MCutoffCells() const {
    return !itsAttr[FIELDSOLVER::P3MRCUTCELLS].defaultUsed();
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

double FieldSolverCmd::getBHTheta() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::BHTHETA]);
}
double FieldSolverCmd::getBHSoftening() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::BHSOFTENING]);
}
bool FieldSolverCmd::getBHLeafBasedSoftening() const {
    return Attributes::getBool(itsAttr[FIELDSOLVER::BHLEAFH]);
}
bool FieldSolverCmd::getBHDirectSum() const {
    return Attributes::getBool(itsAttr[FIELDSOLVER::BHDIRECT]);
}
double FieldSolverCmd::getBHRelativeSoftening() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::BHSOFTENINGREL]);
}
double FieldSolverCmd::getBHBoxMargin() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::BHBOXMARGIN]);
}
double FieldSolverCmd::getBHBoxSlack() const {
    return Attributes::getReal(itsAttr[FIELDSOLVER::BHBOXSLACK]);
}
bool FieldSolverCmd::getBHStatistics() const {
    return Attributes::getBool(itsAttr[FIELDSOLVER::BHSTATS]);
}
std::string FieldSolverCmd::getBHPrecision() const {
    return Attributes::getString(itsAttr[FIELDSOLVER::BHPRECISION]);
}
std::string FieldSolverCmd::getBHRelativistic() const {
    return Attributes::getString(itsAttr[FIELDSOLVER::BHRELATIVISTIC]);
}
void FieldSolverCmd::setBHTheta(const double theta) {
    Attributes::setReal(itsAttr[FIELDSOLVER::BHTHETA], theta);
}
void FieldSolverCmd::setBHSoftening(const double softening) {
    Attributes::setReal(itsAttr[FIELDSOLVER::BHSOFTENING], softening);
}
void FieldSolverCmd::setBHLeafBasedSoftening(const bool enabled) {
    Attributes::setBool(itsAttr[FIELDSOLVER::BHLEAFH], enabled);
}
void FieldSolverCmd::setBHDirectSum(const bool enabled) {
    Attributes::setBool(itsAttr[FIELDSOLVER::BHDIRECT], enabled);
}
void FieldSolverCmd::setBHRelativeSoftening(const double factor) {
    Attributes::setReal(itsAttr[FIELDSOLVER::BHSOFTENINGREL], factor);
}
void FieldSolverCmd::setBHBoxMargin(const double margin) {
    Attributes::setReal(itsAttr[FIELDSOLVER::BHBOXMARGIN], margin);
}
void FieldSolverCmd::setBHBoxSlack(const double slack) {
    Attributes::setReal(itsAttr[FIELDSOLVER::BHBOXSLACK], slack);
}
void FieldSolverCmd::setBHStatistics(const bool enabled) {
    Attributes::setBool(itsAttr[FIELDSOLVER::BHSTATS], enabled);
}
void FieldSolverCmd::setBHPrecision(const std::string& precision) {
    Attributes::setPredefinedString(itsAttr[FIELDSOLVER::BHPRECISION], precision);
}
void FieldSolverCmd::setBHRelativistic(const std::string& mode) {
    Attributes::setPredefinedString(itsAttr[FIELDSOLVER::BHRELATIVISTIC], mode);
}

FieldSolverCmdType FieldSolverCmd::getFieldSolverCmdType() const {
    static const std::map<std::string, FieldSolverCmdType> types = {
            {"NONE", FieldSolverCmdType::NONE}, {"FFT", FieldSolverCmdType::FFT},
            {"P3M", FieldSolverCmdType::P3M},   {"OPEN", FieldSolverCmdType::OPEN},
            {"CG", FieldSolverCmdType::CG},     {"FFT2D5", FieldSolverCmdType::FFT2D5},
            {"BH", FieldSolverCmdType::BH}};
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

    if (getType() == "BH") {
        os << "* BHTHETA      " << getBHTheta() << '\n'
           << "* BHSOFTENING  " << getBHSoftening() << " [m]" << '\n'
           << "* BHLEAFH      " << (getBHLeafBasedSoftening() ? "TRUE" : "FALSE") << '\n'
           << "* BHDIRECT     " << (getBHDirectSum() ? "TRUE" : "FALSE") << '\n'
           << "* BHSOFTENINGREL " << getBHRelativeSoftening() << " [spacings]" << '\n'
           << "* BHBOXMARGIN  " << getBHBoxMargin() << '\n'
           << "* BHBOXSLACK   " << getBHBoxSlack() << '\n'
           << "* BHSTATS      " << (getBHStatistics() ? "TRUE" : "FALSE") << '\n'
           << "* BHPRECISION  " << getBHPrecision() << '\n'
           << "* BHRELATIVISTIC " << getBHRelativistic() << endl;
    }

    if (getType() == "P3M") {
        if (hasP3MCutoff()) {
            const double cutoff = getP3MCutoff();
            os << "* RCUT         " << cutoff << " [m]" << '\n'
               << "* ALPHA        " << 2.0 / cutoff << " [1/m]" << endl;
        } else {
            os << "* RCUT_CELLS   " << getP3MCutoffCells() << " [largest mesh spacing]" << '\n'
               << "* ALPHA        2/RCUT, updated with the mesh" << endl;
        }
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
