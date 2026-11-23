//
// Class DumpEMFields
//   The DUMPEMFIELDS command: a request to write the fields of the lattice on a
//   user-defined Cartesian grid.
//
// Copyright (c) 2017, Chris Rogers
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
#include "BasicActions/DumpEMFields.h"
#include <cmath>
#include <filesystem>
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

std::vector<opalx::run::EMFieldDumpDescription> DumpEMFields::requests_m;

DumpEMFields::DumpEMFields()
    : Action(SIZE, "DUMPEMFIELDS",
             "The \"DUMPEMFIELDS\" statement dumps a field map to a user-defined "
             "field file, for checking that fields are generated correctly. "
             "The fields are written out on a grid in space and time.") {
    // would be nice if "steps" could be integer
    itsAttr[FILE_NAME] =
            Attributes::makeString("FILE_NAME", "Name of the file to which field data is dumped");

    itsAttr[COORDINATE_SYSTEM] = Attributes::makePredefinedString(
            "COORDINATE_SYSTEM", "Coordinate system of the grid; only CARTESIAN is supported",
            {"CARTESIAN"}, "CARTESIAN");

    itsAttr[X_START] = Attributes::makeReal("X_START", "Start point in the grid in x [m]");

    itsAttr[DX] = Attributes::makeReal("DX", "Grid step size in x [m]");

    itsAttr[X_STEPS] = Attributes::makeReal("X_STEPS", "Number of steps in x");

    itsAttr[Y_START] =
            Attributes::makeReal("Y_START", "Start point in the grid in y (vertical) [m]");

    itsAttr[DY] = Attributes::makeReal("DY", "Grid step size in y [m]");

    itsAttr[Y_STEPS] = Attributes::makeReal("Y_STEPS", "Number of steps in y");

    itsAttr[Z_START] = Attributes::makeReal("Z_START", "Start point in the grid in z [m]");

    itsAttr[DZ] = Attributes::makeReal("DZ", "Grid step size in z [m]");

    itsAttr[Z_STEPS] = Attributes::makeReal("Z_STEPS", "Number of steps in z");

    itsAttr[T_START] = Attributes::makeReal("T_START", "Start point in the grid in time [ns]");

    itsAttr[DT] = Attributes::makeReal("DT", "Grid step size in time [ns]");

    itsAttr[T_STEPS] = Attributes::makeReal("T_STEPS", "Number of steps in time");

    registerOwnership(AttributeHandler::STATEMENT);
}

void DumpEMFields::execute() { requests_m.push_back(describe()); }

const std::vector<opalx::run::EMFieldDumpDescription>& DumpEMFields::requests() {
    return requests_m;
}

opalx::run::EMFieldDumpDescription DumpEMFields::describe() const {
    opalx::run::EMFieldDumpDescription dump;
    const int start[]   = {X_START, Y_START, Z_START, T_START};
    const int spacing[] = {DX, DY, DZ, DT};
    const int steps[]   = {X_STEPS, Y_STEPS, Z_STEPS, T_STEPS};
    for (int i = 0; i < 4; ++i) {
        dump.start[i]       = Attributes::getReal(itsAttr[start[i]]);
        dump.spacing[i]     = Attributes::getReal(itsAttr[spacing[i]]);
        const double nSteps = Attributes::getReal(itsAttr[steps[i]]);
        checkInt(nSteps, itsAttr[steps[i]].getName());
        dump.steps[i] = static_cast<int>(nSteps);
    }

    const std::string fileName = Attributes::getString(itsAttr[FILE_NAME]);
    dump.fileName =
            std::filesystem::path(fileName).is_absolute()
                    ? fileName
                    : Util::combineFilePath(
                              {OpalData::getInstance()->getAuxiliaryOutputDirectory(), fileName});
    return dump;
}

void DumpEMFields::checkInt(double value, const std::string& name, const double tolerance) {
    value += tolerance;  // prevent rounding error
    if (std::abs(std::floor(value) - value) > 2 * tolerance) {
        throw OpalException(
                "DumpEMFields::checkInt",
                "Value for " + name + " should be an integer but a real value was found");
    }
    if (std::floor(value) < 0.5) {
        throw OpalException("DumpEMFields::checkInt", "Value for " + name + " should be 1 or more");
    }
}

void DumpEMFields::print(std::ostream& os) const {
    os << "* ************* D U M P  E M  F I E L D S ****************************************** "
       << std::endl;
    os << "* File name: '" << Attributes::getString(itsAttr[FILE_NAME]) << "'\n";
    os << "* Coordinate system: " << Attributes::getString(itsAttr[COORDINATE_SYSTEM]) << '\n'
       << "* X_START   = " << Attributes::getReal(itsAttr[X_START]) << " [m]\n"
       << "* DX        = " << Attributes::getReal(itsAttr[DX]) << " [m]\n"
       << "* X_STEPS   = " << Attributes::getReal(itsAttr[X_STEPS]) << '\n'
       << "* Y_START   = " << Attributes::getReal(itsAttr[Y_START]) << " [m]\n"
       << "* DY        = " << Attributes::getReal(itsAttr[DY]) << " [m]\n"
       << "* Y_STEPS   = " << Attributes::getReal(itsAttr[Y_STEPS]) << '\n'
       << "* Z_START   = " << Attributes::getReal(itsAttr[Z_START]) << " [m]\n"
       << "* DZ        = " << Attributes::getReal(itsAttr[DZ]) << " [m]\n"
       << "* Z_STEPS   = " << Attributes::getReal(itsAttr[Z_STEPS]) << '\n';
    os << "* T_START   = " << Attributes::getReal(itsAttr[T_START]) << " [ns]\n"
       << "* DT        = " << Attributes::getReal(itsAttr[DT]) << " [ns]\n"
       << "* T_STEPS   = " << Attributes::getReal(itsAttr[T_STEPS]) << '\n';
    os << "* ********************************************************************************** "
       << std::endl;
}
