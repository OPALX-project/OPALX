//
// Function writeEMFieldDumps
//   Writes the external fields of the tracked line on the grids of the DUMPEMFIELDS requests.
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
#include "Structure/EMFieldDumpWriter.h"

#include <fstream>

#include "Elements/OpalBeamline.h"
#include "Fields/Interpolation/NDGrid.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"

extern Inform* gmsg;

namespace {

    void writeHeader(std::ofstream& fout, int numPoints) {
        fout << numPoints << "\n";
        fout << 1 << "  x [m]\n";
        fout << 2 << "  y [m]\n";
        fout << 3 << "  z [m]\n";
        fout << 4 << "  t [ns]\n";
        fout << 5 << "  Bx [T]\n";
        fout << 6 << "  By [T]\n";
        fout << 7 << "  Bz [T]\n";
        fout << 8 << "  Ex [MV/m]\n";
        fout << 9 << "  Ey [MV/m]\n";
        fout << 10 << " Ez [MV/m]\n";
        fout << 0 << std::endl;
    }

    void writeFieldLine(
            const ElementSet& elements, const Vector_t<double, 3>& point, const double& time,
            std::ofstream& fout) {
        Vector_t<double, 3> E{};
        Vector_t<double, 3> B{};
        // The grid time is in ns; the elements take seconds.
        OpalBeamline::applyOnHost(elements, point, time * Units::ns2s, E, B);
        fout << point[0] << " " << point[1] << " " << point[2] << " " << time << " ";
        fout << B[0] << " " << B[1] << " " << B[2] << " ";
        fout << E[0] * Units::Vpm2MVpm << " " << E[1] * Units::Vpm2MVpm << " "
             << E[2] * Units::Vpm2MVpm << "\n";
    }

    void print(Inform& os, const opalx::run::EMFieldDumpDescription& dump) {
        os << "* ************* D U M P  E M  F I E L D S "
              "****************************************** "
           << endl;
        os << "* File name: '" << dump.fileName << "'\n";
        const char* axes[] = {"X", "Y", "Z", "T"};
        for (int i = 0; i < 4; ++i) {
            const char* unit = i < 3 ? " [m]\n" : " [ns]\n";
            os << "* " << axes[i] << "_START   = " << dump.start[i] << unit << "* D" << axes[i]
               << "        = " << dump.spacing[i] << unit << "* " << axes[i]
               << "_STEPS   = " << dump.steps[i] << '\n';
        }
        os << "* "
              "********************************************************************************** "
           << endl;
    }

}  // namespace

void writeEMFieldDumps(
        const std::vector<opalx::run::EMFieldDumpDescription>& dumps, const ElementSet& elements) {
    for (const opalx::run::EMFieldDumpDescription& dump : dumps) {
        print(*gmsg << level5, dump);

        std::array<int, 4> steps      = dump.steps;
        std::array<double, 4> spacing = dump.spacing;
        std::array<double, 4> start   = dump.start;
        interpolation::NDGrid grid(4, steps.data(), spacing.data(), start.data());

        std::ofstream fout(dump.fileName, std::ofstream::out);
        if (!fout.good()) {
            throw OpalException(
                    "writeEMFieldDumps", "Failed to open DumpEMFields file " + dump.fileName);
        }
        writeHeader(fout, grid.end().toInteger());
        std::vector<double> position(4);
        Vector_t<double, 3> point(0., 0., 0.);
        for (interpolation::Mesh::Iterator it = grid.begin(); it < grid.end(); ++it) {
            it.getPosition(&position[0]);
            for (size_t i = 0; i < 3; ++i) {
                point[i] = position[i];
            }
            writeFieldLine(elements, point, position[3], fout);
        }
        if (!fout.good()) {
            throw OpalException(
                    "writeEMFieldDumps", "Something went wrong during writing " + dump.fileName);
        }
    }
}
