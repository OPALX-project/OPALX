//
// Class G4BL3DMagnetoStatic
//   Reader for G4beamline `grid` field maps.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#include "Fields/G4BL3DMagnetoStatic.h"
#include "Fields/G4BLMapSyntax.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Units.h"
#include "Utilities/GeneralOpalException.h"

#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using G4BLMapSyntax::firstToken;
using G4BLMapSyntax::keyValues;
using G4BLMapSyntax::nextLine;

G4BL3DMagnetoStatic::G4BL3DMagnetoStatic(std::string aFilename)
    : Fieldmap(aFilename),
      xbegin_m(0.0),
      xend_m(0.0),
      ybegin_m(0.0),
      yend_m(0.0),
      zbegin_m(0.0),
      zend_m(-1e-3),
      hx_m(0.0),
      hy_m(0.0),
      hz_m(0.0),
      num_gridpx_m(0),
      num_gridpy_m(0),
      num_gridpz_m(0),
      fieldScale_m(1.0) {
    Type        = TG4BL3DMagnetoStatic;
    normalize_m = false;  // absolute Tesla, so the element's scale is a plain multiplier

    readHeaderOnly();
}

G4BL3DMagnetoStatic::~G4BL3DMagnetoStatic() { freeMap(); }

void G4BL3DMagnetoStatic::readHeaderOnly() {
    std::ifstream in(Filename_m.c_str());
    if (!in.good()) {
        noFieldmapWarning();
        throw GeneralOpalException(
                "G4BL3DMagnetoStatic::readHeaderOnly",
                "Could not open fieldmap '" + Filename_m + "'");
    }

    std::string line;
    if (!nextLine(in, line)) {
        throw GeneralOpalException(
                "G4BL3DMagnetoStatic::readHeaderOnly", "Fieldmap '" + Filename_m + "' is empty");
    }

    // ---- optional `param` line ------------------------------------------
    fieldScale_m        = 1.0;
    std::string keyword = firstToken(line);
    if (keyword == "param") {
        const std::map<std::string, std::string> pairs = keyValues(line);
        auto optional = [&](const std::string& key, double fallback) -> double {
            const auto entry = pairs.find(key);
            if (entry == pairs.end()) {
                return fallback;
            }
            try {
                return std::stod(entry->second);
            } catch (const std::exception&) {
                throw GeneralOpalException(
                        "G4BL3DMagnetoStatic::readHeaderOnly", "Could not read " + key + "='"
                                                                       + entry->second + "' in '"
                                                                       + Filename_m + "'");
            }
        };

        // G4beamline evaluates B * normB * current_deck / current_file. Fold the file's own
        // half in here, so the element's scale is exactly the deck-side `current=` and a
        // scale of 1 matches a bare `fieldmap ...` placement.
        const double normB   = optional("normB", 1.0);
        const double current = optional("current", 1.0);
        if (current == 0.0) {
            throw GeneralOpalException(
                    "G4BL3DMagnetoStatic::readHeaderOnly",
                    "The 'param' line of '" + Filename_m
                            + "' sets current=0, which G4beamline would divide by");
        }
        fieldScale_m = normB / current;

        if (!nextLine(in, line)) {
            throw GeneralOpalException(
                    "G4BL3DMagnetoStatic::readHeaderOnly",
                    "Fieldmap '" + Filename_m + "' has no section after its 'param' line");
        }
        keyword = firstToken(line);
    }

    // ---- `grid` line -----------------------------------------------------
    if (keyword != "grid") {
        throw GeneralOpalException(
                "G4BL3DMagnetoStatic::readHeaderOnly",
                "Expected a 'grid' section in '" + Filename_m + "', found '" + keyword + "'");
    }

    const std::map<std::string, std::string> grid = keyValues(line);
    auto require                                  = [&](const std::string& key) -> double {
        const auto entry = grid.find(key);
        if (entry == grid.end()) {
            throw GeneralOpalException(
                    "G4BL3DMagnetoStatic::readHeaderOnly",
                    "The 'grid' line of '" + Filename_m + "' has no " + key);
        }
        try {
            return std::stod(entry->second);
        } catch (const std::exception&) {
            throw GeneralOpalException(
                    "G4BL3DMagnetoStatic::readHeaderOnly",
                    "Could not read " + key + "='" + entry->second + "' in '" + Filename_m + "'");
        }
    };

    const double X0 = require("X0");
    const double Y0 = require("Y0");
    const double Z0 = require("Z0");
    const double dX = require("dX");
    const double dY = require("dY");
    const double dZ = require("dZ");
    num_gridpx_m    = static_cast<int>(require("nX"));
    num_gridpy_m    = static_cast<int>(require("nY"));
    num_gridpz_m    = static_cast<int>(require("nZ"));

    if (num_gridpx_m < 2 || num_gridpy_m < 2 || num_gridpz_m < 2) {
        throw GeneralOpalException(
                "G4BL3DMagnetoStatic::readHeaderOnly",
                "'" + Filename_m + "' needs at least 2 grid points in each of x, y and z");
    }
    if (dX <= 0.0 || dY <= 0.0 || dZ <= 0.0) {
        throw GeneralOpalException(
                "G4BL3DMagnetoStatic::readHeaderOnly",
                "'" + Filename_m + "' has a non-positive dX, dY or dZ");
    }

    hx_m     = dX * Units::mm2m;
    hy_m     = dY * Units::mm2m;
    hz_m     = dZ * Units::mm2m;
    xbegin_m = X0 * Units::mm2m;
    ybegin_m = Y0 * Units::mm2m;
    zbegin_m = Z0 * Units::mm2m;
    xend_m   = xbegin_m + (num_gridpx_m - 1) * hx_m;
    yend_m   = ybegin_m + (num_gridpy_m - 1) * hy_m;
    zend_m   = zbegin_m + (num_gridpz_m - 1) * hz_m;

    // The section after `grid` has to be `data`. Reject the alternatives by name so the
    // message says what is wrong rather than "unexpected line".
    if (!nextLine(in, line)) {
        throw GeneralOpalException(
                "G4BL3DMagnetoStatic::readHeaderOnly",
                "Fieldmap '" + Filename_m + "' has no section after its 'grid' line");
    }
    const std::string section = firstToken(line);
    if (section != "data") {
        std::string reason;
        if (section.compare(0, 6, "extend") == 0) {
            reason = "symmetry declarations ('" + section + "') are not supported";
        } else if (section == "points") {
            reason = "the 'points' section is not supported, only 'data'";
        } else {
            reason = "expected a 'data' section, found '" + section + "'";
        }
        throw GeneralOpalException(
                "G4BL3DMagnetoStatic::readHeaderOnly",
                "In fieldmap '" + Filename_m + "': " + reason);
    }
}

void G4BL3DMagnetoStatic::readMap() {
    if (FieldstrengthBz_m.extent(0) != 0) {
        return;
    }

    const size_t size = static_cast<size_t>(num_gridpx_m) * num_gridpy_m * num_gridpz_m;
    FieldstrengthBx_m = Kokkos::DualView<double*>("FieldstrengthBx", size);
    FieldstrengthBy_m = Kokkos::DualView<double*>("FieldstrengthBy", size);
    FieldstrengthBz_m = Kokkos::DualView<double*>("FieldstrengthBz", size);

    auto Bx = FieldstrengthBx_m.view_host();
    auto By = FieldstrengthBy_m.view_host();
    auto Bz = FieldstrengthBz_m.view_host();

    std::ifstream in(Filename_m.c_str());
    if (!in.good()) {
        throw GeneralOpalException(
                "G4BL3DMagnetoStatic::readMap", "Could not reopen fieldmap '" + Filename_m + "'");
    }

    // Skip forward to just past the `data` line. readHeaderOnly() has already checked that
    // the file is well formed this far, so this only has to find the marker.
    std::string line;
    while (nextLine(in, line)) {
        if (firstToken(line) == "data") {
            break;
        }
    }

    // Rows carry their own coordinates, so they are placed by those rather than by their
    // position in the file. `seen` catches a repeated or a missing row: without it, one bad
    // row would shift the rest of the map silently.
    std::vector<bool> seen(size, false);
    size_t filled = 0;

    // Snapping a coordinate to its grid index. A quarter of a cell is generous for files
    // written with six digits and tight enough to catch a genuinely off-grid row.
    auto index1D = [](double value, double begin, double h, int n, const char* axis,
                      const std::string& file, size_t row) -> int {
        const double f = (value - begin) / h;
        const int i    = static_cast<int>(std::lround(f));
        if (i < 0 || i >= n || std::abs(f - i) > 0.25) {
            throw GeneralOpalException(
                    "G4BL3DMagnetoStatic::readMap",
                    "Row " + std::to_string(row) + " of '" + file + "' has " + axis + " = "
                            + std::to_string(value)
                            + " m, which is not on the grid declared by the 'grid' line");
        }
        return i;
    };

    size_t row = 0;
    while (nextLine(in, line)) {
        const std::string label = firstToken(line);
        if (!label.empty() && (std::isalpha(static_cast<unsigned char>(label[0])) != 0)) {
            std::string reason;
            if (label.compare(0, 6, "extend") == 0) {
                reason = "symmetry declarations ('" + label + "') are not supported";
            } else {
                reason = "unexpected line '" + line + "' inside the data block";
            }
            throw GeneralOpalException(
                    "G4BL3DMagnetoStatic::readMap", "In fieldmap '" + Filename_m + "': " + reason);
        }

        ++row;
        std::istringstream values(line);
        double x = 0.0, y = 0.0, z = 0.0, bx = 0.0, by = 0.0, bz = 0.0;
        if (!(values >> x >> y >> z >> bx >> by >> bz)) {
            throw GeneralOpalException(
                    "G4BL3DMagnetoStatic::readMap",
                    "Row " + std::to_string(row) + " of '" + Filename_m
                            + "' does not hold at least six numbers (x y z Bx By Bz)");
        }

        // G4beamline also writes a nine-column form, x y z Bx By Bz Ex Ey Ez -- the muE4
        // quadrupole map qsm01a_210_track.g4blmap is one. This element is magnetostatic, so
        // an electric field cannot be carried: accept the columns only when they are zero,
        // and say so plainly rather than dropping a real field in silence.
        double ex = 0.0, ey = 0.0, ez = 0.0;
        if (values >> ex) {
            if (!(values >> ey >> ez)) {
                throw GeneralOpalException(
                        "G4BL3DMagnetoStatic::readMap",
                        "Row " + std::to_string(row) + " of '" + Filename_m
                                + "' has between seven and eight numbers; the G4beamline grid "
                                  "format is six (x y z Bx By Bz) or nine with Ex Ey Ez");
            }
            if (ex != 0.0 || ey != 0.0 || ez != 0.0) {
                throw GeneralOpalException(
                        "G4BL3DMagnetoStatic::readMap",
                        "Row " + std::to_string(row) + " of '" + Filename_m
                                + "' carries a non-zero electric field. This reader is "
                                  "magnetostatic and would silently discard it.");
            }
            std::string extra;
            if (values >> extra) {
                throw GeneralOpalException(
                        "G4BL3DMagnetoStatic::readMap", "Row " + std::to_string(row) + " of '"
                                                                + Filename_m
                                                                + "' has more than nine numbers");
            }
        }

        const int ix = index1D(x * Units::mm2m, xbegin_m, hx_m, num_gridpx_m, "x", Filename_m, row);
        const int iy = index1D(y * Units::mm2m, ybegin_m, hy_m, num_gridpy_m, "y", Filename_m, row);
        const int iz = index1D(z * Units::mm2m, zbegin_m, hz_m, num_gridpz_m, "z", Filename_m, row);

        const size_t index = (static_cast<size_t>(ix) * num_gridpy_m + iy) * num_gridpz_m + iz;
        if (seen[index]) {
            throw GeneralOpalException(
                    "G4BL3DMagnetoStatic::readMap",
                    "Row " + std::to_string(row) + " of '" + Filename_m
                            + "' repeats a grid point already given earlier in the file");
        }
        seen[index] = true;
        ++filled;

        Bx(index) = fieldScale_m * bx;
        By(index) = fieldScale_m * by;
        Bz(index) = fieldScale_m * bz;
    }

    if (filled != size) {
        throw GeneralOpalException(
                "G4BL3DMagnetoStatic::readMap",
                "Fieldmap '" + Filename_m + "' holds " + std::to_string(filled) + " of the "
                        + std::to_string(size) + " grid points its 'grid' line declares");
    }

    FieldstrengthBx_m.modify<typename decltype(FieldstrengthBx_m)::host_mirror_space>();
    FieldstrengthBx_m.sync<typename decltype(FieldstrengthBx_m)::t_dev::device_type>();
    FieldstrengthBy_m.modify<typename decltype(FieldstrengthBy_m)::host_mirror_space>();
    FieldstrengthBy_m.sync<typename decltype(FieldstrengthBy_m)::t_dev::device_type>();
    FieldstrengthBz_m.modify<typename decltype(FieldstrengthBz_m)::host_mirror_space>();
    FieldstrengthBz_m.sync<typename decltype(FieldstrengthBz_m)::t_dev::device_type>();
}

void G4BL3DMagnetoStatic::freeMap() {
    FieldstrengthBx_m = Kokkos::DualView<double*>();
    FieldstrengthBy_m = Kokkos::DualView<double*>();
    FieldstrengthBz_m = Kokkos::DualView<double*>();
}

void G4BL3DMagnetoStatic::applyField(std::shared_ptr<ParticleContainer_t> pc, double scale) {
    // Members copied to locals; the device kernel must not capture `this`.
    const double xbegin = xbegin_m, xend = xend_m;
    const double ybegin = ybegin_m, yend = yend_m;
    const double zbegin = zbegin_m, zend = zend_m;
    const double hx = hx_m, hy = hy_m, hz = hz_m;
    const int nx = num_gridpx_m, ny = num_gridpy_m, nz = num_gridpz_m;

    Kokkos::View<const double*> Bx_device = FieldstrengthBx_m.view_device();
    Kokkos::View<const double*> By_device = FieldstrengthBy_m.view_device();
    Kokkos::View<const double*> Bz_device = FieldstrengthBz_m.view_device();

    auto Rview          = pc->R.getView();
    auto Bview          = pc->B.getView();
    const size_t nLocal = pc->getLocalNum();

    Kokkos::parallel_for(
            "G4BL3DMagnetoStatic::applyField", nLocal, KOKKOS_LAMBDA(const size_t i) {
                const Vector_t<double, 3>& R = Rview(i);
                if (R(0) >= xbegin && R(0) < xend && R(1) >= ybegin && R(1) < yend && R(2) >= zbegin
                    && R(2) < zend) {
                    Vector_t<double, 3> tmpB = 0.0;
                    computeField(
                            R, tmpB, Bx_device, By_device, Bz_device, xbegin, ybegin, zbegin, hx,
                            hy, hz, nx, ny, nz);
                    Bview(i) += scale * tmpB;
                }
            });
}

bool G4BL3DMagnetoStatic::getFieldstrength(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B) const {
    if (!isInside(R)) {
        return true;
    }
    computeField(
            R, B, Kokkos::View<const double*>(FieldstrengthBx_m.view_host()),
            Kokkos::View<const double*>(FieldstrengthBy_m.view_host()),
            Kokkos::View<const double*>(FieldstrengthBz_m.view_host()), xbegin_m, ybegin_m,
            zbegin_m, hx_m, hy_m, hz_m, num_gridpx_m, num_gridpy_m, num_gridpz_m);
    return false;
}

bool G4BL3DMagnetoStatic::getFieldDerivative(
        const Vector_t<double, 3>& /*R*/, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/,
        const DiffDirection& /*dir*/) const {
    throw GeneralOpalException("G4BL3DMagnetoStatic::getFieldDerivative", "not implemented");
}

void G4BL3DMagnetoStatic::getFieldDimensions(double& zBegin, double& zEnd) const {
    zBegin = zbegin_m;
    zEnd   = zend_m;
}

void G4BL3DMagnetoStatic::getFieldDimensions(
        double& xIni, double& xFinal, double& yIni, double& yFinal, double& zIni,
        double& zFinal) const {
    xIni   = xbegin_m;
    xFinal = xend_m;
    yIni   = ybegin_m;
    yFinal = yend_m;
    zIni   = zbegin_m;
    zFinal = zend_m;
}

void G4BL3DMagnetoStatic::swap() {
    // The grid line carries its own axes, so there is nothing to swap.
}

void G4BL3DMagnetoStatic::getInfo(Inform* msg) {
    (*msg) << Filename_m << " (G4beamline grid, 3D magnetostatic); x= " << xbegin_m << " .. "
           << xend_m << " m; y= " << ybegin_m << " .. " << yend_m << " m; z= " << zbegin_m << " .. "
           << zend_m << " m; " << num_gridpx_m << " x " << num_gridpy_m << " x " << num_gridpz_m
           << " points;" << endl;
}

double G4BL3DMagnetoStatic::getFrequency() const {
    throw GeneralOpalException("G4BL3DMagnetoStatic::getFrequency", "not implemented");
    return 0.0;
}

void G4BL3DMagnetoStatic::setFrequency(double /*freq*/) {
    throw GeneralOpalException("G4BL3DMagnetoStatic::setFrequency", "not implemented");
}
