/* 
 *  Copyright (c) 2012, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are met: 
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer. 
 *  2. Redistributions in binary form must reproduce the above copyright notice, 
 *     this list of conditions and the following disclaimer in the documentation 
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to 
 *     endorse or promote products derived from this software without specific 
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include <math.h>

#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "Utilities/LogicalError.h"

#include "Fields/SectorMagneticFieldMap/Mesh.h"
#include "Fields/SectorMagneticFieldMap/ThreeDGrid.h"
#include "Fields/SectorMagneticFieldMap/Interpolator3dGridTo3d.h"
#include "Fields/SectorMagneticFieldMap.h"

// allow a fairly generous phi tolerance - we don't care about phi much and
// calculation can be flaky due to ascii truncation of double and conversions
// from polar to Cartesian
const double SectorMagneticFieldMap::fractionalBBPhiTolerance_m(1e-2);
std::map<std::string, SectorMagneticFieldMap*> SectorMagneticFieldMap::_fields;

SectorMagneticFieldMap::SectorMagneticFieldMap(std::string file_name,
                                               std::string symmetry,
                                               double length_units,
                                               double field_units)
       : SectorField(file_name), interpolator_m(NULL), symmetry_m(dipole),
         units_m(6, 1.), filename_m(file_name), phiOffset_m(0.) {
    units_m[0] *= length_units;
    units_m[1] *= length_units;
    units_m[2] *= length_units;
    units_m[3] *= field_units;
    units_m[4] *= field_units;
    units_m[5] *= field_units;
    setSymmetry(symmetry);
    if (_fields.find(file_name) == _fields.end()) {
        readMap();
        _fields[file_name] = new SectorMagneticFieldMap(*this);
    } else {
        SectorMagneticFieldMap* tgt = _fields[file_name];
        if (symmetry_m != tgt->symmetry_m || units_m != tgt->units_m ||
            filename_m != tgt->filename_m) {
            throw(LogicalError(
                    "SectorMagneticFieldMap::SectorMagneticFieldMap",
                    "Attempt to construct different SectorFieldMaps with same "+
                    std::string("file but different settings")
                  ));
        }
        setInterpolator(tgt->interpolator_m->clone());
    }
    phiOffset_m = interpolator_m->getMesh()->minZ();
}

SectorMagneticFieldMap::SectorMagneticFieldMap
                                        (const SectorMagneticFieldMap& field)
    : SectorField(field), interpolator_m(NULL), symmetry_m(field.symmetry_m),
      units_m(field.units_m),
      filename_m(field.filename_m), phiOffset_m(field.phiOffset_m) {
    Interpolator3dGridTo3d* interpolator = NULL;
    if (field.interpolator_m != NULL) {
        interpolator = field.interpolator_m->clone();
    }
    setInterpolator(interpolator);
}

SectorMagneticFieldMap::~SectorMagneticFieldMap() {
    delete interpolator_m;
}

Interpolator3dGridTo3d* SectorMagneticFieldMap::getInterpolator() {
    return interpolator_m;
}

void SectorMagneticFieldMap::setInterpolator
                                        (Interpolator3dGridTo3d* interpolator) {
    if (interpolator_m != NULL) {
        delete interpolator_m;
    }
    interpolator_m = interpolator;
    if (interpolator_m != NULL) {
        ThreeDGrid* grid = interpolator_m->getMesh();
        double phiTol = (grid->maxZ()-grid->minZ())/fractionalBBPhiTolerance_m;
        SectorField::setPolarBoundingBox
                                    (grid->minX(), grid->minY(), grid->minZ(),
                                     grid->maxX(), grid->maxY(), grid->maxZ(),
                                     0., 0., phiTol);
    }
    delete interpolator_m->clone(); // why?
}

std::string SectorMagneticFieldMap::getSymmetry() const {
    return SymmetryToString(symmetry_m);
}

void SectorMagneticFieldMap::setSymmetry(std::string name) {
    symmetry_m = StringToSymmetry(name);
}

void SectorMagneticFieldMap::readMap() {
    setInterpolator(IO::readMap(filename_m, units_m, symmetry_m));
}

void SectorMagneticFieldMap::freeMap() {
    setInterpolator(NULL);
}

SectorMagneticFieldMap::symmetry SectorMagneticFieldMap::StringToSymmetry
                                                             (std::string sym) {
    if (sym == "None") {
        return none;
    }
    if (sym == "Dipole") {
        return dipole;
    }
    throw(LogicalError(
            "SectorMagneticFieldMap::StringToSymmetry",
            "Didn't recognise symmetry type "+sym
          ));
}

std::string SectorMagneticFieldMap::SymmetryToString
                                        (SectorMagneticFieldMap::symmetry sym) {
    if (sym == none) {
        return "None";
    }
    if (sym == dipole) {
        return "Dipole";
    }
    throw(LogicalError(
            "SectorMagneticFieldMap::SymmetryToString",
            "Didn't recognise symmetry type "+sym
          ));
}

bool SectorMagneticFieldMap::getFieldstrengthPolar
                    (const Vector_t &R_p, Vector_t &E_p, Vector_t &B_p) const {
    // vector_t::operator[i] const returns by value, not by const reference
    // so we need to make an array here
    double R_temp[3] = {R_p[0], R_p[1], R_p[2]};
    R_temp[2] += phiOffset_m;
    // bounding box
    if (!isInBoundingBox(R_temp))
        return true;
    interpolator_m->function(R_temp, &B_p[0]);
    SectorField::convertToPolar(R_temp, &(B_p[0]));
    return false;
}

bool SectorMagneticFieldMap::getFieldstrength
                    (const Vector_t &R_c, Vector_t &E_c, Vector_t &B_c) const {
    // coordinate transform; field is in the x-z plane but OPAL-CYCL assumes
    // x-y plane; make sure to keep right handed coordinate system
    // I also do a rotation through 180 deg about x to make the bend in the
    // correct sense but this should be implemented in the end by OpalRing (i.e.
    // should be able to do off-midplane transformations)
    double R_temp[3] = {+R_c[1], -R_c[2], -R_c[0]};
    double B_temp[3] = {0., 0., 0.};
    SectorField::convertToPolar(R_temp);
    R_temp[2] += phiOffset_m;
    if (!isInBoundingBox(R_temp)) {
        // std::cerr << "SectorMagneticFieldMap R_c " << R_c << " R_temp "
        //           << R_temp[0] << " " << R_temp[1] << " " << R_temp[2]
        //           << std::endl;
        return true;
    }
    interpolator_m->function(R_temp, B_temp);
    // and here we transform back
    B_c(0) = -B_temp[2];
    B_c(1) = +B_temp[0];
    B_c(2) = -B_temp[1];
    // std::cerr << "\nSMFM::getFieldStrength R_c: " << R_c
    //          << " R_p: " << R_temp[0] << " " << R_temp[1] << " " << 180.*R_temp[2]/3.142
    //         << " B: " << B_c[0] << " " << B_c[1] << " " << B_c[2] << std::endl;
    return false;
}


void SectorMagneticFieldMap::clearFieldCache() {
    for (std::map<std::string, SectorMagneticFieldMap*>::iterator it =
                                   _fields.begin(); it != _fields.end(); ++it) {
        delete (*it).second;
    }
    _fields = std::map<std::string, SectorMagneticFieldMap*>();
}

void SectorMagneticFieldMap::getInfo(Inform* msg) {
   std::vector<double> bbmin = SectorField::getPolarBoundingBoxMin();
   std::vector<double> bbmax = SectorField::getPolarBoundingBoxMax();
   (*msg) << Filename_m << " (3D Sector Magnetostatic); "
          << "zini= " << bbmin[2] << " m; zfinal= " << bbmax[2] << " m;"
          << "rini= " << bbmin[0] << " m; rfinal= " << bbmax[0] << " m;"
          << endl;
}

void SectorMagneticFieldMap::print(std::ostream& out) {
   std::vector<double> bbmin = SectorField::getPolarBoundingBoxMin();
   std::vector<double> bbmax = SectorField::getPolarBoundingBoxMax();
   out << Filename_m << " (3D Sector Magnetostatic); "
       << "zini= " << bbmin[2] << " m; zfinal= " << bbmax[2] << " m;"
       << "rini= " << bbmin[0] << " m; rfinal= " << bbmax[0] << " m;"
       << std::endl;
}

bool SectorMagneticFieldMap::getFieldDerivative(const Vector_t &R, Vector_t &E,
                               Vector_t &B, const DiffDirection &dir) const {
    throw(LogicalError("SectorMagneticFieldMap::getFieldDerivative",
                       "Field map derivatives not implemented"));
}

void SectorMagneticFieldMap::swap() {}

const double SectorMagneticFieldMap::IO::floatTolerance_m = 1e-3;
const int SectorMagneticFieldMap::IO::sortOrder_m[3] = {0, 1, 2};

Interpolator3dGridTo3d* SectorMagneticFieldMap::IO::readMap
                          (std::string file_name, std::vector<double> units, 
                                        SectorMagneticFieldMap::symmetry sym) {
    try {
        INFOMSG("Opening sector field map " << file_name << endl);
        // get raw data
        std::vector< std::vector<double> > field_points = readLines
                                                             (file_name, units);
        // build grid
        ThreeDGrid* grid = generateGrid(field_points, sym);
        // build interpolator (convert grid to useful units)
        return getInterpolator(field_points, grid, sym);
    } catch(std::exception& exc) {
        throw(LogicalError(
                     "SectorMagneticFieldMap::IO::ReadMap",
                     "Failed to read file "+file_name+" with "+(&exc)->what()
               ));
    }
    return NULL;
}

Interpolator3dGridTo3d* SectorMagneticFieldMap::IO::getInterpolator
                        (const std::vector< std::vector<double> > field_points,
                         ThreeDGrid* grid,
                         SectorMagneticFieldMap::symmetry sym) {
    int y_start = 0;
    if (sym == SectorMagneticFieldMap::dipole) {
        y_start = (grid->ySize()+1)/2-1;
    }
    // build field arrays
    double *** Bx, *** By, *** Bz;
    int index = 0;
    Bx = new double**[grid->xSize()];
    By = new double**[grid->xSize()];
    Bz = new double**[grid->xSize()];
    for (int i = 0; i < grid->xSize(); ++i) {
        Bx[i] = new double*[grid->ySize()];
        By[i] = new double*[grid->ySize()];
        Bz[i] = new double*[grid->ySize()];
        for (int j = 0; j < grid->ySize(); ++j) {
            Bx[i][j] = new double[grid->zSize()];
            By[i][j] = new double[grid->zSize()];
            Bz[i][j] = new double[grid->zSize()];
            for (int k = 0; k < grid->zSize() && j >= y_start; ++k) {
                Bx[i][j][k] = field_points[index][3];
                By[i][j][k] = field_points[index][4];
                Bz[i][j][k] = field_points[index][5];
                ++index;
            }
        }
    }

    // extend field array downwards if field is dipole
    if (sym == SectorMagneticFieldMap::dipole) {
        for (int i = 0; i < grid->xSize(); ++i) {
            for (int j = 0; j < y_start; ++j) {
                for (int k = 0; k < grid->zSize(); ++k) {
                    Bx[i][j][k] = -Bx[i][grid->ySize()-j-1][k];
                    By[i][j][k] = By[i][grid->ySize()-j-1][k];
                    Bz[i][j][k] = -Bz[i][grid->ySize()-j-1][k];
                    ++index;
                }
            }
        }
    }

    return new Interpolator3dGridTo3d(grid, Bx, By, Bz);
}

bool SectorMagneticFieldMap::IO::comparator(std::vector<double> field_item1,
                     std::vector<double> field_item2) {
    const int* order = sortOrder_m;
    if (fabs(field_item1[order[0]] - field_item2[order[0]]) > floatTolerance_m) {
        return field_item1[order[0]] < field_item2[order[0]];
    }
    if (fabs(field_item1[order[1]] - field_item2[order[1]]) > floatTolerance_m) {
        return field_item1[order[1]] < field_item2[order[1]];
    }
    return field_item1[order[2]] < field_item2[order[2]];
}

std::vector< std::vector<double> > SectorMagneticFieldMap::IO::readLines
                             (std::string file_name, std::vector<double> units) {
    std::vector< std::vector<double> > field_points;
    std::string line;
    if (units.size() != 6) {
        throw(LogicalError(
                     "SectorMagneticFieldMap::IO::ReadLines",
                     "Units should be of length 6"
               ));
    }

    std::ifstream fin(file_name.c_str());
    if (!fin || !fin.is_open()) {
        throw(LogicalError(
                     "SectorMagneticFieldMap::IO::ReadLines",
                     "Failed to open file "+file_name
               ));
    }
    // skip header lines
    INFOMSG("Opened "+file_name << endl);
    for (size_t i = 0; i < 8; ++i) {
        std::getline(fin, line);
    }
    // read in field map
    int line_number = 0;
    while (fin) {
        std::vector<double> field(6);
        fin >> field[0] >> field[1] >> field[2] >> field[3] >> field[4]
            >> field[5];
        if (fin) {
            for (size_t i = 0; i < 6; ++i) {
                field[i] *= units[i];
            }
            field_points.push_back(field);
            line_number++;
        }
    }
    INFOMSG("Read " << line_number << " lines" << endl);

    // convert coordinates to polar; nb we leave field as cartesian
    for (size_t i = 0; i < field_points.size(); ++i) {
        SectorField::convertToPolar(&field_points[i][0]);
    }

    // force check sort order
    std::sort(field_points.begin(), field_points.end(),
                                        SectorMagneticFieldMap::IO::comparator);
    return field_points;
}

bool SectorMagneticFieldMap::IO::floatGreaterEqual(double in1, double in2) {
    in2 = in2*(1+floatTolerance_m)+floatTolerance_m;
    return in1 > in2;
}

ThreeDGrid* SectorMagneticFieldMap::IO::generateGrid
                       (const std::vector< std::vector<double> > field_points,
                        SectorMagneticFieldMap::symmetry sym) {
    std::vector<double> r_grid(1, field_points[0][0]);
    std::vector<double> y_grid(1, field_points[0][1]);
    std::vector<double> phi_grid(1, field_points[0][2]);
    for (size_t i = 0; i < field_points.size(); ++i) {
        if (floatGreaterEqual(field_points[i][0], r_grid.back())) {
            r_grid.push_back(field_points[i][0]);
        }
        if (floatGreaterEqual(field_points[i][1], y_grid.back())) {
            y_grid.push_back(field_points[i][1]);
        }
        if (floatGreaterEqual(field_points[i][2], phi_grid.back())) {
            phi_grid.push_back(field_points[i][2]);
        }
    }

    // reflect about y if symmetry is dipole
    int y_end = static_cast<int>(y_grid.size());
    if (sym == SectorMagneticFieldMap::dipole) {
        for (int i = 0; i < y_end-1; ++i) {
            y_grid.insert(y_grid.begin(), 2.*y_grid[i]-y_grid[2*i+1]);
        }
        y_end = y_grid.size();
    }
    INFOMSG("Grid size (r, y, phi) = ("
            << r_grid.size() << ", " << y_grid.size() << ", " << phi_grid.size()
            << ")" << endl);
    INFOMSG("Grid min (r [mm], y [mm], phi [rad]) = ("
            << r_grid[0] << ", " << y_grid[0] << ", " << phi_grid[0]
            << ")" << endl);
    INFOMSG("Grid max (r [mm], y [mm], phi [rad]) = ("
            << r_grid.back() << ", " << y_grid.back() << ", " << phi_grid.back()
            << ")" << endl);

    ThreeDGrid* grid = new ThreeDGrid(r_grid, y_grid, phi_grid);
    return grid;
}

