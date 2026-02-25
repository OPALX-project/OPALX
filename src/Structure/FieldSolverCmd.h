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
#ifndef OPAL_FieldSolver_HH
#define OPAL_FieldSolver_HH

#include <string>
#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Attributes/Attributes.h"
#include "PartBunch/BCHandler.hpp"

#include "Ippl.h"

class BinningCmd;

enum class FieldSolverCmdType : short { NONE = -1, FFT = 0, OPEN = 1, CG = 2 };

// The attributes of class FieldSolverCmd.
namespace FIELDSOLVER {
    enum {
        TYPE,      // The field solver name
        BINS,      // Name of BINNING definition or NONE
        NX,        // mesh size in x
        NY,        // mesh size in y
        NZ,        // mesh size in z
        PARFFTX,   // parallelized grid in x
        PARFFTY,   // parallelized grid in y
        PARFFTZ,   // parallelized grid in z
        BCFFTX,    // boundary condition in x [FFT + AMR_MG only]
        BCFFTY,    // boundary condition in y [FFT + AMR_MG only]
        BCFFTZ,    // boundary condition in z [FFT + AMR_MG only]
        GREENSF,   // holds greensfunction to be used [FFT + P3M only]
        BBOXINCR,  // how much the boundingbox is increased
        SIZE
    };
}

class FieldSolverCmd : public Definition {
public:
    /// Exemplar constructor.
    FieldSolverCmd();

    virtual ~FieldSolverCmd();

    /// Make clone.
    virtual FieldSolverCmd* clone(const std::string& name);

    /// Find named FieldSolverCmd.
    static FieldSolverCmd* find(const std::string& name);

    std::string getType();
    std::string getBinsName() const;
    BinningCmd* getBinningCmd() const;

    /// Returns solver boundary conditions handler object.
    BCHandler<3> constructBCHandler() const;

    /// Return meshsize
    double getNX() const;

    /// Return meshsize
    double getNY() const;

    /// Return meshsize
    double getNZ() const;

    void setNX(double);

    void setNY(double);

    void setNZ(double);

    double getBoxIncr() const;

    /// Update the field solver data.
    virtual void update();

    /// Execute (init) the field solver data.
    virtual void execute();

    bool hasValidSolver();

    void setFieldSolverCmdType();
    FieldSolverCmdType getFieldSolverCmdType() const;

    ippl::Vector<bool, 3> getDomDec() const;

    Inform& printInfo(Inform& os) const;

private:
    // Not implemented.
    FieldSolverCmd(const FieldSolverCmd&);
    void operator=(const FieldSolverCmd&);

    // Clone constructor.
    FieldSolverCmd(const std::string& name, FieldSolverCmd* parent);

    std::string fsName_m;
    FieldSolverCmdType fsType_m;
};

inline FieldSolverCmdType FieldSolverCmd::getFieldSolverCmdType() const {
    return fsType_m;
}
inline ippl::Vector<bool, 3> FieldSolverCmd::getDomDec() const {
    return ippl::Vector<bool, 3>(
        Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTX]),
        Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTY]),
        Attributes::getBool(itsAttr[FIELDSOLVER::PARFFTZ]));
}

inline Inform& operator<<(Inform& os, const FieldSolverCmd& fs) {
    return fs.printInfo(os);
}

#endif  // OPAL_FieldSolverCmd_HH
