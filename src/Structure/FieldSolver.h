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
#include "Ippl.h"

enum class FieldSolverType : short { NONE = -1 };

class FieldSolver : public Definition {
public:
    /// Exemplar constructor.
    FieldSolver();

    virtual ~FieldSolver();

    /// Make clone.
    virtual FieldSolver* clone(const std::string& name);

    /// Find named FieldSolver.
    static FieldSolver* find(const std::string& name);

    std::string getType();

    /// Return meshsize
    double getMX() const;

    /// Return meshsize
    double getMY() const;

    /// Return meshsize
    double getMZ() const;

    /// Store emittance for mode 1.
    void setMX(double);

    /// Store emittance for mode 2.
    void setMY(double);

    /// Store emittance for mode 3.
    void setMZ(double);

    /// Update the field solver data.
    virtual void update();

    /// Execute (init) the field solver data.
    virtual void execute();

    bool hasValidSolver();

    void setFieldSolverType();
    FieldSolverType getFieldSolverType() const;

    Inform& printInfo(Inform& os) const;

private:
    // Not implemented.
    FieldSolver(const FieldSolver&);
    void operator=(const FieldSolver&);

    // Clone constructor.
    FieldSolver(const std::string& name, FieldSolver* parent);

    std::string fsName_m;
    FieldSolverType fsType_m;
};

inline FieldSolverType FieldSolver::getFieldSolverType() const {
    return fsType_m;
}

inline Inform& operator<<(Inform& os, const FieldSolver& fs) {
    return fs.printInfo(os);
}

#endif  // OPAL_FieldSolver_HH
