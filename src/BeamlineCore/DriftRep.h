//
// Class DriftRep
//   Representation for a drift space.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPALX_DriftRep_HH
#define OPALX_DriftRep_HH

#include "AbsBeamline/Drift.h"

class DriftRep : public Drift {
public:
    /// Constructor with given name.
    explicit DriftRep(const std::string& name);

    DriftRep();
    DriftRep(const DriftRep&);
    virtual ~DriftRep();

private:
    // Not implemented.
    void operator=(const DriftRep&);
};

#endif  // OPALX_DriftRep_HH
