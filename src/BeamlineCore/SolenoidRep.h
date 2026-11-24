//
// Class SolenoidRep
//   Representation for a solenoid magnet.
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
#ifndef OPALX_SolenoidRep_HH
#define OPALX_SolenoidRep_HH

#include "AbsBeamline/Solenoid.h"

class SolenoidRep : public Solenoid {
public:
    /// Constructor with given name.
    explicit SolenoidRep(const std::string& name);

    SolenoidRep();
    SolenoidRep(const SolenoidRep&);
    virtual ~SolenoidRep();

private:
    // Not implemented.
    void operator=(const SolenoidRep&);
};

#endif  // OPALX_SolenoidRep_HH
