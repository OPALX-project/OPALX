//
// Class MonitorRep
//   Representation for an orbit position monitor.
//   The base class observes both planes.
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
#ifndef OPALX_MonitorRep_HH
#define OPALX_MonitorRep_HH

#include "AbsBeamline/Monitor.h"

class MonitorRep : public Monitor {
public:
    /// Constructor with given name.
    explicit MonitorRep(const std::string& name);

    MonitorRep();
    MonitorRep(const MonitorRep&);
    virtual ~MonitorRep();

private:
    // Not implemented.
    void operator=(const MonitorRep&);
};

#endif  // OPALX_MonitorRep_HH
