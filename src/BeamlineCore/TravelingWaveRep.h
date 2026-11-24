//
// Class TravelingWaveRep
//   Representation for a traveling wave.
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
#ifndef OPALX_TravelingWaveRep_HH
#define OPALX_TravelingWaveRep_HH

#include "AbsBeamline/TravelingWave.h"

class TravelingWaveRep : public TravelingWave {
public:
    /// Constructor with given name.
    explicit TravelingWaveRep(const std::string& name);

    TravelingWaveRep();
    TravelingWaveRep(const TravelingWaveRep&);
    virtual ~TravelingWaveRep();

private:
    // Not implemented.
    void operator=(const TravelingWaveRep&);
};

#endif  // OPALX_TravelingWaveRep_HH
