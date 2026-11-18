//
// Class OpalCavity
//   The RFCAVITY element.
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
#ifndef OPAL_OpalCavity_HH
#define OPAL_OpalCavity_HH

#include <memory>

#include "Elements/OpalElement.h"

class OpalWake;

class OpalCavity : public OpalElement {
public:
    /// The attributes of class OpalCavity.
    enum {
        VOLT = COMMON,    // Peak field [MV/m].
        DVOLT,            // The peak field error [MV/m].
        GEOMETRY,         // geometry of boundary
        FREQ,             // The RF frequency.
        LAG,              // The phase lag.
        DLAG,             // The phase lag error.
        FMAPFN,           // The filename of the fieldmap
        FAST,             // Faster but less accurate
        APVETO,           // Do not use this cavity in the Autophase procedure
        DESIGNENERGY,     // The mean kinetic energy at exit
        PHASE_MODEL,      // time dependent phase
        AMPLITUDE_MODEL,  // time dependent amplitude
        FREQUENCY_MODEL,  // time dependent frequency
        SIZE
    };

    /// Exemplar constructor.
    OpalCavity();

    virtual ~OpalCavity();

    /// Build the runtime element of one occurrence from its attribute values.
    static std::shared_ptr<ElementBase> build(const ObjectValues& values, double p0);

private:
    /// Set the fields of @p element from @p values.
    static void fill(ElementBase& element, const ObjectValues& values);

    // Not implemented.
    OpalCavity(const OpalCavity&);
    void operator=(const OpalCavity&);
};

#endif  // OPAL_OpalCavity_HH
