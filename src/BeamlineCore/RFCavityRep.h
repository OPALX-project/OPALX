//
// Class RFCavityRep
//   Representation for a RF cavity.
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
#ifndef OPALX_RFCavityRep_HH
#define OPALX_RFCavityRep_HH

#include "AbsBeamline/RFCavity.h"
#include "BeamlineGeometry/Geometry.h"

class RFCavityRep : public RFCavity {
public:
    /// Constructor with given name.
    explicit RFCavityRep(const std::string& name);

    RFCavityRep();
    RFCavityRep(const RFCavityRep&);
    virtual ~RFCavityRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase* clone() const;

    /// Get geometry.
    //  Return the element geometry.
    //  Version for non-constant object.
    virtual Geometry& getGeometry();

    /// Get geometry.
    //  Return the element geometry
    //  Version for constant object.
    virtual const Geometry& getGeometry() const;

private:
    // Not implemented.
    void operator=(const RFCavityRep&);

    /// The cavity's geometry.
    Geometry geometry;
};

#endif  // OPALX_RFCavityRep_HH
