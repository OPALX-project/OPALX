//
// Class CollimatorRep
//   Representation for a transverse-aperture collimator.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPALX_CollimatorRep_HH
#define OPALX_CollimatorRep_HH

#include "AbsBeamline/Collimator.h"
#include "BeamlineGeometry/Geometry.h"

class CollimatorRep : public Collimator {
public:
    /// Constructor with given name.
    explicit CollimatorRep(const std::string& name);

    CollimatorRep();
    CollimatorRep(const CollimatorRep&);
    virtual ~CollimatorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase* clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns nullptr.
    virtual Channel* getChannel(const std::string& aKey, bool = false);

    /// Get geometry.
    //  Version for non-constant object.
    virtual Geometry& getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const Geometry& getGeometry() const;

private:
    // Not implemented.
    void operator=(const CollimatorRep&);

    /// The geometry.
    Geometry geometry;
};

#endif  // OPALX_CollimatorRep_HH
