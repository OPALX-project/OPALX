//
// Class MultipoleRep
//   Representation for a general multipole.
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
#ifndef OPALX_MultipoleRep_HH
#define OPALX_MultipoleRep_HH

#include "AbsBeamline/Multipole.h"
#include "BeamlineGeometry/Geometry.h"

class MultipoleRep : public Multipole {
public:
    /// Constructor with given name.
    explicit MultipoleRep(const std::string& name);

    MultipoleRep();
    MultipoleRep(const MultipoleRep&);
    virtual ~MultipoleRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase* clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b]r and returns it.
    //  If the attribute does not exist, it returns nullptr.
    virtual Channel* getChannel(const std::string& aKey, bool = false);

    /// Get geometry.
    //  Return the element geometry.
    //  Version for non-constant object.
    virtual Geometry& getGeometry();

    /// Get geometry.
    //  Return the element geometry
    //  Version for constant object.
    virtual const Geometry& getGeometry() const;

private:
    /// Multipole geometry.
    Geometry geometry;

    // Not implemented.
    void operator=(const MultipoleRep&);
};

#endif  // OPALX_MultipoleRep_HH
