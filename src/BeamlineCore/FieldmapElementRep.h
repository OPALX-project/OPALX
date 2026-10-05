//
// Class FieldmapElementRep
//   Representation for a field-map-driven element.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_FieldmapElementRep_HH
#define OPALX_FieldmapElementRep_HH

#include "AbsBeamline/FieldmapElement.h"
#include "BeamlineGeometry/Geometry.h"

class FieldmapElementRep : public FieldmapElement {
public:
    /// Constructor with given name.
    explicit FieldmapElementRep(const std::string& name);

    FieldmapElementRep();
    FieldmapElementRep(const FieldmapElementRep&);
    virtual ~FieldmapElementRep();

    /// Return an identical deep copy of the element.
    virtual ElementBase* clone() const override;

    /// Construct a read/write channel for the attribute [b]aKey[/b], or nullptr.
    virtual Channel* getChannel(const std::string& aKey, bool create = false) override;

    /// Get geometry. Version for non-constant object.
    virtual Geometry& getGeometry() override;

    /// Get geometry. Version for constant object.
    virtual const Geometry& getGeometry() const override;

private:
    void operator=(const FieldmapElementRep&) = delete;

    /// The geometry. Straight, with its length taken from the map in initialise().
    Geometry geometry;
};

#endif  // OPALX_FieldmapElementRep_HH
