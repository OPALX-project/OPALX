//
// Class BeamBeamRep
//   Geometric representation of a collective beam-beam element.
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
#ifndef OPALX_BeamBeamRep_HH
#define OPALX_BeamBeamRep_HH

#include "AbsBeamline/BeamBeam.h"
#include "BeamlineGeometry/Geometry.h"

/**
 * @brief Concrete BeamBeam placement with straight-line geometry.
 *
 * Owns the geometric length and supports cloning/channel access. Collective
 * fields are evaluated by the per-run BeamBeamInteraction, not by this geometry.
 */
class BeamBeamRep : public BeamBeam {
public:
    /// Constructor with given name.
    /// @param name Placed element name.
    explicit BeamBeamRep(const std::string& name);

    /// Construct an unnamed placement.
    BeamBeamRep();
    /// Copy the element and its geometry, not per-run interaction state.
    /// @param right Source placement.
    BeamBeamRep(const BeamBeamRep& right);
    virtual ~BeamBeamRep();

    /// @return An identical deep copy of the element; ownership transfers to the caller.
    ElementBase* clone() const override;

    /// Construct a read/write channel for aKey via the base element's channel interface.
    /// @param aKey Attribute name.
    /// @param create Whether the base interface may create a missing channel.
    /// @return A caller-owned channel, or nullptr if the attribute is unknown.
    Channel* getChannel(const std::string& aKey, bool create = false) override;

    /// @return Borrowed mutable geometry, valid for the element lifetime.
    Geometry& getGeometry() override;

    /// @return Borrowed immutable geometry, valid for the element lifetime.
    const Geometry& getGeometry() const override;

private:
    // Not implemented.
    void operator=(const BeamBeamRep&);

    /// The geometry.
    Geometry geometry;
};

#endif  // OPALX_BeamBeamRep_HH
