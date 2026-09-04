//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef OPALX_MultipoleTRep_HH
#define OPALX_MultipoleTRep_HH

#include "AbsBeamline/MultipoleT.h"
#include "BeamlineGeometry/Geometry.h"

/**
 * @class MultipoleTRep
 * @brief Representation of a combined-function multipole with a tanh fringe.
 *
 * The representation owns the body geometry: a straight body, or a planar arc
 * when the deck gives a bend angle. The field itself lives in MultipoleT.
 */
class MultipoleTRep : public MultipoleT {
public:
    MultipoleTRep();
    explicit MultipoleTRep(const std::string& name);
    MultipoleTRep(const MultipoleTRep&);
    ~MultipoleTRep() override;

    ElementBase* clone() const override;
    Channel* getChannel(const std::string& aKey, bool create = false) override;

    Geometry& getGeometry() override;
    const Geometry& getGeometry() const override;

private:
    Geometry geometry_m{Geometry::makeStraight(0.0)};
};

#endif  // OPALX_MultipoleTRep_HH
