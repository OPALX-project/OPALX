//
// Class OpalQuadrupole
//   The QUADRUPOLE element.
//
//   A convenience wrapper around MULTIPOLE that accepts a single
//   quadrupole strength K1 instead of the full KN vector.
//   Internally the element is represented as a MultipoleRep with
//   KN = {0, K1}.
//
// Copyright (c) 2024 – 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_OpalQuadrupole_HH
#define OPAL_OpalQuadrupole_HH

#include "Elements/OpalElement.h"

class OpalQuadrupole : public OpalElement {
public:
    /// The attributes of class OpalQuadrupole.
    enum {
        K1 = COMMON,  // Normalised quadrupole strength in m^(-2).
        DK1,          // Normalised quadrupole strength error in m^(-2).
        K1S,          // Normalised skew quadrupole strength in m^(-2).
        DK1S,         // Normalised skew quadrupole strength error in m^(-2).
        SIZE
    };

    /// Exemplar constructor.
    OpalQuadrupole();

    virtual ~OpalQuadrupole();

    /// Make clone.
    virtual OpalQuadrupole* clone(const std::string& name);

    /// Print the object.
    virtual void print(std::ostream&) const;

    /// Update the embedded OPALX multipole.
    virtual void update();

private:
    // Not implemented.
    OpalQuadrupole(const OpalQuadrupole&);
    void operator=(const OpalQuadrupole&);

    // Clone constructor.
    OpalQuadrupole(const std::string& name, OpalQuadrupole* parent);
};

#endif  // OPAL_OpalQuadrupole_HH
