//
// Class OpalSBend
//   The SBEND element.
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
#ifndef OPAL_OpalSBend_HH
#define OPAL_OpalSBend_HH

#include "Elements/OpalBend.h"

class OpalSBend : public OpalBend {
public:
    /// Exemplar constructor.
    OpalSBend();

    virtual ~OpalSBend();

    /// Make clone.
    virtual OpalSBend* clone(const std::string& name);

    /**
     * @brief Update the embedded CLASSIC sector bend.
     *
     * The OPALX-native `SBEND` follows historical OPAL normalization where
     * the input `L` is interpreted as the chord used to compute the design
     * radius:
     * \f[
     *   R = \frac{L}{2\sin(\theta/2)}, \qquad
     *   L_\mathrm{ref} = R\theta, \qquad
     *   k_0 = \frac{1}{R}.
     * \f]
     * This aligns the hard-edge body, Enge fringe support, and peak dipole
     * field with old OPAL for identical `SBEND, L, ANGLE` input.
     */
    virtual void update();

private:
    // Not implemented.
    OpalSBend(const OpalSBend&);
    void operator=(const OpalSBend&);

    // Clone constructor.
    OpalSBend(const std::string& name, OpalSBend* parent);
};

#endif  // OPAL_OpalSBend_HH
