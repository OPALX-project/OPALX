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

#ifndef OPAL_OPALMULTIPOLET_HH
#define OPAL_OPALMULTIPOLET_HH

#include "Elements/OpalElement.h"

/**
 * @class OpalMultipoleT
 * @brief Parser element for a combined-function multipole with a tanh fringe.
 *
 * @c TP gives the mid-plane field profile in tesla (dipole first, then the
 * gradient, and so on), sign included. A non-zero @c ANGLE makes the body a
 * planar arc of that bend angle, with @c L the arc length; the deck has to keep
 * the sign of @c TP[0] consistent with the bend direction, exactly as it would
 * for the field of any other magnet. Everything else about placement, apertures
 * and misalignment comes from OpalElement.
 */
class OpalMultipoleT final : public OpalElement {
public:
    // The attributes of class OpalMultipoleT
    enum {
        TP = COMMON,  // Transverse field profile
        LFRINGE,      // Length of the left (entrance) end field
        RFRINGE,      // Length of the right (exit) end field
        MAXFORDER,    // Number of terms in the vertical field expansion
        ANGLE,        // Bend angle of a curved magnet
        HAPERT,       // Aperture width
        VAPERT,       // Aperture height
        SCALING_MODEL,  // Name of a time dependence scaling the field
        SIZE            // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalMultipoleT();

    /** Inherited copy constructor */
    OpalMultipoleT* clone(const std::string& name) override;

    /** Destructor, nothing special */
    ~OpalMultipoleT() override = default;

    /** Update the MultipoleT with new parameters from UI parser */
    void update() override;

    void print(std::ostream& os) const override;

    // Not implemented.
    OpalMultipoleT(const OpalMultipoleT&) = delete;
    void operator=(const OpalMultipoleT&) = delete;

private:
    /** Constants */
    static constexpr double DefaultMAXFORDER = 3.0;

    // Clone constructor.
    OpalMultipoleT(const std::string& name, OpalMultipoleT* parent);
};

#endif  // OPAL_OpalMultipoleT_HH
