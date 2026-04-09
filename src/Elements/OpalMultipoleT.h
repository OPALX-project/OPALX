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

class OpalMultipoleT final : public OpalElement {
public:
    // The attributes of class OpalMultipoleT
    enum {
        TP = COMMON,  // Transverse field components
        // Attributes for a straight multipole
        RFRINGE,    // Length of right fringe field
        LFRINGE,    // Length of left fringe field
        HAPERT,     // Aperture horizontal dimension
        VAPERT,     // Aperture vertical dimension
        MAXFORDER,  // Maximum order in the field expansion
        ROTATION,   // Rotation angle about central axis for skew elements
        EANGLE,     // Entrance angle
        BBLENGTH,   // Length within which field is calculated
        // Further attributes for a constant radius curved multipole
        ANGLE,      // Bending angle of a sector magnet
        MAXXORDER,  // Maximum order in x in polynomial expansions
        // Further attributes for a variable radius multipole
        VARRADIUS,    // Variable radius flag
        ENTRYOFFSET,  // Longitudinal offset from standard entrance point
        // Time dependence
        SCALING_MODEL,  // Name of a time dependence object
        // Misalignments
        MISALIGN_H,           // Horizontal misalignment [m]
        MISALIGN_V,           // Vertical misalignment [m]
        MISALIGN_S,           // Longitudinal misalignment [m]
        MISALIGN_ROLL,        // Roll misalignment [rad] about the longitudinal axis
        MISALIGN_PITCH,       // Pitch misalignment [rad] about the horizontal axis
        MISALIGN_YAW,         // Yaw misalignment [rad] about the vertical axis
        MISALIGN_AXISOFFSET,  // Vertical offset of rotation axes [m]
        SIZE                  // size of the enum
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
    static constexpr double MinimumMAXFORDER = 1.0;
    static constexpr double MaximumMAXFORDER = 9.0;
    static constexpr double DefaultMAXXORDER = 20.0;

    // Clone constructor.
    OpalMultipoleT(const std::string& name, OpalMultipoleT* parent);
};

#endif  // OPAL_OpalMultipoleT_HH
