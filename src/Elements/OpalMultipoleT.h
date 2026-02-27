/*
 *  Copyright (c) 2017, Titus Dascalu
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */


#ifndef OPAL_OPALMULTIPOLET_HH
#define OPAL_OPALMULTIPOLET_HH

#include "Elements/OpalElement.h"

class OpalMultipoleT final: public OpalElement {

public:
    // The attributes of class OpalMultipoleT
    enum {
        TP = COMMON,     // Transverse field components
        // Attributes for a straight multipole
        RFRINGE,         // Length of right fringe field
        LFRINGE,         // Length of left fringe field
        HAPERT,          // Aperture horizontal dimension
        VAPERT,          // Aperture vertical dimension
        MAXFORDER,       // Maximum order in the field expansion
        ROTATION,        // Rotation angle about central axis for skew elements
        EANGLE,          // Entrance angle
        BBLENGTH,        // Length within which field is calculated
        // Further attributes for a constant radius curved multipole
        ANGLE,           // Bending angle of a sector magnet
        MAXXORDER,       // Maximum order in x in polynomial expansions
        // Further attributes for a variable radius multipole
        VARRADIUS,       // Variable radius flag
        ENTRYOFFSET,     // Longitudinal offset from standard entrance point
        // Time dependence
        SCALING_MODEL,   // Name of a time dependence object
        // Misalignments
        MISALIGN_H,      // Horizontal misalignment [m]
        MISALIGN_V,      // Vertical misalignment [m]
        MISALIGN_S,      // Longitudinal misalignment [m]
        MISALIGN_ROLL,   // Roll misalignment [rad] about the longitudinal axis
        MISALIGN_PITCH,  // Pitch misalignment [rad] about the horizontal axis
        MISALIGN_YAW,    // Yaw misalignment [rad] about the vertical axis
        MISALIGN_AXISOFFSET, // Vertical offset of rotation axes [m]
        SIZE             // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalMultipoleT();

    /** Inherited copy constructor */
    OpalMultipoleT* clone(const std::string& name) override;

    /** Update the MultipoleT with new parameters from UI parser */
    void update() override;

    void print(std::ostream& os) const override;

private:
    // Not implemented.
    OpalMultipoleT(const OpalMultipoleT&) = delete;
    void operator=(const OpalMultipoleT&) = delete;

    /** Default value for maximum series expansion order */
    static constexpr double DefaultMAXFORDER = 3.0;
    static constexpr double MinimumMAXFORDER = 1.0;
    static constexpr double MaximumMAXFORDER = 20.0;
    static constexpr double DefaultMAXXORDER = 20.0;

    // Clone constructor.
    OpalMultipoleT(const std::string& name, OpalMultipoleT* parent);
};

#endif // OPAL_OpalMultipoleT_HH
