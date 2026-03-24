/*
 *  Copyright (c) 2025, Jon Thompson
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

#ifndef OPALX_MULTIPOLETCONFIG_H
#define OPALX_MULTIPOLETCONFIG_H

#include <Kokkos_Core.hpp>

/** A structure containing the magnet configuration parameters */
struct MultipoleTConfig {
    /** Number of terms in z expansion used in calculating field components */
    unsigned int maxFOrder_m{3};
    /** Highest order of polynomial expansions in x */
    unsigned int maxXOrder_m{20};
    /** List of transverse profile coefficients */
    static constexpr unsigned int NumPoles = 6;
    Kokkos::Array<double, NumPoles> transverseProfile_m{};
    unsigned int transverseProfileMaxOrder_m{1};
    /** Magnet parameters */
    double length_m{1.0};
    double entranceAngle_m{0.0};
    double rotation_m{0.0};
    double bendAngle_m{0.0};
    bool variableRadius_m{false};
    double entryOffset_m{0.0};
    /** Define the zone in which the field is calculated */
    double verticalAperture_m{0.5};
    double horizontalAperture_m{0.5};
    double boundingBoxLength_m{0.0};
    /** Fringe field parameters **/
    double fringeS0_m;
    double fringeLambdaLeft_m;
    double fringeLambdaRight_m;
};

#endif  // OPALX_MULTIPOLETCONFIG_H