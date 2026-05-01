//
// Cubic Spline Interpolation to replace GSL spline
//
// Copyright (c) 2023, Paul Scherrer Institute, Villigen PSI, Switzerland
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

#ifndef ABSBEAMLINE_MULTIPOLET_CONFIG_H
#define ABSBEAMLINE_MULTIPOLET_CONFIG_H

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