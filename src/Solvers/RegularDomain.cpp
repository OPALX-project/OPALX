//
// Class RegularDomain
//   Base class for simple domains that do not change the x-y shape in
//   longitudinal direction.
//
// Copyright (c) 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Solvers/RegularDomain.h"

RegularDomain::RegularDomain(const IntVector_t& nr,
                             const Vector_t& hr,
                             const std::string& interpl)
    : IrregularDomain(nr, hr, interpl)
    , nxy_m(nr[0] * nr[1])
{ }

void RegularDomain::resizeMesh(Vector_t& origin, Vector_t& hr, const Vector_t& rmin,
                                const Vector_t& rmax, double dh)
{
    Vector_t mymax = Vector_t(0.0, 0.0, 0.0);

    // apply bounding box increment dh, i.e., "BBOXINCR" input argument
    double zsize = rmax[2] - rmin[2];

    setMinMaxZ(rmin[2] - zsize * (1.0 + dh),
               rmax[2] + zsize * (1.0 + dh));

    origin = Vector_t(getXRangeMin(), getYRangeMin(), getMinZ());
    mymax  = Vector_t(getXRangeMax(), getYRangeMax(), getMaxZ());

    for (int i = 0; i < 3; ++i)
        hr[i] = (mymax[i] - origin[i]) / nr_m[i];
}