//
// Class RectangularDomain
//   This class provides a rectangular beam pipe. The mesh adapts to the bunch
//   in longitudinal direction.
//
// Copyright (c) 2008,        Yves Ineichen, ETH Zürich,
//               2013 - 2015, Tülin Kaman, Paul Scherrer Institut, Villigen PSI, Switzerland
//               2017 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the master thesis
// "A Parallel Multigrid Solver for Beam Dynamics"
// and the paper
// "A fast parallel Poisson solver on irregular domains applied to beam dynamics simulations"
// (https://doi.org/10.1016/j.jcp.2010.02.022)
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
#ifdef HAVE_SAAMG_SOLVER

#include "Solvers/RectangularDomain.h"
#include "Utilities/OpalException.h"

RectangularDomain::RectangularDomain(double a, double b, IntVector_t nr, Vector_t hr)
    : RegularDomain(nr, hr, "CONSTANT")
{
    setRangeMin(Vector_t(-a, -b, getMinZ()));
    setRangeMax(Vector_t( a,  b, getMaxZ()));
}

void RectangularDomain::compute(Vector_t hr, NDIndex<3> /*localId*/){
    setHr(hr);
}

void RectangularDomain::constantInterpolation(int x, int y, int z, StencilValue_t& value,
                                              double &scaleFactor) const
{
    scaleFactor = 1.0;
    value.west  = -1.0 / (hr[0] * hr[0]);
    value.east  = -1.0 / (hr[0] * hr[0]);
    value.north = -1.0 / (hr[1] * hr[1]);
    value.south = -1.0 / (hr[1] * hr[1]);
    value.front = -1.0 / (hr[2] * hr[2]);
    value.back  = -1.0 / (hr[2] * hr[2]);

    value.center = 2.0 / (hr[0] * hr[0])
                 + 2.0 / (hr[1] * hr[1])
                 + 2.0 / (hr[2] * hr[2]);

    if (!isInside(x + 1, y, z))
        value.east = 0.0; //we are a right boundary point

    if (!isInside(x - 1, y, z))
        value.west = 0.0; //we are a left boundary point

    if (!isInside(x, y + 1, z))
        value.north = 0.0; //we are a upper boundary point

    if (!isInside(x, y - 1, z))
        value.south = 0.0; //we are a lower boundary point

    //dirichlet
    if (z == 0)
        value.front = 0.0;

    if (z == nr[2] - 1)
        value.back = 0.0;
}

#endif //#ifdef HAVE_SAAMG_SOLVER