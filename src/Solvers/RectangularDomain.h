//
// Class RectangularDomain
//   :FIXME: add brief description
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
#ifndef RECTANGULAR_DOMAIN_H
#define RECTANGULAR_DOMAIN_H

#include <vector>
#include <string>
#include "IrregularDomain.h"

class RectangularDomain : public IrregularDomain {

public:
    using IrregularDomain::StencilIndex_t;
    using IrregularDomain::StencilValue_t;

    /**
     * \param a is the longer side a of the rectangle
     * \param b is the shorter side b of the rectangle
     *
     */
    RectangularDomain(double a, double b, Vector_t nr, Vector_t hr);

    /// calculates intersection with the beam pipe
    void compute(Vector_t hr, NDIndex<3> /*localId*/);

    /// returns number of nodes in xy plane (here independent of z coordinate)
    int getNumXY(int z);

    /// returns index of neighbours at (x,y,z)
    using IrregularDomain::getNeighbours;
    /// queries if a given (x,y,z) coordinate lies inside the domain
    inline bool isInside(int x, int y, int /*z*/) {
        double xx = (x - (nr_m[0] - 1) / 2.0) * hr_m[0];
        double yy = (y - (nr_m[1] - 1) / 2.0) * hr_m[1];
        return (xx <= getXRangeMax() && yy < getYRangeMax());
    }

private:
    /// number of nodes in the xy plane (for this case: independent of the z coordinate)
    int nxy_m;

    /// conversion from (x,y,z) to index on the 3D grid
    inline int getIdx(int x, int y, int z) {
        if(isInside(x, y, z) && x >= 0 && y >= 0 && z >= 0)
            return y * nr_m[0] + x + z * nxy_m;
        else
            return -1;
    }
    /// conversion from a 3D index to (x,y,z)
    inline void getCoord(int idx, int &x, int &y, int &z) {
        int ixy = idx % nxy_m;
        int inr = nr_m[0];
        x = ixy % inr;
        y = (ixy - x) / nr_m[0];
        z = (idx - ixy) / nxy_m;
    }

    void constantInterpolation(int x, int y, int z, StencilValue_t& value,
                               double &scaleFactor) override;
};

#endif

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:

