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
#ifndef REGULAR_DOMAIN_H
#define REGULAR_DOMAIN_H

#include "Solvers/IrregularDomain.h"

class RegularDomain : public IrregularDomain {

public:

    RegularDomain(const IntVector_t<double, 3>& nr,
                  const Vector_t<double, 3>& hr,
                  const std::string& interpl);

    int getNumXY() const override {
        return nxy_m;
    }

    void setNumXY(int nxy) { nxy_m = nxy; }

    void resizeMesh(Vector_t<double, 3>& origin, Vector_t<double, 3>& hr, const Vector_t<double, 3>& rmin,
                    const Vector_t<double, 3>& rmax, double dh) override;

protected:
    /// function to handle the open boundary condition in longitudinal direction
    void robinBoundaryStencil(int z, double &F, double &B, double &C) const;

private:
    void constantInterpolation(int x, int y, int z, StencilValue_t& value,
                               double &scaleFactor) const override;

    /// number of nodes in the xy plane (for this case: independent of the z coordinate)
    int nxy_m;
};

#endif