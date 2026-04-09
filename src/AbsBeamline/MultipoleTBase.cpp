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

#include "MultipoleTBase.h"
#include "MultipoleT.h"

MultipoleTBase::MultipoleTBase(MultipoleT* element) : element_m(element) {}

void MultipoleTBase::generateTanhCoefficients(const unsigned int numDerivatives) {
    const auto numCoefficients = numDerivatives + 2;
    Kokkos::resize(tanhCoefficientsGpu_m, numDerivatives + 1, numCoefficients);
    tanhCoefficientsHost_m = Kokkos::create_mirror_view(tanhCoefficientsGpu_m);
    // Zero initialize
    for (unsigned int n = 0; n <= numDerivatives; ++n) {
        for (unsigned int k = 0; k < numDerivatives; ++k) {
            tanhCoefficientsHost_m(n, k) = 0.0;
        }
    }
    // 0th derivative: P0(t) = t
    tanhCoefficientsHost_m(0, 1) = 1.0;
    // Build higher derivatives iteratively
    for (unsigned int n = 1; n <= numDerivatives; ++n) {
        for (unsigned int k = 0; k < numCoefficients; ++k) {
            double val = 0.0;
            if (k + 1 < numCoefficients) {
                val += (k + 1) * tanhCoefficientsHost_m(n - 1, k + 1);
            }
            if (k >= 1) {
                val -= (k - 1) * tanhCoefficientsHost_m(n - 1, k - 1);
            }
            tanhCoefficientsHost_m(n, k) = val;
        }
    }
    Kokkos::deep_copy(tanhCoefficientsGpu_m, tanhCoefficientsHost_m);
}
