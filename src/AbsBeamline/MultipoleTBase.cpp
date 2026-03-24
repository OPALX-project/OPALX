/*
 *  Copyright (c) 2017, Titus Dascalu
 *  Copyright (c) 2018, Martin Duy Tat
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

#include "MultipoleTBase.h"
#include "MultipoleT.h"

MultipoleTBase::MultipoleTBase(MultipoleT* element) : element_m(element) {}

void MultipoleTBase::generateTanhCoefficients(const unsigned int numDerivatives) {
    const auto numCoefficients = numDerivatives + 2;
    Kokkos::resize(tanhCoefficientsDevice_m, numDerivatives + 1, numCoefficients);
    Kokkos::resize(tanhCoefficientsHost_m, numDerivatives + 1, numCoefficients);
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
    Kokkos::deep_copy(tanhCoefficientsDevice_m, tanhCoefficientsHost_m);
}
