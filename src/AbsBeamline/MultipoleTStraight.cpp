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

#include "MultipoleTStraight.h"
#include "MultipoleT.h"

MultipoleTStraight::MultipoleTStraight(MultipoleT* element) : MultipoleTBase(element) {}

void MultipoleTStraight::initialise() {
    straightGeometry_m.setElementLength(element_m->getLength());
    generateTanhCoefficients(element_m->getMaxFOrder() * 2 + 1);
}

void MultipoleTStraight::getField(
    const Kokkos::View<Vector_t<double, 3>*>& R, Kokkos::View<Vector_t<double, 3>*>& /*E*/,
    Kokkos::View<Vector_t<double, 3>*>& B, const double scaling) {
    // Local variables that are copied into the kernel
    const auto config           = element_m->getConfig();
    const auto tanhCoefficients = tanhCoefficientsDevice_m;
    // Kernel launch over all particles
    Kokkos::parallel_for(
        "MultipoleTStraight::getField()", ippl::getRangePolicy(R), KOKKOS_LAMBDA(const int i) {
            computeBField(R(i), B(i), scaling, config, tanhCoefficients);
        });
}

bool MultipoleTStraight::getField(
    const Vector_t<double, 3>& R, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B,
    const double scaling) {
    return computeBField(R, B, scaling, element_m->getConfig(), tanhCoefficientsHost_m);
}

Vector_t<double, 3> MultipoleTStraight::localCartesianToOpalCartesian(
    const Vector_t<double, 3>& r) {
    return {r[0], r[1], r[2] + element_m->getLength() / 2.0};
}
