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

#include "MultipoleTStraight.h"
#include "MultipoleT.h"

MultipoleTStraight::MultipoleTStraight(MultipoleT* element) : MultipoleTBase(element) {}

void MultipoleTStraight::initialise() {
    straightGeometry_m.setElementLength(element_m->getLength());
    generateTanhCoefficients(element_m->getMaxFOrder() * 2 + 1);
}

void MultipoleTStraight::getField(
        const Kokkos::View<Vector_t<double, 3>*> R, Kokkos::View<Vector_t<double, 3>*> /*E*/,
        const Kokkos::View<Vector_t<double, 3>*> B, const double scaling, const size_t count) {
    // Local variables that are copied into the kernel
    const auto config           = element_m->getConfig();
    const auto tanhCoefficients = tanhCoefficientsGpu_m;
    // Kernel launch over all particles
    Kokkos::parallel_for(
            "MultipoleTStraight::getField()", count, KOKKOS_LAMBDA(const size_t i) {
                computeBField(R(i), B(i), scaling, config, tanhCoefficients);
            });
}

bool MultipoleTStraight::getField(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B,
        const double scaling) {
    return computeBField(R, B, scaling, element_m->getConfig(), tanhCoefficientsHost_m);
}
