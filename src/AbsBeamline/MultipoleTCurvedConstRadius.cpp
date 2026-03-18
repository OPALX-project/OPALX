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


#include "MultipoleTCurvedConstRadius.h"
#include "MultipoleT.h"

MultipoleTCurvedConstRadius::MultipoleTCurvedConstRadius(MultipoleT* element) :
    MultipoleTBase(element),
    planarArcGeometry_m(1.0, 1.0) {
}

void MultipoleTCurvedConstRadius::initialise() {
    planarArcGeometry_m.setElementLength(element_m->getLength());
    planarArcGeometry_m.setCurvature(element_m->getBendAngle() / element_m->getLength());
    generateTanhCoefficients(element_m->getMaxFOrder() * 2 + 1);
}

void MultipoleTCurvedConstRadius::getField(const Kokkos::View<Vector_t<double, 3>*>& R,
        Kokkos::View<Vector_t<double, 3>*>& /*E*/, Kokkos::View<Vector_t<double, 3>*>& B,
        const double scaling) {
    // Local variables that are copied into the kernel
    const auto config = element_m->getConfig();
    const auto tanhCoefficients = tanhCoefficientsDevice_m;
    // Kernel launch over all particles
    Kokkos::parallel_for("MultipoleTCurvedConstRadius::getField()", ippl::getRangePolicy(R),
            KOKKOS_LAMBDA(const int i) {
                computeBField(R(i), B(i), scaling, config, tanhCoefficients);
            });
}

void MultipoleTCurvedConstRadius::getField(const Vector_t<double, 3>& R,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B,
        const double scaling) {
    computeBField(R, B, scaling, element_m->getConfig(), tanhCoefficientsHost_m);
}

void MultipoleTCurvedConstRadius::transformCoords(Vector_t<double, 3>& R) {
    if(element_m->getBendAngle() != 0.0) {
        const double radius = element_m->getLength() / element_m->getBendAngle();
        const double alpha = std::atan(R[2] / (R[0] + radius));
        if(alpha != 0.0) {
            R[0] = R[2] / std::sin(alpha) - radius;
            R[2] = radius * alpha;
        }
    }
    R[2] += element_m->getLength() / 2.0; // Magnet origin at the center rather than entry
}

void MultipoleTCurvedConstRadius::transformBField(Vector_t<double, 3>& B, const Vector_t<double, 3>& R) {
    const double theta = R[2] * element_m->getBendAngle() / element_m->getLength();
    const double Bx = B[0];
    const double Bs = B[2];
    B[0] = Bx * std::cos(theta) - Bs * std::sin(theta);
    B[2] = Bx * std::sin(theta) + Bs * std::cos(theta);
}

Vector_t<double, 3> MultipoleTCurvedConstRadius::localCartesianToOpalCartesian(
        const Vector_t<double, 3>& r) {
    Vector_t<double, 3> result = r;
    if(element_m->getBendAngle() != 0.0) {
        const double radius = element_m->getLength() / element_m->getBendAngle();
        const double theta = element_m->getBendAngle() / 2.0;
        const double ds = radius * std::sin(theta);
        const double dx = radius * (1 - std::cos(theta));
        result[0] = -dx;
        result[2] = ds;
    }
    return result;
}

double MultipoleTCurvedConstRadius::getScaleFactor(const double x, double /*s*/) {
    return 1 + x * element_m->getBendAngle() / element_m->getLength();
}

double MultipoleTCurvedConstRadius::getFn(unsigned int /*n*/, double /*x*/, double /*s*/) {
    return 0.0;
}
