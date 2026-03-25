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


#ifndef CLASSIC_MULTIPOLET_CURVED_CONST_RADIUS_H
#define CLASSIC_MULTIPOLET_CURVED_CONST_RADIUS_H

/** ---------------------------------------------------------------------
  *
  * MultipoleT defines a curved combined function magnet with constant radius\n
  * of curvature, (up to arbitrary multipole component) with fringe fields
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * $Author: Titus Dascalu, Martin Duy Tat, Chris Rogers
  *
  * ---------------------------------------------------------------------
  *
  * The field is obtained from the scalar potential \n
  *     @f[ V = f_0(x,s) z + f_1 (x,s) \frac{z^3}{3!} + f_2 (x,s) \frac{z^5}{5!}
  *      + ...  @f] \n
  *     (x,z,s) -> Frenet-Serret local coordinates along the magnet \n
  *     z -> vertical component \n
  *     assume mid-plane symmetry \n 
  *     set field on mid-plane -> @f$ B_z = f_0(x,s) = T(x) \cdot S(s) @f$ \n
  *     T(x) -> transverse profile; this is a polynomial describing
  *             the field expansion on the mid-plane inside the magnet
  *             (not in the fringe field);
  *             1st term is the dipole strength, 2nd term is the 
  *             quadrupole gradient * x, etc. \n
  *          -> when setting the magnet, one gives the multipole
  *             coefficients of this polynomial (i.e. dipole strength,  
  *             quadrupole gradient, etc.) \n
  * \n
  * ------------- example ----------------------------------------------- \n
  *     Setting a combined function magnet with dipole, quadrupole and 
  *     sextupole components: \n
  *     @f$ T(x) = B_0 + B_1 \cdot x + B_2 \cdot x^2 @f$\n
  *     user gives @f$ B_0, B_1, B_2 @f$ \n
  * ------------- example end ------------------------------------------- \n
  * \n
  *     S(s) -> fringe field \n
  *     recursion -> @f$ f_n (x,s) = (-1)^n \cdot \sum_{i=0}^{n} C_n^i 
  *     \cdot T^{(2i)} \cdot S^{(2n-2i)} @f$ \n
  *     for curved magnets the above recursion is more complicated \n
  *     @f$ C_n^i @f$ -> binomial coeff; 
  *     @f$ T^{(n)} @f$ -> n-th derivative
  *
  * ---------------------------------------------------------------------
  */

#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "AbsBeamline/MultipoleTBase.h"

class MultipoleTCurvedConstRadius final : public MultipoleTBase {
public:
    /** Constructor */
    explicit MultipoleTCurvedConstRadius(MultipoleT* element);
    /** Initialise the element */
    void initialise() override;
    /** Return the cell geometry */
    BGeometryBase* getGeometry() override { return &planarArcGeometry_m; }
    /** Return the cell geometry */
    const BGeometryBase* getGeometry() const override { return &planarArcGeometry_m; }
    /** Return the field for an array of points */
    void getField(const Kokkos::View<Vector_t<double, 3>*>& R,
            Kokkos::View<Vector_t<double, 3>*>& E, Kokkos::View<Vector_t<double, 3>*>& B,
            double scaling, size_t count) override;
    bool getField(const Vector_t<double, 3>& R,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B,
            double scaling) override;
    Vector_t<double, 3> localCartesianToOpalCartesian(const Vector_t<double, 3>& r) override;

private:
    /** Geometry */
    PlanarArcGeometry planarArcGeometry_m;

    // Helpers
    KOKKOS_INLINE_FUNCTION static Vector_t<double, 3> toMagnetCoords(const Vector_t<double, 3>& R,
            const MultipoleTConfig& config);
    template<class ViewType>
    KOKKOS_INLINE_FUNCTION static bool computeBField(const Vector_t<double, 3>& R,
            Vector_t<double, 3>& B, double scaling, const MultipoleTConfig& config,
            const ViewType& tanhCoefficients);
};

KOKKOS_INLINE_FUNCTION Vector_t<double, 3> MultipoleTCurvedConstRadius::toMagnetCoords(
        const Vector_t<double, 3>& R, const MultipoleTConfig& config) {
    // Skew and entry angle
    Vector_t<double, 3> result = rotateFrame(R, config);
    // Go to local Frenet-Serret coordinates
    if(config.bendAngle_m != 0.0) {
        const double radius = config.length_m / config.bendAngle_m;
        const double rMinusX = radius - R[0];
        const double alpha = Kokkos::hypot(rMinusX, R[2]);
        result[0] = alpha - radius;
        result[1] = R[1];
        result[2] = radius * Kokkos::atan2(R[2], rMinusX);
    }
    // Magnet origin at the center rather than entry
    result[2] -= config.length_m / 2.0;
    return result;
}

template<class ViewType>
KOKKOS_INLINE_FUNCTION bool MultipoleTCurvedConstRadius::computeBField(const Vector_t<double, 3>& R,
        Vector_t<double, 3>& B, const double scaling, const MultipoleTConfig& config,
        const ViewType& tanhCoefficients) {
    const Vector_t<double, 3> RPrime = toMagnetCoords(R, config);
    const bool insideAperture =
            Kokkos::abs(RPrime[1]) <= config.verticalAperture_m / 2.0 &&
            Kokkos::abs(RPrime[0]) <= config.horizontalAperture_m / 2.0;
    const bool insideBoundingBox = config.boundingBoxLength_m == 0.0 ||
            Kokkos::abs(RPrime[2]) <= config.boundingBoxLength_m / 2.0;
    if(insideAperture && insideBoundingBox) {
        Vector_t<double, 3> myB{};
        Kokkos::Array<double, MaxDerivatives> dt{};
        Kokkos::Array<double, MaxDerivatives> ds{};
        Kokkos::Array<double, MaxPowerInteger> rhoPowers{};
        Kokkos::Array<double, MaxPowerInteger> hsPowers{};
        calcTransverseDerivatives(config.transverseProfile_m, config.maxFOrder_m * 2 + 1, RPrime[0],
                dt);
        calcFringeDerivatives(config.fringeS0_m, config.fringeLambdaLeft_m,
                config.fringeLambdaRight_m,
                RPrime[2], tanhCoefficients, ds);
        const double rho = config.length_m / config.bendAngle_m;
        calcPowers(rho, config.maxFOrder_m, rhoPowers);
        const double hs = 1 + RPrime[0] / rho;
        calcPowers(hs, 2 * config.maxFOrder_m, hsPowers);
        for(unsigned int n = 0; n <= config.maxFOrder_m; n++) {
            double innerSumX{};
            double innerSumZ{};
            double innerSumS{};
            for(unsigned int i = 0; i <= n; i++) {
                for(unsigned int j = 0; j <= n - i; j++) {
                    const double k = factorial(n) /
                            (factorial(i) * factorial(j) * factorial(n - i - j));
                    const double rhoHs = 1.0 / rhoPowers[i] / hsPowers[2 * n - i - 2 * j];
                    const double dtx = dt[i + 2 * j + 1] -
                            (2 * n - i - 2 * j) / rho / hs * dt[i + 2 * j];
                    innerSumX += k * rhoHs * dtx * ds[2 * n - 2 * i - 2 * j];
                    innerSumZ += k * rhoHs * dt[i + 2 * j] * ds[2 * n - 2 * i - 2 * j];
                    innerSumS += k * rhoHs * dt[i + 2 * j] * ds[2 * n - 2 * i - 2 * j + 1];
                }
            }
            const double negOnePowN = powerInteger(-1.0, n);
            const double xszk = powerInteger(RPrime[1], 2 * n + 1) / factorial(2 * n + 1) *
                    negOnePowN;
            const double zzk = powerInteger(RPrime[1], 2 * n) / factorial(2 * n) * negOnePowN;
            myB[0] += innerSumX * xszk;
            myB[1] += innerSumZ * zzk;
            myB[2] += innerSumS * xszk;
        }
        B += myB * scaling;
    }
    return !insideAperture;
}

#endif
