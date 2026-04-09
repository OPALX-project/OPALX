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

#ifndef ABSBEAMLINE_MULTIPOLET_BASE_H
#define ABSBEAMLINE_MULTIPOLET_BASE_H

#include "PartBunch/PartBunch.h"

/** ---------------------------------------------------------------------
 *
 * MultipoleTBase is a base class for a straight or curved combined \n
 * function magnet (up to arbitrary multipole component) with fringe fields.\n
 *
 * ---------------------------------------------------------------------
 *
 * Class category: AbsBeamline \n
 * $Author: Titus Dascalu, Martin Duy Tat, Chris Rogers
 *
 * ---------------------------------------------------------------------
 *
 * The field is obtained from the scalar potential \n
 *     @f[ V = f_0(x,s) z + f_1 (x,s) \frac{z^3}{3!} + f_2 (x,s) \frac{z^5}{5!} +
 *     ...  @f] \n
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

#include "AbsBeamline/MultipoleTConfig.h"
class MultipoleT;
class BeamlineVisitor;
class BGeometryBase;

class MultipoleTBase {
public:
    /** Constructor */
    explicit MultipoleTBase(MultipoleT* element);
    /** Destructor */
    virtual ~MultipoleTBase() = default;

protected:
    MultipoleT* element_m;
    Kokkos::View<double**> tanhCoefficientsGpu_m;
    Kokkos::View<double**>::host_mirror_type tanhCoefficientsHost_m;

public:
    /** Initialise the element */
    virtual void initialise() = 0;
    /** Return the cell geometry */
    virtual BGeometryBase* getGeometry() = 0;
    /** Return the field for an array of points */
    virtual void getField(
            Kokkos::View<Vector_t<double, 3>*> /*R*/, Kokkos::View<Vector_t<double, 3>*> /*E*/,
            Kokkos::View<Vector_t<double, 3>*> /*B*/, double /*scaling*/, size_t /*count*/) = 0;

    /** Return the field for a single point */
    virtual bool getField(
            const Vector_t<double, 3>& /*R*/, Vector_t<double, 3>& /*E*/,
            Vector_t<double, 3>& /*B*/, double /*scaling*/) = 0;

    // Constants
    static constexpr size_t MaxFactorial         = 20;
    static constexpr size_t MaxPowerInteger      = 20;
    static constexpr unsigned int MaxDerivatives = 20;

    /** Helper function that returns factorial of n for n<=20 on both host and GPU */
    KOKKOS_INLINE_FUNCTION static double factorial(unsigned int n);
    /** Helper function that returns x^n for 0<=n<=20 on both host and GPU */
    KOKKOS_INLINE_FUNCTION static double powerInteger(double x, unsigned int n);
    /** Helper function that calculates transverse derivatives for multipole fields */
    KOKKOS_INLINE_FUNCTION static void calcTransverseDerivatives(
            const Kokkos::Array<double, MultipoleTConfig::NumPoles>& poles,
            unsigned int numDerivatives, double x,
            Kokkos::Array<double, MaxDerivatives>& derivatives);
    template <class ViewType>
    KOKKOS_INLINE_FUNCTION static void calcFringeDerivatives(
            const double& s0, const double& lambdaLeft, const double& lambdaRight, double s,
            const ViewType& tanhCoefficients, Kokkos::Array<double, MaxDerivatives>& derivatives);
    void generateTanhCoefficients(unsigned int numDerivatives);
    KOKKOS_INLINE_FUNCTION static Vector_t<double, 3> rotateFrame(
            const Vector_t<double, 3>& R, const MultipoleTConfig& config);
    KOKKOS_INLINE_FUNCTION static void calcPowers(
            double value, unsigned int maxPower, Kokkos::Array<double, MaxPowerInteger>& powers);
};

KOKKOS_INLINE_FUNCTION
double MultipoleTBase::factorial(const unsigned int n) {
    static constexpr double factorialTable[MaxFactorial + 1] = {
            1.0,
            1.0,
            2.0,
            6.0,
            24.0,
            120.0,
            720.0,
            5040.0,
            40320.0,
            362880.0,
            3628800.0,
            39916800.0,
            479001600.0,
            6227020800.0,
            87178291200.0,
            1307674368000.0,
            20922789888000.0,
            355687428096000.0,
            6402373705728000.0,
            121645100408832000.0,
            2432902008176640000.0};
    return factorialTable[n];
}

KOKKOS_INLINE_FUNCTION
double MultipoleTBase::powerInteger(double x, unsigned int n) {
    double result = 1.0;
    while (n > 0) {
        if (n & 1) {
            result *= x;
        }
        x *= x;
        n >>= 1;
    }
    return result;
}

KOKKOS_INLINE_FUNCTION
void MultipoleTBase::calcTransverseDerivatives(
        const Kokkos::Array<double, MultipoleTConfig::NumPoles>& poles,
        const unsigned int numDerivatives, const double x,
        Kokkos::Array<double, MaxDerivatives>& derivatives) {
    Kokkos::Array<double, MultipoleTConfig::NumPoles> coefficients = poles;
    for (unsigned int i = 0; i < numDerivatives; ++i) {
        // Calculate the value of this derivative
        derivatives[i] = 0;
        for (unsigned int j = 0; j < MultipoleTConfig::NumPoles; ++j) {
            derivatives[i] += coefficients[j] * powerInteger(x, j);
        }
        // Differentiate for the next derivative
        for (unsigned int j = 0; j < MultipoleTConfig::NumPoles - 1; ++j) {
            coefficients[j] = coefficients[j + 1] * static_cast<double>(j + 1);
        }
        coefficients[MultipoleTConfig::NumPoles - 1] = 0.0;
    }
}

template <class ViewType>
KOKKOS_INLINE_FUNCTION void MultipoleTBase::calcFringeDerivatives(
        const double& s0, const double& lambdaLeft, const double& lambdaRight, const double s,
        const ViewType& tanhCoefficients, Kokkos::Array<double, MaxDerivatives>& derivatives) {
    const double tLeft                 = Kokkos::tanh((s + s0) / lambdaLeft);
    const double tRight                = Kokkos::tanh((s - s0) / lambdaRight);
    double lambdaLeftN                 = 1.0;
    double lambdaRightN                = 1.0;
    const unsigned int numDerivatives  = tanhCoefficients.extent(0);
    const unsigned int numCoefficients = tanhCoefficients.extent(1);
    for (unsigned int i = 0; i < numDerivatives; ++i) {
        // Evaluate the polynomial for both left and right
        double tLeftN    = 1.0;
        double tRightN   = 1.0;
        double leftTerm  = 0.0;
        double rightTerm = 0.0;
        for (unsigned int j = 0; j < numCoefficients; ++j) {
            const auto coeff = tanhCoefficients(i, j);
            leftTerm += coeff * tLeftN;
            rightTerm += coeff * tRightN;
            tLeftN *= tLeft;
            tRightN *= tRight;
        }
        // Combine the left and right terms
        derivatives[i] = (leftTerm / lambdaLeftN - rightTerm / lambdaRightN) / 2;
        // Lambda powers for next derivative
        lambdaLeftN *= lambdaLeft;
        lambdaRightN *= lambdaRight;
    }
}

KOKKOS_INLINE_FUNCTION
Vector_t<double, 3> MultipoleTBase::rotateFrame(
        const Vector_t<double, 3>& R, const MultipoleTConfig& config) {
    Vector_t<double, 3> R1(3), R2(3);
    /** Apply two 2D rotation matrices to coordinate vector
     * Rotate around central axis => skew fields
     * Rotate azimuthally => entrance angle
     */
    // 1st rotation
    R1[0] = R[0] * std::cos(config.rotation_m) + R[1] * std::sin(config.rotation_m);
    R1[1] = -R[0] * std::sin(config.rotation_m) + R[1] * std::cos(config.rotation_m);
    R1[2] = R[2];
    // 2nd rotation
    R2[0] = R1[2] * std::sin(config.entranceAngle_m) + R1[0] * std::cos(config.entranceAngle_m);
    R2[1] = R1[1];
    R2[2] = R1[2] * std::cos(config.entranceAngle_m) - R1[0] * std::sin(config.entranceAngle_m);
    return R2;
}

KOKKOS_INLINE_FUNCTION void MultipoleTBase::calcPowers(
        const double value, const unsigned int maxPower,
        Kokkos::Array<double, MaxPowerInteger>& powers) {
    powers[0] = 1;
    for (unsigned int i = 1; i <= maxPower; i++) {
        powers[i] = powers[i - 1] * value;
    }
}

#endif
