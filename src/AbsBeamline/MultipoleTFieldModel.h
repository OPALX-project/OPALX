// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_MultipoleTFieldModel_HH
#define OPALX_MultipoleTFieldModel_HH

#include <Kokkos_Core.hpp>
#include "VectorMath.h"

/**
 * @namespace MultipoleTFieldModel
 * @brief Stateless field math for the combined-function multipole with a tanh
 * fringe, callable inside Kokkos kernels.
 *
 * The field comes from the scalar potential
 * @f[
 *   V = f_0(x,s)\,y + f_1(x,s)\,\frac{y^3}{3!} + f_2(x,s)\,\frac{y^5}{5!} + \dots
 * @f]
 * in the local Frenet-Serret coordinates @f$(x, y, s)@f$ of the body, with
 * mid-plane symmetry. On the mid-plane
 * @f[
 *   B_y = f_0(x,s) = T(x)\,S(s),
 * @f]
 * where the transverse profile @f$T(x) = \sum_k b_k x^k@f$ is the polynomial the
 * user gives (@c TP: dipole, gradient, and so on) and @f$S(s)@f$ is the tanh
 * fringe
 * @f[
 *   S(s) = \frac{1}{2}\left[\tanh\frac{s + s_0}{\lambda_L}
 *                         - \tanh\frac{s - s_0}{\lambda_R}\right],
 *   \qquad s_0 = \frac{L}{2},
 * @f]
 * measured from the centre of the body. For a straight body the higher terms
 * follow
 * @f[
 *   f_n = (-1)^n \sum_{i=0}^{n} \binom{n}{i}\, T^{(2i)}\, S^{(2n-2i)},
 * @f]
 * and for a curved body the equivalent trinomial sum that carries the Frenet
 * scale factor @f$h_s = 1 + x/\rho@f$.
 *
 * @note With @f$\lambda = 0@f$ on a side that side is a hard edge: the tanh
 * becomes a step, so @f$S = 1@f$ inside the body and 0 outside with all
 * derivatives zero. This mirrors @c HGAP = 0 on an @c SBEND.
 */
namespace MultipoleTFieldModel {

    /// Number of transverse profile coefficients: dipole (b_0) to decapole (b_5).
    constexpr unsigned int NumPoles = 6;

    /// Highest supported order of the vertical expansion (MAXFORDER).
    constexpr unsigned int MaxFOrder = 9;

    /// Length of the derivative arrays; the recursions index up to 2 * MaxFOrder + 2.
    constexpr unsigned int MaxDerivatives = 20;

    /// Factorial table bound (and the highest power powers() may be asked for).
    constexpr unsigned int MaxFactorial = 20;

    /// Fringe reach in fringe lengths: tanh has fallen to 6e-6 of the flat top at 6 lambda.
    constexpr double FringeReach = 6.0;

    /// @brief Factorial of n for n <= MaxFactorial, on host and device.
    /// @note No bounds check; callers cap the order (see MaxFOrder).
    KOKKOS_INLINE_FUNCTION double factorial(const unsigned int n) {
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

    /// @brief x^n by repeated squaring, on host and device.
    KOKKOS_INLINE_FUNCTION double powerInteger(double x, unsigned int n) {
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

    /// @brief Fill powers[0 .. maxPower] with value^i.
    KOKKOS_INLINE_FUNCTION void powers(
            const double value, const unsigned int maxPower,
            Kokkos::Array<double, MaxDerivatives>& out) {
        out[0] = 1.0;
        for (unsigned int i = 1; i <= maxPower; ++i) {
            out[i] = out[i - 1] * value;
        }
    }

    /**
     * @brief The transverse profile T(x) and its first @p count - 1 derivatives.
     *
     * @param profile Coefficients b_k of T(x) = sum_k b_k x^k.
     * @param count Number of entries to fill (T, T', T'', ...).
     * @param x Transverse coordinate in the local Frenet frame.
     * @param out Filled with T^(i)(x) for i = 0 .. count - 1.
     */
    KOKKOS_INLINE_FUNCTION void transverseDerivatives(
            const Kokkos::Array<double, NumPoles>& profile, const unsigned int count,
            const double x, Kokkos::Array<double, MaxDerivatives>& out) {
        Kokkos::Array<double, NumPoles> coefficients = profile;
        for (unsigned int i = 0; i < count; ++i) {
            out[i] = 0.0;
            for (unsigned int j = 0; j < NumPoles; ++j) {
                out[i] += coefficients[j] * powerInteger(x, j);
            }
            // Differentiate the coefficient array in place for the next order.
            for (unsigned int j = 0; j < NumPoles - 1; ++j) {
                coefficients[j] = coefficients[j + 1] * static_cast<double>(j + 1);
            }
            coefficients[NumPoles - 1] = 0.0;
        }
    }

    /**
     * @brief Build the polynomial coefficients of the tanh derivatives.
     *
     * Row n holds the coefficients of the polynomial @f$P_n@f$ with
     * @f$\frac{d^n}{du^n}\tanh u = P_n(\tanh u)@f$, from @f$P_0(t) = t@f$ and
     * @f$P_n = (1 - t^2)\,P_{n-1}'@f$. Host only; the element copies the table to
     * the device once.
     *
     * @param numDerivatives Highest derivative order needed; the table gets
     *        @p numDerivatives + 1 rows and @p numDerivatives + 2 columns.
     * @param table Any host-accessible 2-D table, already sized, with
     *        @c operator()(row, column) and @c extent().
     */
    template <class HostTable>
    void tanhDerivativeTable(const unsigned int numDerivatives, HostTable& table) {
        const unsigned int numCoefficients = numDerivatives + 2;
        for (unsigned int n = 0; n <= numDerivatives; ++n) {
            for (unsigned int k = 0; k < numCoefficients; ++k) {
                table(n, k) = 0.0;
            }
        }
        // P0(t) = t
        table(0, 1) = 1.0;
        for (unsigned int n = 1; n <= numDerivatives; ++n) {
            for (unsigned int k = 0; k < numCoefficients; ++k) {
                double value = 0.0;
                if (k + 1 < numCoefficients) {
                    value += (k + 1) * table(n - 1, k + 1);
                }
                if (k >= 1) {
                    value -= (k - 1) * table(n - 1, k - 1);
                }
                table(n, k) = value;
            }
        }
    }

    /**
     * @brief The fringe S(s) and its derivatives, from the tanh coefficient table.
     *
     * @param s0 Half the body length; the two tanh steps sit at s = -s0 and s = +s0.
     * @param lambdaLeft Fringe length of the entrance side; 0 gives a hard edge there.
     * @param lambdaRight Fringe length of the exit side; 0 gives a hard edge there.
     * @param s Longitudinal coordinate measured from the centre of the body.
     * @param tanhTable Table from tanhDerivativeTable(); its extents set how many
     *        derivatives are produced.
     * @param out Filled with S^(i)(s) for i = 0 .. tanhTable.extent(0) - 1.
     */
    template <class TableView>
    KOKKOS_INLINE_FUNCTION void fringeDerivatives(
            const double s0, const double lambdaLeft, const double lambdaRight, const double s,
            const TableView& tanhTable, Kokkos::Array<double, MaxDerivatives>& out) {
        const bool hardLeft  = lambdaLeft <= 0.0;
        const bool hardRight = lambdaRight <= 0.0;
        // A hard edge is the step the tanh becomes as lambda -> 0: +-1 with no derivatives.
        const double tLeft  = hardLeft ? ((s + s0 >= 0.0) ? 1.0 : -1.0)
                                       : Kokkos::tanh((s + s0) / lambdaLeft);
        const double tRight = hardRight ? ((s - s0 >= 0.0) ? 1.0 : -1.0)
                                        : Kokkos::tanh((s - s0) / lambdaRight);

        const unsigned int numDerivatives  = tanhTable.extent(0);
        const unsigned int numCoefficients = tanhTable.extent(1);
        double lambdaLeftN                 = 1.0;
        double lambdaRightN                = 1.0;
        for (unsigned int i = 0; i < numDerivatives; ++i) {
            double tLeftN    = 1.0;
            double tRightN   = 1.0;
            double leftTerm  = 0.0;
            double rightTerm = 0.0;
            for (unsigned int j = 0; j < numCoefficients; ++j) {
                const double coefficient = tanhTable(i, j);
                leftTerm += coefficient * tLeftN;
                rightTerm += coefficient * tRightN;
                tLeftN *= tLeft;
                tRightN *= tRight;
            }
            // A hard-edge side contributes only its step value (i == 0), no derivatives.
            const double left  = (hardLeft && i > 0) ? 0.0 : leftTerm / lambdaLeftN;
            const double right = (hardRight && i > 0) ? 0.0 : rightTerm / lambdaRightN;
            out[i]             = (left - right) / 2.0;
            if (!hardLeft) {
                lambdaLeftN *= lambdaLeft;
            }
            if (!hardRight) {
                lambdaRightN *= lambdaRight;
            }
        }
    }

    /**
     * @brief Pure-value inputs for one field evaluation inside a Kokkos kernel.
     */
    struct FieldInputs {
        /// Transverse profile coefficients b_k [T m^-k], already multiplied by the time scale.
        Kokkos::Array<double, NumPoles> profile;
        /// Number of terms in the vertical expansion (MAXFORDER).
        unsigned int maxFOrder;
        /// Body length: the arc length for a curved body, the straight length otherwise.
        double bodyLength;
        /// Design-path curvature h = angle / L, signed; 0 for a straight body.
        double curvature;
        /// Entrance-side fringe length [m]; 0 = hard edge.
        double fringeLeft;
        /// Exit-side fringe length [m]; 0 = hard edge.
        double fringeRight;
    };

    /**
     * @brief Field of a straight combined-function multipole, in the body frame.
     *
     * @param x Transverse coordinate.
     * @param y Vertical coordinate.
     * @param s Longitudinal coordinate measured from the centre of the body.
     */
    template <class TableView>
    KOKKOS_INLINE_FUNCTION Vector_t<double, 3> straightField(
            const double x, const double y, const double s, const FieldInputs& in,
            const TableView& tanhTable) {
        Kokkos::Array<double, MaxDerivatives> dt{};
        Kokkos::Array<double, MaxDerivatives> ds{};
        transverseDerivatives(in.profile, in.maxFOrder * 2 + 2, x, dt);
        fringeDerivatives(
                0.5 * in.bodyLength, in.fringeLeft, in.fringeRight, s, tanhTable, ds);

        Vector_t<double, 3> B(0.0);
        for (unsigned int n = 0; n <= in.maxFOrder; ++n) {
            double sumX = 0.0;
            double sumY = 0.0;
            double sumS = 0.0;
            for (unsigned int i = 0; i <= n; ++i) {
                const double binomial = factorial(n) / (factorial(i) * factorial(n - i));
                sumX += binomial * dt[2 * i + 1] * ds[2 * n - 2 * i];
                sumY += binomial * dt[2 * i] * ds[2 * n - 2 * i];
                sumS += binomial * dt[2 * i] * ds[2 * n - 2 * i + 1];
            }
            const double sign  = powerInteger(-1.0, n);
            const double oddY  = powerInteger(y, 2 * n + 1) / factorial(2 * n + 1) * sign;
            const double evenY = powerInteger(y, 2 * n) / factorial(2 * n) * sign;
            B(0) += sumX * oddY;
            B(1) += sumY * evenY;
            B(2) += sumS * oddY;
        }
        return B;
    }

    /**
     * @brief Field of a curved combined-function multipole, in the local Frenet basis.
     *
     * The curvature is signed and @p x is the Frenet transverse coordinate that
     * continues the entrance frame's x, so @f$\rho = 1/h@f$ is signed too and
     * @f$h_s = 1 + x/\rho@f$ is positive on both sides of the design orbit.
     *
     * @param x Transverse coordinate in the local Frenet frame.
     * @param y Vertical coordinate.
     * @param s Arc length measured from the centre of the body.
     */
    template <class TableView>
    KOKKOS_INLINE_FUNCTION Vector_t<double, 3> curvedField(
            const double x, const double y, const double s, const FieldInputs& in,
            const TableView& tanhTable) {
        Kokkos::Array<double, MaxDerivatives> dt{};
        Kokkos::Array<double, MaxDerivatives> ds{};
        Kokkos::Array<double, MaxDerivatives> rhoPowers{};
        Kokkos::Array<double, MaxDerivatives> hsPowers{};
        transverseDerivatives(in.profile, in.maxFOrder * 2 + 2, x, dt);
        fringeDerivatives(
                0.5 * in.bodyLength, in.fringeLeft, in.fringeRight, s, tanhTable, ds);

        const double rho = 1.0 / in.curvature;
        const double hs  = 1.0 + x / rho;
        powers(rho, in.maxFOrder, rhoPowers);
        powers(hs, 2 * in.maxFOrder, hsPowers);

        Vector_t<double, 3> B(0.0);
        for (unsigned int n = 0; n <= in.maxFOrder; ++n) {
            double sumX = 0.0;
            double sumY = 0.0;
            double sumS = 0.0;
            for (unsigned int i = 0; i <= n; ++i) {
                for (unsigned int j = 0; j <= n - i; ++j) {
                    const double trinomial =
                            factorial(n) / (factorial(i) * factorial(j) * factorial(n - i - j));
                    const double scale = 1.0 / rhoPowers[i] / hsPowers[2 * n - i - 2 * j];
                    const double dtx =
                            dt[i + 2 * j + 1] - (2 * n - i - 2 * j) / rho / hs * dt[i + 2 * j];
                    sumX += trinomial * scale * dtx * ds[2 * n - 2 * i - 2 * j];
                    sumY += trinomial * scale * dt[i + 2 * j] * ds[2 * n - 2 * i - 2 * j];
                    sumS += trinomial * scale * dt[i + 2 * j] * ds[2 * n - 2 * i - 2 * j + 1];
                }
            }
            const double sign  = powerInteger(-1.0, n);
            const double oddY  = powerInteger(y, 2 * n + 1) / factorial(2 * n + 1) * sign;
            const double evenY = powerInteger(y, 2 * n) / factorial(2 * n) * sign;
            B(0) += sumX * oddY;
            B(1) += sumY * evenY;
            B(2) += sumS * oddY;
        }
        return B;
    }

    /**
     * @brief Field at a point in bend coordinates, in the local body basis.
     *
     * @param arc (radial offset, vertical, arc length) from
     *        GeometryHelper::toBendArcCoords(); the arc length runs from 0 at the
     *        entrance face to @c bodyLength at the exit face.
     * @param in Field inputs, curvature included.
     * @param tanhTable Table from tanhDerivativeTable().
     * @return @f$(B_x, B_y, B_s)@f$ in the local Frenet basis at that arc position.
     *         The caller rotates it into the entrance frame with
     *         GeometryHelper::rotateArcFieldToEntry().
     *
     * @note toBendArcCoords() measures the radial offset away from the centre of
     *       curvature for either sign of the curvature, while the Frenet x has to
     *       continue the entrance frame's x, so the sign of the curvature is
     *       applied here.
     */
    template <class TableView>
    KOKKOS_INLINE_FUNCTION Vector_t<double, 3> field(
            const Vector_t<double, 3>& arc, const FieldInputs& in, const TableView& tanhTable) {
        const double sFromCentre = arc(2) - 0.5 * in.bodyLength;
        if (Kokkos::abs(in.curvature) <= 1.0e-15) {
            return straightField(arc(0), arc(1), sFromCentre, in, tanhTable);
        }
        const double x = (in.curvature >= 0.0) ? arc(0) : -arc(0);
        return curvedField(x, arc(1), sFromCentre, in, tanhTable);
    }

}  // namespace MultipoleTFieldModel

#endif  // OPALX_MultipoleTFieldModel_HH
