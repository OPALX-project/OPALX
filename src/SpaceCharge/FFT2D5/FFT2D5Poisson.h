/*
 *  Copyright (c) 2025, Jon Thompson
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

#ifndef OPALX_FFT2D5SOLVER_H
#define OPALX_FFT2D5SOLVER_H

// clang-format off
#include "Ippl.h"
#include "FFT/FFT.h"
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
// clang-format on

class FFT2D5Poisson {
public:
    using Vector2_t       = ippl::Vector<double, 2>;
    using Mesh2_t         = ippl::UniformCartesian<double, 2>;
    using ScalarField2_t  = ippl::Field<double, 2, Mesh2_t, Cell>;
    using VectorField2_t  = ippl::Field<Vector2_t, 2, Mesh2_t, Cell>;
    using Layout2_t       = ippl::FieldLayout<2>;
    using Vector3_t       = ippl::Vector<double, 3>;
    using FFT_t           = ippl::FFT<ippl::RCTransform, ScalarField2_t>;
    using NDIndex2_t      = ippl::NDIndex<2>;
    using ComplexField2_t = FFT_t::ComplexField;

    // API
    FFT2D5Poisson();
    void solve(ScalarField2_t& chargeDensity, VectorField2_t& electricField);

    /**
     * @brief Free-space potential of a uniformly charged rectangular source cell.
     *
     * Returns the cell average of -log(r)/(2*pi), including its finite value at
     * the origin. The convolution separately multiplies by the source-cell area.
     * Positive hx and hy are required. The arbitrary reference length in log(r)
     * changes only the potential gauge, never the electric field.
     */
    KOKKOS_INLINE_FUNCTION static double cellAveragedGreensFunction(
            const double x, const double y, const double hx, const double hy) {
        // Scale first so physical units do not enter cancellation of the four
        // rectangle-corner primitives. Restore the logarithmic gauge afterward.
        const double scale = hx > hy ? hx : hy;
        const double u     = x / scale;
        const double v     = y / scale;
        const double a     = 0.5 * hx / scale;
        const double b     = 0.5 * hy / scale;
        double average;
        if (u * u + v * v >= 64.0) {
            // Away from the singularity, 4-point Gauss integration in each
            // direction avoids subtraction of nearly equal corner primitives.
            // At r >= 8*max(hx,hy), its log-kernel truncation error is < 1e-13.
            constexpr double nodes[4] = {
                    -0.8611363115940525752, -0.3399810435848562648, 0.3399810435848562648,
                    0.8611363115940525752};
            constexpr double weights[4] = {
                    0.3478548451374538574, 0.6521451548625461426, 0.6521451548625461426,
                    0.3478548451374538574};
            average = 0.0;
            for (unsigned i = 0; i < 4; ++i) {
                for (unsigned j = 0; j < 4; ++j) {
                    const double dx = u + a * nodes[i];
                    const double dy = v + b * nodes[j];
                    average += weights[i] * weights[j] * Kokkos::log(dx * dx + dy * dy) / 8.0;
                }
            }
        } else {
            average = (logPrimitive(u + a, v + b) - logPrimitive(u - a, v + b)
                       - logPrimitive(u + a, v - b) + logPrimitive(u - a, v - b))
                      / (4.0 * a * b);
        }
        return -(Kokkos::log(scale) + average) / (2.0 * Kokkos::numbers::pi_v<double>);
    }

    // Accessors for unit tests
    const ScalarField2_t& getGreensFn() const { return greensFn_m; }
    const ComplexField2_t& getGreensFnTr() const { return greensFnTr_m; }
    const ScalarField2_t& getDoubledRho() const { return doubledRho_m; }

    // Constants
    static constexpr size_t X = 0;
    static constexpr size_t Y = 1;
    static constexpr size_t Z = 2;

    // Helpers
    void initialiseViews(const ScalarField2_t& chargeDensity);
    void makeGreensFn();
    void placeChargeDensity(const ScalarField2_t& chargeDensity);
    void performConvolution();
    /// Copy the physical potential and its halo from the doubled free-space convolution.
    void extractPotential(ScalarField2_t& potential);
    /// Centered gradient; potential must include a valid one-cell halo (as extractPotential does).
    void determineField(const ScalarField2_t& potential, VectorField2_t& electricField);

private:
    // Mixed primitive of log(sqrt(x*x+y*y)), extended continuously to the axes.
    // Using atan2 on absolute coordinates preserves the odd symmetry in each
    // coordinate without branch jumps across the source-cell singularity.
    KOKKOS_INLINE_FUNCTION static double logPrimitive(const double x, const double y) {
        const double a = Kokkos::abs(x);
        const double b = Kokkos::abs(y);
        if (a == 0.0 || b == 0.0) {
            return 0.0;
        }
        const double value = 0.5 * a * b * Kokkos::log(a * a + b * b) - 1.5 * a * b
                             + 0.5 * a * a * Kokkos::atan2(b, a)
                             + 0.5 * b * b * Kokkos::atan2(a, b);
        return (x < 0.0) != (y < 0.0) ? -value : value;
    }

    size_t nx2_{};
    size_t ny2_{};
    Vector2_t dr_{};
    NDIndex2_t doubledDomain_m;
    NDIndex2_t complexDomain_m;
    Mesh2_t doubledMesh_m;
    Mesh2_t complexMesh_m;
    Layout2_t doubledLayout_m;
    Layout2_t complexLayout_m;
    ScalarField2_t doubledRho_m;
    ScalarField2_t greensFn_m;
    ComplexField2_t doubledRhoTr_m;
    ComplexField2_t greensFnTr_m;
    std::unique_ptr<FFT_t> fft_m;
    ippl::ParameterList params_m;
};

#endif  // OPALX_FFT2D5SOLVER_H
