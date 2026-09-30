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
    void extractPotential(ScalarField2_t& potential);
    void determineField(const ScalarField2_t& potential, VectorField2_t& electricField);

private:
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
