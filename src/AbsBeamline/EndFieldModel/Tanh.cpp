/*
 *  Copyright (c) 2017, Chris Rogers
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

#include <cmath>

#include "Utilities/GSLCompat.h"
#include "Utilities/GeneralOpalException.h"

#include "AbsBeamline/EndFieldModel/CompactVector.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"

namespace endfieldmodel {

    Tanh* Tanh::clone() const { return new Tanh(*this); }

    Tanh::~Tanh() {
        config_m.lambdaPower_m = Kokkos::View<double*>();  // deallocation
    }

    std::ostream& Tanh::print(std::ostream& out) const {
        out << "Tanh end field model with centre length: " << getX0()
            << " end length: " << getLambda();
        return out;
    }

    std::vector<std::vector<std::vector<int> > > Tanh::tdi_m;

    Tanh::Tanh(double x0, double lambda, int max_index) {
        setX0(x0);
        setLambda(lambda);
        setMaximumDerivative(max_index);
        initialise();
    }

    Tanh::Tanh(const Tanh& rhs) : config_m(rhs.config_m) { initialise(); }

    double Tanh::function(double x, int n) const {
        auto xDevice = Kokkos::View<double*>("efm::Tanh::xDevice", 1);
        auto xHost   = Kokkos::create_mirror_view(xDevice);
        xHost(0)     = x;
        Kokkos::deep_copy(xDevice, xHost);
        auto yDevice = Kokkos::View<double**>("efm::Tanh::yDevice", 1, n + 1);
        auto lHost   = Kokkos::create_mirror_view(config_m.lambdaPower_m);
        Kokkos::deep_copy(lHost, config_m.lambdaPower_m);

        function(xDevice, n, yDevice);
        auto yHost = Kokkos::create_mirror_view(yDevice);
        Kokkos::deep_copy(yHost, yDevice);

        double value = yHost(0, n);
        return value;
    }

    void Tanh::setMaximumDerivative(size_t n) {
        if (n > config_m.maxDerivative_m) {
            throw GeneralOpalException(
                    "Tanh::setMaximumDerivative",
                    "GPU-compatible tanh derivatives are limited to order 21");
        }
        setTanhDiffIndices(n);
        for (size_t i = 0; i < config_m.coefficientCount_m; ++i) {
            config_m.coefficients_m[i] = 0;
        }
        for (size_t derivative = 0; derivative <= n; ++derivative) {
            for (const auto& term : tdi_m[derivative]) {
                config_m.coefficients_m[derivative * (config_m.maxDerivative_m + 2) + term[1]] =
                        term[0];
            }
        }
    }

    void Tanh::setTanhDiffIndices(size_t n) {
        tdi_m.reserve(n + 1);
        if (tdi_m.size() == 0) {
            tdi_m.push_back(std::vector<std::vector<int> >(1, std::vector<int>(2)));
            tdi_m[0][0][0] = 1;  // 1*tanh(x) - third index is redundant
            tdi_m[0][0][1] = 1;
        }
        for (size_t i = tdi_m.size(); i < n + 1; ++i) {
            tdi_m.push_back(std::vector<std::vector<int> >());
            for (size_t j = 0; j < tdi_m[i - 1].size(); ++j) {
                int value = tdi_m[i - 1][j][1];
                if (value != 0) {
                    std::vector<int> new_vec(tdi_m[i - 1][j]);
                    new_vec[0] *= value;
                    new_vec[1] -= 1;
                    tdi_m[i].push_back(new_vec);
                    std::vector<int> new_vec2(tdi_m[i - 1][j]);
                    new_vec2[0] *= -value;
                    new_vec2[1] += 1;
                    tdi_m[i].push_back(new_vec2);
                }
            }
            tdi_m[i] = compactVector(tdi_m[i]);
        }
    }

    std::vector<std::vector<int> > Tanh::getTanhDiffIndices(size_t n) {
        setTanhDiffIndices(n);
        return tdi_m[n];
    }

    void Tanh::rescale(double scaleFactor) {
        config_m.x0_m *= scaleFactor;
        config_m.lambda_m *= scaleFactor;
        initialise();
    }

    void Tanh::initialise() {
        config_m.lambdaPower_m =
                Kokkos::View<double*>("efm::Tanh::lambdaPower", config_m.maxDerivative_m);
        auto lHost = Kokkos::create_mirror_view(config_m.lambdaPower_m);
        lHost[0]   = 1.0;
        for (size_t i = 1; i < config_m.maxDerivative_m; ++i) {
            lHost(i) = lHost(i - 1) * config_m.lambda_m;
        }
        Kokkos::deep_copy(config_m.lambdaPower_m, lHost);
    }

    void Tanh::finalise() { config_m.lambdaPower_m = Kokkos::View<double*>(); }

}  // namespace endfieldmodel
