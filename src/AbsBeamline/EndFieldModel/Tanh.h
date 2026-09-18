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

#ifndef ENDFIELDMODEL_TANH_H_
#define ENDFIELDMODEL_TANH_H_

#include <iostream>
#include <vector>

#include "AbsBeamline/EndFieldModel/EndFieldModel.h"

namespace endfieldmodel {

/** Calculate the Tanh function (e.g. for multipole end fields).
 *
 *  DoubleTanh function is given by\n
 *  \f$T(x) = (tanh( (x+x0)/\lambda )-tanh( (x-x0)/\lambda ))/2\f$\n
 *  The derivatives of tanh(x) are given by\n
 *  \f$d^p tanh(x)/dx^p = \sum_q I_{pq} tanh^{q}(x)\f$\n
 *  where \f$I_{pq}\f$ are calculated using some recursion relation. Using these
 *  expressions, one can calculate a recursion relation for higher order
 *  derivatives and hence calculate analytical derivatives at arbitrary order.
 */
struct TanhConfig {
    double x0_m = 0.;
    double lambda_m = 0.;
    static constexpr size_t maxDerivative_m = 21;
    static constexpr size_t coefficientCount_m =
            (maxDerivative_m + 1) * (maxDerivative_m + 2);
    Kokkos::Array<int, coefficientCount_m> coefficients_m{};
};

class Tanh : public EndFieldModel {
public:
    Tanh(double x0, double lambda, int max_index);

    /** Default constructor (initialises x0 and lambda to 0) */
    Tanh() = default;

    /** Copy constructor */
    Tanh(const Tanh& rhs) = default;

    /** Destructor (no mallocs so does nothing) */
    ~Tanh() = default;
    /** Inherited copy constructor. */
    Tanh* clone() const override;

    endfieldmodel::Tanh& operator=(const endfieldmodel::Tanh& rhs) = default;

    double function(double x, int n) const;

    /** Device-callable evaluation using coefficients prepared on the host. */
    static KOKKOS_INLINE_FUNCTION double functionDevice(const TanhConfig& config, double x, int n);

    void function(const Kokkos::View<double*>& xView,
                  const int& maxDerivative,
                  Kokkos::View<double**>& derivatives);


    static void functionHost(const TanhConfig& config,
                        const Kokkos::View<double*>& xView,
                        const int& maxDerivative,
                        Kokkos::View<double**>& derivatives);

    /** Returns the value of tanh((x+x0)/lambda) or its \f$n^{th}\f$ derivative. */
    double getTanh(double x, int n) const;

    /** Returns the value of tanh((x-x0)/lambda) or its \f$n^{th}\f$ derivative. */
    double getNegTanh(double x, int n) const;

    /** Get all the tanh differential indices \f$I_{pq}\f$.
     *
     *  Returns vector of vector of ints where p indexes the differential and
     * q indexes the tanh power - so
     */
    static std::vector<std::vector<int> > getTanhDiffIndices(size_t n);

    /** Nominal flat top length is twice x0 (one x0 in each direction) */
    double getCentreLength() const override { return config_m.x0_m * 2.0; }

    /** Return nominal fringe field length */
    double getEndLength() const override { return config_m.lambda_m; }

    /** Set the value of tanh differential indices to nth order differentials. */
    static void setTanhDiffIndices(size_t n);

    /** Return lambda (end length) */
    inline double getLambda() const { return config_m.lambda_m; }

    /** Return x0 (half the flat top length) */
    inline double getX0() const { return config_m.x0_m; }

    /** Set lambda (end length) */
    inline void setLambda(double lambda) { config_m.lambda_m = lambda; }

    /** Set x0 (flat top length) */
    inline void setX0(double x0) { config_m.x0_m = x0; }

    /** Set the maximum derivative prior to tracking */
    void setMaximumDerivative(size_t n) override;

    /** Rescale the endfield */
    void rescale(double scalefactor);

    /** Create a double tanh function
     *
     *  Here x0 is the centre length and lambda is the end length. max_index is
     *  used to set up for differentiation - don't try to calculate
     *  higher differentials than exist in max_index.
     *
     *  This is a thin wrapper for the Tanh providing interface
     *  to EndFieldModel
     */

    /** GPU aware version of the function
     *
     */
    void function(const Kokkos::View<double*>& xView,
                  const int maxDerivative,
                  Kokkos::View<double**>& derivatives) const override;

    /** Print summary of the Tanh model to out */
    std::ostream& print(std::ostream& out) const override;

    /** Return the trivially-copyable data used inside device kernels. */
    TanhConfig getConfig() const { return config_m; }

private:
    TanhConfig config_m;

    /** _tdi indexes powers of tanh in d^n tanh/dx^n as sum of powers of tanh
     *
     *  For some reason we index as n, +a, -a, but the third index is redundant
     */
    static std::vector<std::vector<std::vector<int> > > tdi_m;
};

inline void Tanh::function(const Kokkos::View<double*>& xView,  const int maxDerivative, Kokkos::View<double**>& derivatives) const {
    Tanh::functionHost(config_m, xView, maxDerivative, derivatives);
}

KOKKOS_INLINE_FUNCTION
double Tanh::functionDevice(const TanhConfig& config, double x, int n) {
    const double tanhPositive = Kokkos::tanh((x + config.x0_m) / config.lambda_m);
    const double tanhNegative = Kokkos::tanh((x - config.x0_m) / config.lambda_m);
    double positivePower = 1.;
    double negativePower = 1.;
    double result = 0.;
    double lambdaPower = 1.;
    for (int i = 0; i < n; ++i) {
        lambdaPower *= config.lambda_m;
    }
    for (int power = 0; power <= n + 1; ++power) {
        const int coefficient = config.coefficients_m[n * (config.maxDerivative_m + 2) + power];
        result += static_cast<double>(coefficient) * (positivePower - negativePower);
        positivePower *= tanhPositive;
        negativePower *= tanhNegative;
    }
    return result / (2. * lambdaPower);
}

inline void Tanh::function(const Kokkos::View<double*>& xView,
              const int& maxDerivative,
              Kokkos::View<double**>& derivatives) {
    functionHost(config_m, xView, maxDerivative, derivatives);
}

inline void Tanh::functionHost(
                        const TanhConfig& config,
                        const Kokkos::View<double*>& xView,
                        const int& maxDerivative,
                        Kokkos::View<double**>& derivatives) {
    const size_t count = xView.size();
    Kokkos::parallel_for(
        "Tanh::function", count, KOKKOS_LAMBDA(const size_t i) {
            for (int order = 0; order < maxDerivative; ++order) {
                double x = xView(i);
                derivatives(i, order) = Tanh::functionDevice(config, x, order);
            }
        }
    );
}

}  // namespace endfieldmodel


#endif
