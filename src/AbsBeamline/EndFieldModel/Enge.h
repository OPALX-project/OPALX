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

#ifndef ENDFIELDMODEL_ENGE_H_
#define ENDFIELDMODEL_ENGE_H_

#include <iostream>
#include <vector>

#include "AbsBeamline/EndFieldModel/CompactVector.h"
#include "AbsBeamline/EndFieldModel/EndFieldModel.h"

namespace endfieldmodel {

/** Enge class is a symmetric Enge function for analytical field models
 *
 *  Enge function is
 *
 *  \f$f(x) = 1/(1+exp(h(x-x0)))+1/(1+exp(h(-x-x0)))-1\f$.
 *
 *  where h is a polynomial in x/lambda with polynomial coefficients a. Also
 *  use g(x) = 1+exp(h)
 */


struct EngeConfig {
    std::vector<double> a_m;
    double lambda_m = 0.0;
    double x0_m = 0.0;

    constexpr int max_derivative_m = 12;
    /** Indexes the derivatives of enge in terms of g */
    static std::vector<std::vector<std::vector<int> > > q_m;
    /** Indexes the derivatives of g in terms of h */
    static std::vector<std::vector<std::vector<int> > > h_m;
};

class Enge : public EndFieldModel {
public:
    /** Default constructor */
    Enge() { setEngeDiffIndices(10); }
    /** Builds Enge function with parameters a_0, a_1, ..., lambda and x0.
     *
     *  Note that this class is in the inner loop of tracking, so many function
     *  calls are _not_ checked for correct indexing. Call setMaximumDerivative
     *  before use.
     */
    Enge(std::vector<double> a, double x0, double lambda);

    /** Destructor - no mallocs, so does nothing */
    ~Enge() = default;

    /** Inheritable copy constructor - no mallocs, so does nothing */
    [[nodiscard]] Enge* clone() const;

    /** Rescale so Enge(x) -> Enge(scaleFactor*x)
     *
     *  Sets x0 to scaleFactor*x0 and lambda to scaleFactor*lambda
     */
    void rescale(double scaleFactor);

    /** Return the value of enge(x+x0) + enge(-x-x0) at some point x */
    [[nodiscard]] inline double function(double x, int n) const;

    /** GPU-aware version of function
     *
     *  @param xView: list of x values
     *  @param n: order of the derivative
     *  @param values: fill values with the derivatives. Should have same
     *         number of rows as xView and number of columns should be n
     */
    inline void function(const Kokkos::View<double*>& xView,
                  const int n,
                  Kokkos::View<double**>& values) const override;

    /** Host side static wrapper */
    static KOKKOS_INLINE_FUNCTION void functionHost(
            const EngeConfig& config,
            const Kokkos::View<double*>& xView,
            const int n,
            Kokkos::View<double**>& values);

    /** Potentially device-side function call
     *
     *  @param x: returns d^n f(x)/dx^n
     *  @param n: the derivative
     */
    static KOKKOS_INLINE_FUNCTION double functionDevice(const EngeConfig& config, double x, int n);

    /** Nominal end length is lambda */
    [[nodiscard]] inline double getEndLength() const;

    /** Nominal centre length is x0/2 */
    [[nodiscard]] inline double getCentreLength() const;

    /** Print human-readable version of enge */
    std::ostream& print(std::ostream& out) const override;

    /** Returns the enge polynomial coefficients (a_i) */
    [[nodiscard]] std::vector<double> getCoefficients() const { return config_m.a_m; }

    /** Sets the enge polynomial coefficients (a_i) */
    void setCoefficients(std::vector<double> a) { config_m.a_m = a; }

    /** Returns the value of lambda */
    [[nodiscard]] double getLambda() const { return config_m.lambda_m; }

    /** Sets the value of lambda */
    inline void setLambda(double lambda) { config_m.lambda_m = lambda; }

    /** Returns the value of x0 */
    [[nodiscard]] double getX0() const { return config_m.x0_m; }

    /** Sets the value of x0 */
    inline void setX0(double x0) { config_m.x0_m = x0; }

    /** Calls setEngeDiffIndices to set the maximum derivative */
    inline void setMaximumDerivative(size_t n);

    /** Get a copy of the Config object */
    inline EngeConfig getConfig() const {return config_m;}

    /** Get a copy of the Config object */
    inline void setConfig(EngeConfig& config) {config_m = config;}

    /** Returns the value of the Enge function or its \f$n^{th}\f$ derivative.
     *
     *  Please call setEngeDiffIndices(n) before calling if n > max_index
     */
    static double getEnge(const EngeConfig& config, double x, int n);

    /** Returns \f$Enge(x-x0) + Enge(-x-x0)-1\f$ and its derivatives */
    static inline double  getDoubleEnge(const EngeConfig& config, double x, int n);

    /** Returns \f$h(x)\f$ or its \f$n^{th}\f$ derivative.
     *
     *  Here \f$h(x) = a_0 + a_1 x/\lambda + a_2 x^2/lambda^2 + \ldots \f$
     *  Please call setEngeDiffIndices(n) before calling if n > max_index
     */
    static double hN(const EngeConfig& config, double x, int n);

    /** Returns \f$g(x)\f$ or its \f$n^{th}\f$ derivative.
     *
     *  Here \f$g(x) = 1+exp(h(x))\f$.
     *  Please call setEngeDiffIndices(n) before calling if n > max_index
     */
    static double gN(const EngeConfig& config, double x, int n);

    /** Recursively calculate the indices for Enge and H
     *
     *  This will calculate the indices for Enge and H that are required to
     *  calculate the differential up to order n.
     */
    static void setEngeDiffIndices(size_t n);

    /** Return the indices for calculating the nth derivative of Enge ito g(x) */
    inline static std::vector<std::vector<int> > getQIndex(int n);

    /** Return the indices for calculating the nth derivative of g(x) ito h(x) */
    inline static std::vector<std::vector<int> > getHIndex(int n);

private:
    Enge(const Enge& enge);
    Enge& operator=(const Enge& enge);
    EngeConfig config_m;

};

void Enge::setMaximumDerivative(size_t n) { Enge::setEngeDiffIndices(n); }

double Enge::function(double x, int n) const { return getDoubleEnge(config_m, x, n); }

std::vector<std::vector<int> > Enge::getQIndex(int n) { return EngeConfig::q_m[n]; }

std::vector<std::vector<int> > Enge::getHIndex(int n) { return EngeConfig::h_m[n]; }

double Enge::getDoubleEnge(const EngeConfig& config, double x, int n) {
    if (n == 0) {
        return -1+(getEnge(config, x - config.x0_m, n) + getEnge(config, -x - config.x0_m, n));
    } else {
        if (n % 2 == 1)
            return + getEnge(config, x - config.x0_m, n) - getEnge(config, -x - config.x0_m, n);
        else
            return + getEnge(config, x - config.x0_m, n) + getEnge(config, -x - config.x0_m, n);
    }
}

double Enge::getCentreLength() const { return config_m.x0_m * 2.0; }

double Enge::getEndLength() const { return config_m.lambda_m; }

void Enge::function(const Kokkos::View<double*>& xView,
              const int n,
              Kokkos::View<double**>& values) const {
    EngeConfig config = config_m;
    Enge::functionHost(config, xView, n, values);
}


void Enge::functionHost(
            const EngeConfig& config,
            const Kokkos::View<double*>& xView,
            const int n,
            Kokkos::View<double**>& values) {
    const size_t count = xView.size();
    Kokkos::parallel_for(
        "Enge::functionHost()", count, KOKKOS_LAMBDA(const size_t i) {
            for (int j = 0; j < n; ++j)
                values(i, j) = functionDevice(config, xView(i), i);
        });
}


double Enge::functionDevice(const EngeConfig& config, double x, int n) {
    return getDoubleEnge(config, x, n);
}


}  // namespace endfieldmodel

#endif
