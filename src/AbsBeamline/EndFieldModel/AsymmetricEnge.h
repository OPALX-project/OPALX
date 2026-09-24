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

#ifndef ENDFIELDMODEL_ASYMMETRICENGE_H_
#define ENDFIELDMODEL_ASYMMETRICENGE_H_

#include <iostream>
#include <memory>
#include <vector>

#include "AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "AbsBeamline/EndFieldModel/Enge.h"

namespace endfieldmodel {

/** Calculate the AsymmetricEnge function (e.g. for multipole end fields).
 *
 *  AsymmetricEnge function is given by\n
 *  \f$T(x) = (tanh( (x+x0)/\lambda )-tanh( (x-x0)/\lambda ))/2\f$\n
 *  The derivatives of tanh(x) are given by\n
 *  \f$d^p tanh(x)/dx^p = \sum_q I_{pq} tanh^{q}(x)\f$\n
 *  where \f$I_{pq}\f$ are calculated using some recursion relation. Using these
 *  expressions, one can calculate a recursion relation for higher order
 *  derivatives and hence calculate analytical derivatives at arbitrary order.
 */
struct AsymmetricEngeConfig {
    EngeConfig engeStart_m;
    EngeConfig engeEnd_m;
};

class AsymmetricEnge : public EndFieldModel {
public:
    /** Default constructor */
    AsymmetricEnge() = default;
    /** Constructor taking enge parameters */
    AsymmetricEnge(
            const std::vector<double> aStart, double x0Start, double lambdaStart,
            const std::vector<double> aEnd, double x0End, double lambdaEnd);

    /** Inheritable copy constructor. We take a deep copy of the engeStart
     *  and engeEnd
     */
    inline AsymmetricEnge* clone() const;

    /** Print a human-readable description of the end field model */
    std::ostream& print(std::ostream& out) const;

    /** Return the value of enge at some point x */
    inline double function(double x, int n) const;

    inline virtual void function(const Kokkos::View<double*>& xView,
                                 const int n,
                                 Kokkos::View<double**>& values) const override;

    /** Host side static wrapper */
    static inline void functionHost(
            const AsymmetricEngeConfig& config,
            const Kokkos::View<double*>& xView,
            const int n,
            Kokkos::View<double**>& values);

    static KOKKOS_INLINE_FUNCTION double functionDevice(
        const AsymmetricEngeConfig& config, double x, int n);


    /** Centre length is the average of x0End and x0Start */
    inline double getCentreLength() const;

    /** End length is the average of lambdaEnd and lambdaStart */
    inline double getEndLength() const;

    /** Return x0Start, offset of the start Enge */
    inline double getX0Start() const;

    /** Set x0Start, offset of the start Enge */
    inline void setX0Start(double x0);

    /** Return x0End, offset of the end Enge */
    inline double getX0End() const;

    /** Set x0End, offset of the end Enge */
    inline void setX0End(double x0);

    /** Return x0Start, offset of the start Enge */
    inline double getLambdaStart() const {return config_m.engeStart_m.lambda_m;}

    /** Return x0End, offset of the end Enge */
    inline double getLambdaEnd() const {return config_m.engeEnd_m.lambda_m;}

    /** Setup the Enge recursion for derivatives */
    inline void setMaximumDerivative(size_t n);

    /** Rescale the Enge to a new length scale */
    void rescale(double scaleFactor);

    AsymmetricEngeConfig getConfig() const {return config_m;}

private:
    AsymmetricEnge(const AsymmetricEnge& rhs) = default;
    AsymmetricEngeConfig config_m;
};

double AsymmetricEnge::function(double x, int n) const {
    Kokkos::View<double*> xView("tmpX", 1);
    Kokkos::deep_copy(xView, x);
    Kokkos::View<double**> valueView("tmpY", 1, n+1);
    function(xView, n, valueView);
    Kokkos::fence("Function calculation");
    double value(0);
    auto element = Kokkos::subview(valueView, 0, n);
    Kokkos::deep_copy(value, element);
    return value;
}

void AsymmetricEnge::function(const Kokkos::View<double*>& xView,
              const int n,
              Kokkos::View<double**>& values) const  {
    functionHost(config_m, xView, n, values);
}

void AsymmetricEnge::functionHost(
            const AsymmetricEngeConfig& config,
            const Kokkos::View<double*>& xView,
            const int n,
            Kokkos::View<double**>& values) {
    const size_t count = xView.extent(0);
    Kokkos::parallel_for(
        "AsymmetricEnge::functionHost()", count, KOKKOS_LAMBDA(const size_t i) {
            for (int j = 0; j < n+1; ++j)
                values(i, j) = functionDevice(config, xView(i), j);
    });
}

double AsymmetricEnge::functionDevice(const AsymmetricEngeConfig& config, double x, int n) {
    if (n == 0) {
        return (Enge::getEnge(config.engeStart_m, -x - config.engeStart_m.x0_m, n) +
                Enge::getEnge(config.engeEnd_m, x - config.engeEnd_m.x0_m, n))-1;
    } else {
        if (n % 2 == 1)
            return -Enge::getEnge(config.engeStart_m, -x - config.engeStart_m.x0_m, n) +
                    Enge::getEnge(config.engeEnd_m, x - config.engeEnd_m.x0_m, n);
        else
            return Enge::getEnge(config.engeStart_m, -x - config.engeStart_m.x0_m, n) +
                   Enge::getEnge(config.engeEnd_m, x - config.engeEnd_m.x0_m, n);
    }
}
double AsymmetricEnge::getX0Start() const { return config_m.engeStart_m.x0_m; }

double AsymmetricEnge::getX0End() const { return config_m.engeEnd_m.x0_m; }

void AsymmetricEnge::setX0Start(double x0) { config_m.engeStart_m.x0_m = x0; }

void AsymmetricEnge::setX0End(double x0) { config_m.engeEnd_m.x0_m = x0; }

AsymmetricEnge* AsymmetricEnge::clone() const { return new AsymmetricEnge(*this); }

void AsymmetricEnge::setMaximumDerivative(size_t n) {
    Enge::setEngeDiffIndices(n, config_m.engeStart_m);
    Enge::setEngeDiffIndices(n, config_m.engeEnd_m);
}

double AsymmetricEnge::getCentreLength() const {
    return config_m.engeStart_m.x0_m + config_m.engeEnd_m.x0_m;
}

double AsymmetricEnge::getEndLength() const {
    return config_m.engeStart_m.lambda_m + config_m.engeEnd_m.lambda_m;
}

}  // namespace endfieldmodel

#endif
