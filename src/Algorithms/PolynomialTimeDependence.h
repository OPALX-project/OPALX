/*
 *  Copyright (c) 2014, Chris Rogers
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

#ifndef CLASSIC_SRC_ALGORITHMS_POLYNOMIALTIMEDEPENDENCE_H_
#define CLASSIC_SRC_ALGORITHMS_POLYNOMIALTIMEDEPENDENCE_H_

#include <vector>
#include "Algorithms/AbstractTimeDependence.h"

class Inform;

/** @class PolynomialTimeDependence
 *
 *  Time dependence that follows a polynomial, like
 *      p_0 + p_1*t + p_2*t^2 + ... + p_i*t^i + ...
 */
class PolynomialTimeDependence : public AbstractTimeDependence {
public:
    /** Constructor
     *
     *  @param ptd the polynomial coefficients p_i; can be of arbitrary length
     *  (user is responsible for issues like floating point precision).
     */
    explicit PolynomialTimeDependence(const std::vector<double>& ptd) : coeffs(ptd) {
    }

    /** Default Constructor makes a 0 length polynomial */
    PolynomialTimeDependence() = default;

    /** Destructor does nothing */
    ~PolynomialTimeDependence() override = default;

    /** Return the polynomial Sum_i p_i t^i; returns 0 if p is of 0 length */
    double getValue(double time) override;

    /** Return the integral of the polynomial from 0 to time */
    double getIntegral(double time) override;

    /** Inheritable copy constructor
     *
     *  @returns new PolynomialTimeDependence that is a copy of this. User owns
     *  returned memory.
     */
    PolynomialTimeDependence* clone() override;

    /** Print the polynomials
     *
     *  @param os "Inform" stream to which the polynomials are printed.
     */
    Inform& print(Inform& os) const;

    /** Return the polynomial coefficients */
    const std::vector<double>& getCoefficients() const { return coeffs; }

private:
    std::vector<double> coeffs;
};

inline Inform& operator<<(Inform& os, const PolynomialTimeDependence& p) {
    return p.print(os);
}

#endif
