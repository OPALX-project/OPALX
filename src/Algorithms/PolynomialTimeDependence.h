//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
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

#ifndef ALGORITHMS_POLYNOMIALTIMEDEPENDENCE_H_
#define ALGORITHMS_POLYNOMIALTIMEDEPENDENCE_H_

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
    explicit PolynomialTimeDependence(const std::vector<double>& ptd) : coeffs(ptd) {}

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

inline Inform& operator<<(Inform& os, const PolynomialTimeDependence& p) { return p.print(os); }

#endif
