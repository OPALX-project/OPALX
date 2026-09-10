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

#ifndef ENDFIELDMODEL_ENDFIELDMODEL_H_
#define ENDFIELDMODEL_ENDFIELDMODEL_H_

#include <iostream>
#include <map>
#include <memory>
#include <vector>
#include "VectorMath.h"

namespace endfieldmodel {

    class EndFieldModel {
    public:
        /** Destructor */
        virtual ~EndFieldModel() { ; }

        /** Stream a human readable description of the end field model to out */
        virtual std::ostream& print(std::ostream& out) const = 0;

        /** Return the value of the function or its n^th derivative
         *
         *  @param x: returns d^n f(x)/dx^n
         *  @param n: the derivative
         */
        [[nodiscard]] virtual double function(double x, int n) const = 0;

        /** GPU-aware version of function
         *
         *  @param x: returns d^n f(x)/dx^n
         *  @param n: the derivative
         */
        virtual void function(Kokkos::View<double*> xView,  const int& n, Kokkos::View<double**> values) = 0;

        /** Return the nominal flat top length of the magnet
         */
        [[nodiscard]] virtual double getCentreLength() const = 0;

        /** Return the nominal end field length of the magnet
         */
        [[nodiscard]] virtual double getEndLength() const = 0;

        /** Inheritable copy constructor - returns a deep copy of the EndFieldModel */
        [[nodiscard]] virtual EndFieldModel* clone() const = 0;

        /** Set the maximum derivative that will be required to be calculated
         *
         *  Some end field models e.g. Enge use recursion relations to calculate
         *  analytically derivatives at high order. By setting the maximum derivative
         *  these models can set up the tables of recursion coefficients at set-up
         *  time which makes the derivative lookup faster.
         */
        virtual void setMaximumDerivative(size_t n) = 0;

        /** Rescale the end field lengths and offsets by a factor x0
         *
         *  If before rescaling the endfieldmodel returns f(x), after rescaling the
         *  endfieldmodel should return f(x*scaleFactor)
         */
        virtual void rescale(double scaleFactor) = 0;
    private:
    };

}  // namespace endfieldmodel

#endif
