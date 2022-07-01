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
#include <memory>

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
class AsymmetricEnge : public EndFieldModel {
    public:
        /** Default constructor */
        AsymmetricEnge();
        /** Constructor taking enge parameters */
        AsymmetricEnge(const std::vector<double> aStart,
                       double x0Start,
                       double lambdaStart, 
                       const std::vector<double> aEnd,
                       double x0End,
                       double lambdaEnd);

        /** Inheritable copy constructor. We take a deep copy of the engeStart
         *  and engeEnd
         */
        inline AsymmetricEnge* clone() const;

        /** Print a human-readable description of the end field model */
        inline std::ostream& print(std::ostream& out) const;

        /** Return the value of enge at some point x */
        inline double function(double x, int n) const;

        /** Get the enge function for the magnet entrance */
        inline std::shared_ptr<Enge> getEngeStart() const;

        /** Set the enge function for the magnet entrance */
        inline void setEngeStart(std::shared_ptr<Enge> eStart);

        /** Get the enge function for the magnet exit */
        inline std::shared_ptr<Enge> getEngeEnd() const;

        /** Set the enge function for the magnet exit */
        inline void setEngeEnd(std::shared_ptr<Enge> eEnd);

        /** Return x0Start, offset of the start Enge */
        inline double getX0Start() const;

        /** Set x0Start, offset of the start Enge */
        inline void setX0Start(double x0);

        /** Return x0End, offset of the end Enge */
        inline double getX0End() const;

        /** Set x0End, offset of the end Enge */
        inline void setX0End(double x0);

        /** Setup the Enge recursion for derivatives */
        inline void setMaximumDerivative(size_t n);

        /** Rescale the Enge to a new length scale */
        void rescale(double scaleFactor);

    private:
        AsymmetricEnge(const AsymmetricEnge& rhs);
        std::shared_ptr<Enge> engeStart_m;
        std::shared_ptr<Enge> engeEnd_m;
        double x0Start_m;
        double x0End_m;
};

std::shared_ptr<Enge> AsymmetricEnge::getEngeStart() const {
    return engeStart_m;
}
std::shared_ptr<Enge> AsymmetricEnge::getEngeEnd() const {
    return engeEnd_m;
}
void AsymmetricEnge::setEngeStart(std::shared_ptr<Enge> enge) {
    engeStart_m = enge;
}
void AsymmetricEnge::setEngeEnd(std::shared_ptr<Enge> enge) {
    engeEnd_m = enge;
}

double AsymmetricEnge::getX0Start() const {
    return x0Start_m;
}

double AsymmetricEnge::getX0End() const {
    return x0End_m;
}

void AsymmetricEnge::setX0Start(double x0) {
    x0Start_m = x0;
}

void AsymmetricEnge::setX0End(double x0) {
    x0End_m = x0;
}

double AsymmetricEnge::function(double x, int n) const {
    return engeStart_m->getEnge(x-x0Start_m, n)+engeEnd_m->getEnge(-x-x0End_m, n);
}

std::ostream& AsymmetricEnge::print(std::ostream& out) const {
    out << "AsymmetricEnge start ";
    engeStart_m->print(out);
    out  << " end ";
    engeStart_m->print(out);
    return out;
}

AsymmetricEnge* AsymmetricEnge::clone() const {
    return new AsymmetricEnge(*this);
}

void AsymmetricEnge::setMaximumDerivative(size_t n) {
    Enge::setEngeDiffIndices(n);
}


}

#endif

