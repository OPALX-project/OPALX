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
class Tanh : public EndFieldModel {
  public:
    /** Create a double tanh function
     *
     *  Here x0 is the centre length and lambda is the end length. max_index is
     *  used to set up for differentiation - don't try to calculate
     *  higher differentials than exist in max_index.
     */
    Tanh(double x0, double lambda, int max_index);

    /** Default constructor (initialises x0 and lambda to 0) */
    Tanh() : _x0(0.), _lambda(0.) {SetTanhDiffIndices(10);}

    /** Destructor (no mallocs so does nothing) */
    ~Tanh() {}

    /** Returns the value of DoubleTanh or its \f$n^{th}\f$ derivative. */
    Tanh* Clone() const;

    /** Double Tanh is given by\n
     *  \f$d(x) = \f$
     */
    double GetDoubleTanh(double x, int n) const;

    /** Returns the value of tanh((x+x0)/lambda) or its \f$n^{th}\f$ derivative. */
    double GetTanh(double x, int n) const;

    /** Returns the value of tanh((x-x0)/lambda) or its \f$n^{th}\f$ derivative. */
    double GetNegTanh(double x, int n) const;

    /** Get all the tanh differential indices \f$I_{pq}\f$.
     *
     *  Returns vector of vector of ints where p indexes the differential and
     * q indexes the tanh power - so
     */
    static std::vector< std::vector<int> > GetTanhDiffIndices(size_t n);

    /** Set the value of tanh differential indices to nth order differentials. */
    static void SetTanhDiffIndices(size_t n);

    /** Return lambda (end length) */
    inline double GetLambda() const {return _lambda;}

    /** Return x0 (flat top length) */
    inline double GetX0() const {return _x0;}

    /** Set lambda (end length) */
    inline void   SetLambda(double lambda) {_lambda = lambda;}

    /** Set x0 (flat top length) */
    inline void   SetX0(double x0)     {_x0 = x0;}
  private:
    double _x0, _lambda;

    /** _tdi indexes powers of tanh in d^n tanh/dx^n as sum of powers of tanh
     * 
     *  For some reason we index as n, +a, -a, but the third index is redundant
     */
    static std::vector< std::vector< std::vector<int> > > _tdi;
};

}

#endif

