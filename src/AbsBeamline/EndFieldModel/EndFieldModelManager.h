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

namespace endfieldmodel {

class Tanh;

class EndFieldModelManager {
    public:
        /** This is horrible. I really want to use an Abstraction but GPU does
         *  not allow it so I have to do if (type == BLAH) { do something } 
         */
        enum EndFieldModelType {kTANH, kENGE, kASYMMETRICENGE};

        EndFieldModelManager() {}
        ~EndFieldModelManager() {}

        /** Look up the EndFieldModel that has a given name
         *
         *  @param name: name of the EndFieldModel
         *
         *  @returns shared_ptr to the appropriate EndFieldModel.
         *  @throws GeneralOpalException if name is not recognised
         */
        template <class EFM>
        std::shared_ptr<EFM> getEndFieldModel(std::string name);

        /** Add a value to the lookup table
         *
         *  @param name: name of the EndFieldModel. If name already exists in the
         *  map, it is overwritten with the new value.
         *  @param efm: shared_ptr to the EndFieldModel.
         */
        template <class EFM>
        void setEndFieldModel(std::string name, std::shared_ptr<EFM> efm);

        /** Get the name corresponding to a given EndFieldModel
         *
         *  @param efm: EndFieldModel to lookup
         *
         *  @returns name corresponding to the EndFieldModel. Note that this
         *  just does a dumb loop over the stored map values; so O(N).
         *  @throws GeneralOpalException if efm is not recognised
         */
        template <class EFM>
        std::string getName(std::shared_ptr<EFM> efm);

    private:
        std::map<std::string, std::shared_ptr<Tanh> > tanhMap_m;
        std::map<std::string, std::shared_ptr<Enge> > engeMap_m;
        std::map<std::string, std::shared_ptr<AsymmetricEnge> > asymmtricEngeMap_m;
    };

    template <>
    std::shared_ptr<Tanh> getEndFieldModel(std::string name) {
        return tanhMap_m(name);
    }


    std::vector<std::vector<int> > CompactVector(std::vector<std::vector<int> > vec);

    /// CompactVector helper function, used for sorting
    bool GreaterThan(std::vector<int> v1, std::vector<int> v2);

    /** Return a == b if a and b are same size and a[i] == b[i] for all i.
     *
     *  The following operations must be defined for TEMP_ITER it:
     *    - ++it prefix increment operator
     *    - (*it) (that is unary *, i.e. dereference operator)
     *    - it1 != it2 not equals operator
     *    - (*it1) != (*it2) not equals operator of dereferenced object
     *
     *  Call like e.g. \n
     *      std::vector<int> a,b;\n
     *      bool test_equal = IterableEquality(a.begin(), a.end(), b.begin(),
     *                        b.end());\n
     *
     *  Can give a segmentation fault if a.begin() is not between a.begin() and
     *  a.end() (inclusive)
     */
    template <class TEMP_ITER>
    bool IterableEquality(TEMP_ITER a_begin, TEMP_ITER a_end, TEMP_ITER b_begin, TEMP_ITER b_end);

    /** Return a == b if a and b are same size and a[i] == b[i] for all i.
     *
     *  The following operations must be defined for TEMP_ITER it:
     *    - ++it prefix increment operator
     *    - (*it) (that is unary *, i.e. dereference operator)
     *    - it1 != it2 not equals operator
     *    - (*it1) != (*it2) not equals operator of dereferenced object
     *
     *  Call like e.g. \n
     *      std::vector<int> a,b;\n
     *      bool test_equal = IterableEquality(a.begin(), a.end(), b.begin(),
     *                        b.end());\n
     *
     *  Can give a segmentation fault if a.begin() is not between a.begin() and
     *  a.end() (inclusive)
     */
    template <class TEMP_ITER>
    bool IterableEquality(TEMP_ITER a_begin, TEMP_ITER a_end, TEMP_ITER b_begin, TEMP_ITER b_end);

    template <class TEMP_CLASS>
    bool IterableEquality(const TEMP_CLASS& a, const TEMP_CLASS& b) {
        return IterableEquality(a.begin(), a.end(), b.begin(), b.end());
    }

    template <class TEMP_ITER>
    bool IterableEquality(TEMP_ITER a_begin, TEMP_ITER a_end, TEMP_ITER b_begin, TEMP_ITER b_end) {
        TEMP_ITER a_it = a_begin;
        TEMP_ITER b_it = b_begin;
        while (a_it != a_end && b_it != b_end) {
            if (*a_it != *b_it) return false;
            ++a_it;
            ++b_it;
        }
        if (a_it != a_end || b_it != b_end) return false;
        return true;
    }

}  // namespace endfieldmodel

#endif
