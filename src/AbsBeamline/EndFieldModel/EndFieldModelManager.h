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

#ifndef ENDFIELDMODEL_ENDFIELDMODELMANAGER_H_
#define ENDFIELDMODEL_ENDFIELDMODELMANAGER_H_

#include <iostream>
#include <map>
#include <memory>
#include <vector>

namespace endfieldmodel {

class EndFieldModel;

/** Singleton class to handle global register of EndFieldModels */
class EndFieldModelManager {
    public:
        EndFieldModelManager()  = default;
        ~EndFieldModelManager() = default;

        /** Return the global EndFieldModelManager.
         *
         *  If it is not initialised, initialise it.
         */
        static std::shared_ptr<EndFieldModelManager> getEFMManager();

        /** Clear the global EndFieldModelManager.
         */
        static void clearEFMManager() {globalEFM_m.reset();}

        /** Look up the EndFieldModel that has a given name
         *
         *  @param name: name of the EndFieldModel
         *
         *  @returns shared_ptr to the appropriate EndFieldModel.
         *
         *  @throws GeneralOpalException if name is not recognised
         */
        std::shared_ptr<EndFieldModel> getEndFieldModel(const std::string& name);

        /** Add a value to the lookup table
         *
         *  @param name: name of the EndFieldModel. If name already exists in the
         *  map, it is overwritten with the new value.
         *  @param efm: shared_ptr to the EndFieldModel.
         */
        void setEndFieldModel(const std::string& name,
                              const std::shared_ptr<EndFieldModel>& efm);

        /** Get the name corresponding to a given EndFieldModel
         *
         *  @param efm: EndFieldModel to lookup
         *
         *  @returns name corresponding to the EndFieldModel. Note that this
         *  just does a dumb loop over the stored map values; so O(N).
         *  @throws GeneralOpalException if efm is not recognised
         */
        template <class EFM>
        std::string getName(const std::shared_ptr<EndFieldModel>& efm);

    private:
        std::map<std::string, std::shared_ptr<EndFieldModel> > efmMap_m;
        static std::shared_ptr<EndFieldModelManager> globalEFM_m;

};

}  // namespace endfieldmodel

#endif
