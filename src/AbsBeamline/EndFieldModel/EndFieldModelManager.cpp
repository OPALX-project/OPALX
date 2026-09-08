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

#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include <algorithm>
#include <map>
#include <sstream>
#include "Utilities/GeneralOpalException.h"

namespace endfieldmodel {
std::shared_ptr<EndFieldModelManager> EndFieldModelManager::globalEFM_m;

std::shared_ptr<EndFieldModel> EndFieldModelManager::getEndFieldModel(const std::string& name) {
    if (efmMap_m.find(name) == efmMap_m.end()) {
        throw OpalException("EndFieldModelManager::getEndFieldModel", "Could not find model '"+name+"'");
    }
    return efmMap_m[name];
}

std::shared_ptr<EndFieldModelManager> EndFieldModelManager::getEFMManager() {
        if (globalEFM_m) {
            return globalEFM_m;
        }
        globalEFM_m = std::make_shared<EndFieldModelManager>();
        return globalEFM_m;
}

void EndFieldModelManager::setEndFieldModel(const std::string& name,
                      const std::shared_ptr<EndFieldModel>& efm) {
    efmMap_m[name] = efm;
}


}  // namespace endfieldmodel
