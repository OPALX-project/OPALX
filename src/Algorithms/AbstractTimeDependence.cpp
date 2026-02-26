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

#include "Algorithms/AbstractTimeDependence.h"
#include <sstream>
#include <utility>
#include "Utilities/GeneralClassicException.h"

std::map<std::string, std::shared_ptr<AbstractTimeDependence> > AbstractTimeDependence::td_map =
    std::map<std::string, std::shared_ptr<AbstractTimeDependence> >();

std::shared_ptr<AbstractTimeDependence> AbstractTimeDependence::getTimeDependence(
    const std::string& name) {
    const auto pos = td_map.find(name);
    if (pos == td_map.end()) {
        throw GeneralClassicException(
            "AbstractTimeDependence::getTimeDependence",
            "Could not find TimeDependence called " + name);
    }
    return pos->second;
}

void AbstractTimeDependence::setTimeDependence(
    const std::string& name, std::shared_ptr<AbstractTimeDependence> time_dep) {
    td_map[name] = std::move(time_dep);
}

std::string AbstractTimeDependence::getName(
    const std::shared_ptr<AbstractTimeDependence>& time_dep) {
    for (auto& [name, dep] : td_map) {
        if (dep == time_dep) {
            return name;
        }
    }
    std::stringstream ss;
    ss << time_dep;
    throw GeneralClassicException(
        "AbstractTimeDependence::getTimeDependence",
        "Could not find TimeDependence with address " + ss.str());
}
