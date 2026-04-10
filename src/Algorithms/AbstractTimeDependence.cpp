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
