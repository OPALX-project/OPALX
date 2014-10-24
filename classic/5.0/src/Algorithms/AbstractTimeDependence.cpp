#include "Algorithms/AbstractTimeDependence.h"
#include "Utilities/OpalException.h"

#include <sstream>

std::map<std::string, AbstractTimeDependence*> AbstractTimeDependence::td_map =
                              std::map<std::string, AbstractTimeDependence*>();

AbstractTimeDependence* AbstractTimeDependence::getTimeDependence(std::string name) {
    if (td_map.find(name) != td_map.end()) {
        return td_map[name];
    } else {
        throw OpalException("AbstractTimeDependence::getTimeDependence",
                            "Could not find TimeDependence called "+name);
    }
}

void AbstractTimeDependence::setTimeDependence(std::string name,
                                               AbstractTimeDependence* time_dep) {
    if (td_map.find(name) != td_map.end()) {
        delete td_map[name];
    }
    td_map[name] = time_dep;
}

std::string AbstractTimeDependence::getName(AbstractTimeDependence* time_dep) {
    typedef std::map<std::string, AbstractTimeDependence*>::iterator iter;
    for (iter i = td_map.begin(); i != td_map.end(); ++i) {
        if (i->second == time_dep)
            return i->first;
    }
    std::stringstream ss;
    ss << time_dep;
    throw OpalException("AbstractTimeDependence::getTimeDependence",
                        "Could not find TimeDependence with address "+ss.str());
}
