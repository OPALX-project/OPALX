#ifndef OPALX_Port_HH
#define OPALX_Port_HH

#include "Algorithms/CoordinateSystemTrafo.h"

#include <string>

/**
 * @class Port
 * @brief Named local frame attached to an element body.
 *
 * A port stores the rigid transform from the canonical body frame of an element
 * to a named interface frame such as `entry`, `body`, or `exit`.
 */
class Port {
public:
    Port() = default;
    Port(const std::string& name, const CoordinateSystemTrafo& bodyToPort)
        : name_m(name), bodyToPort_m(bodyToPort) {}

    const std::string& getName() const { return name_m; }
    const CoordinateSystemTrafo& getBodyToPort() const { return bodyToPort_m; }

private:
    std::string name_m;
    CoordinateSystemTrafo bodyToPort_m;
};

#endif
