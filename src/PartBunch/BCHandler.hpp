#ifndef OPALX_SRC_PARTBUNCH_BCHANDLER_HPP
#define OPALX_SRC_PARTBUNCH_BCHANDLER_HPP

#include <array>
#include <ostream>

#include "OpalException.h"

template <unsigned Dim>
class BCHandler {
public:
    enum BCType { 
        OPEN, 
        PERIODIC 
    };

    // Copy constructor
    BCHandler(const BCHandler& other) = default;

    template <typename... Bcs>
    BCHandler(Bcs... bcs)
    : bcs_m{{static_cast<BCType>(bcs)...}} {
        if (sizeof...(bcs) != Dim) {
            throw OpalException("BCHandler::BCHandler", 
                                "Number of passed BCs does not match dimensionality!");
        }
    }

    bool isAll(BCType bc_type) const {
        for (unsigned int d = 0; d < Dim; ++d) {
            if (bcs_m[d] != bc_type) return false;
        }
        return true;
    }

    bool isAllEqual() const {
        bool base_case = bcs_m[0];
        for (unsigned int d = 1; d < Dim; ++d) {
            if (bcs_m[d] != base_case) return false;
        }
        return true;
    }

    friend std::ostream& operator<<(std::ostream& os, const BCHandler& h) {
        os << "BCHandler<" << Dim << ">: [";
        for (unsigned int d = 0; d < Dim; ++d) {
            switch (h.bcs_m[d]) {
                case OPEN: os << "OPEN"; break;
                case PERIODIC: os << "PERIODIC"; break;
                default: os << "UNKNOWN"; break;
            }
            if (d + 1 < Dim) os << ", ";
        }
        os << "]";
        return os;
    }


private:
    std::array<BCType, Dim> bcs_m;
};

#endif // OPALX_SRC_PARTBUNCH_BCHANDLER_HPP