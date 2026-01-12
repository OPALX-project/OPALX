#ifndef OPALX_SRC_PARTBUNCH_BCHANDLER_HPP
#define OPALX_SRC_PARTBUNCH_BCHANDLER_HPP

#include <array>
#include <ostream>

#include "Utilities/OpalException.h"

/**
 * @brief Handler for boundary conditions per spatial dimension.
 *
 * `BCHandler<Dim>` stores a boundary-condition type for each of the
 * `Dim` dimensions. It offers simple construction, string-to-enum parsing,
 * and some basic functionality used in the field container environment.
 *
 * Usage examples:
 * - Construct for 3D: `BCHandler<3> b(BCHandler<3>::OPEN, 
 *                                     BCHandler<3>::PERIODIC, 
 *                                     BCHandler<3>::OPEN);`
 * - More likely: construct from strToBCType, e.g. from `FieldSolverCmd`:
 *   `BCH_t::strToBCType(Attributes::getString(itsAttr[FIELDSOLVER::BCFFTX]))`
 *
 * @tparam Dim Number of spatial dimensions.
 */
template <unsigned Dim>
class BCHandler {
public:
    /**
     * @brief Supported boundary condition types.
     *
     * - `OPEN`: non-periodic/open boundary.
     * - `PERIODIC`: periodic boundary across the corresponding dimension.
     * 
     * @note This enum needs to be updated should there be other supported BC
     * types in the future (`DIRICHLET`!).
     */
    enum BCType { 
        OPEN, 
        PERIODIC 
    };

    /**
     * @brief Convert a textual boundary-condition name to `BCType` enum.
     *
     * Currently recognised strings are `"OPEN"` and `"PERIODIC"`. The strings
     * come from the pre-defined strings in the `FieldSolverCmd` class. If these
     * don't match, you might get an exception during current `FieldSolver`
     * construction.
     *
     * @param str Input string representing the boundary condition.
     * @return Corresponding `BCType` value.
     * @throws OpalException if the string does not match an implemented BC type.
     */
    static BCType strToBCType(const std::string& str) {
        if (str == "OPEN") {
            return OPEN;
        } else if (str == "PERIODIC") {
            return PERIODIC;
        } else {
            throw OpalException("BCHandler::strToBCType", 
                                "Unknown boundary condition type: " + str);
        }
    }

    /**
     * @brief Defaulted copy constructor.
     */
    BCHandler(const BCHandler& other) = default;

    /**
     * @brief Variadic constructor accepting exactly `Dim` BC values.
     *
     * Each argument should be convertible to `BCType`. The constructor 
     * validates that the number of provided BCs equals `Dim` and throws an
     * `OpalException` otherwise.
     *
     * @tparam Bcs Variadic parameter pack type.
     * @param bcs Boundary condition values for each dimension.
     * @throws OpalException when the number of provided BCs != `Dim`.
     */
    template <typename... Bcs>
    BCHandler(Bcs... bcs)
    : bcs_m{{static_cast<BCType>(bcs)...}} {
        if (sizeof...(bcs) != Dim) {
            throw OpalException("BCHandler::BCHandler", 
                                "Number of passed BCs does not match dimensionality!");
        }
    }

    /**
     * @brief Return `true` if every stored BC equals `bc_type`.
     *
     * @param bc_type Boundary condition value to compare against.
     * @return `true` if all dimensions have `bc_type`, otherwise `false`.
     */
    bool isAll(BCType bc_type) const {
        for (unsigned int d = 0; d < Dim; ++d) {
            if (bcs_m[d] != bc_type) return false;
        }
        return true;
    }

    /**
     * @brief Check whether all stored BCs are equal (all the same value).
     *
     * This function is necesary, since current IPPL does not support mixed BCs,
     * so this consistency check validates user input. 
     *
     * @return `true` when every element in the internal array equals the first 
     * element.
     */
    bool isAllEqual() const {
        bool base_case = bcs_m[0];
        for (unsigned int d = 1; d < Dim; ++d) {
            if (bcs_m[d] != base_case) return false;
        }
        return true;
    }

    /**
     * @brief Stream output helper for debugging and logging.
     *
     * Produces a compact, human-readable representation such as:
     * `* BCHandler<3>: [OPEN, PERIODIC, OPEN]`.
     */
    friend std::ostream& operator<<(std::ostream& os, const BCHandler& h) {
        os << "* BCHandler<" << Dim << ">: [";
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
    /// Internal storage of boundary conditions, one per dimension.
    std::array<BCType, Dim> bcs_m;
};

#endif // OPALX_SRC_PARTBUNCH_BCHANDLER_HPP