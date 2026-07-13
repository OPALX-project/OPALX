#ifndef OPALX_SRC_PARTBUNCH_BCHANDLER_HPP
#define OPALX_SRC_PARTBUNCH_BCHANDLER_HPP

#include <array>
#include <ostream>

#include "Ippl.h"
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
 *
 * @note This class is currently not designed to be used on device.
 */
template <unsigned Dim>
class BCHandler {
public:
    /**
     * @brief Supported boundary condition types.
     *
     * - `OPEN`: non-periodic/open boundary.
     * - `PERIODIC`: periodic boundary across the corresponding dimension.
     * - `DIRICHLET`: fixed potential (zero for now) at the boundary.
     *
     * @note This enum needs to be updated should there be other supported BC
     * types in the future (`DIRICHLET`!).
     */
    enum BCType { OPEN, PERIODIC, DIRICHLET };

    /**
     * @brief Convert a textual boundary-condition name to `BCType` enum.
     *
     * Currently recognised strings are `"OPEN"`, `"PERIODIC"`, and
     * `"DIRICHLET"`. The strings come from the pre-defined strings in the
     * `FieldSolverCmd` class. If these don't match, you might get an exception
     * during current `FieldSolver` construction.
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
        } else if (str == "DIRICHLET") {
            return DIRICHLET;
        } else {
            throw OpalException(
                    "BCHandler::strToBCType", "Unknown boundary condition type: " + str);
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
    BCHandler(Bcs... bcs) : bcs_m{{static_cast<BCType>(bcs)...}} {
        if (sizeof...(bcs) != Dim) {
            throw OpalException(
                    "BCHandler::BCHandler", "Number of passed BCs does not match dimensionality!");
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
        int base_case = bcs_m[0];
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
                case OPEN:
                    os << "OPEN";
                    break;
                case PERIODIC:
                    os << "PERIODIC";
                    break;
                case DIRICHLET:
                    os << "DIRICHLET";
                    break;
                default:
                    os << "UNKNOWN";
                    break;
            }
            if (d + 1 < Dim) os << ", ";
        }
        os << "]";
        return os;
    }

    /**
     * @brief Create an IPPL BConds container from the stored boundary conditions.
     *
     * This function constructs an `ippl::BConds<Field, Dim>` container and
     * populates it with appropriate face boundary condition objects based on the
     * boundary conditions stored in this handler. For each of the `2 * Dim` faces:
     * - If the corresponding dimension has `PERIODIC` BC, a `PeriodicFace` is created.
     * - If the corresponding dimension has `OPEN` BC, a `NoBcFace` is created.
     *
     * The resulting container can be directly used with IPPL field solvers via
     * the field's `setFieldBC()` method.
     *
     * @tparam Field The IPPL field type to which boundary conditions will be applied.
     * @return An `ippl::BConds<Field, Dim>` container populated with appropriate
     *         face boundary condition objects.
     *
     * @throws OpalException if an unsupported BC type is encountered.
     *
     * @example
     * ```cpp
     * BCHandler<3> handler(BCHandler<3>::PERIODIC,
     *                       BCHandler<3>::OPEN,
     *                       BCHandler<3>::PERIODIC);
     * auto bc_container = handler.toIPPLBConds<Field<double, 3>>();
     * field_ptr->setFieldBC(bc_container);
     * ```
     */
    template <typename Field>
    ippl::BConds<Field, Dim> toIPPLBConds() const {
        typedef ippl::BConds<Field, Dim> bc_type;
        bc_type bc_container;

        for (unsigned int face = 0; face < 2 * Dim; ++face) {
            unsigned int d = face / 2;  // Dimension associated with this face

            if (bcs_m[d] == PERIODIC) {
                bc_container[face] = std::make_shared<ippl::PeriodicFace<Field>>(face);
            } else if (bcs_m[d] == DIRICHLET) {
                bc_container[face] = std::make_shared<ippl::ZeroFace<Field>>(face);
            } else if (bcs_m[d] == OPEN) {
                bc_container[face] = std::make_shared<ippl::NoBcFace<Field>>(face);
            } else {
                throw OpalException("BCHandler::toIPPLBConds", "Unsupported BCType encountered.");
            }
        }

        return bc_container;
    }

private:
    /// Internal storage of boundary conditions, one per dimension.
    std::array<BCType, Dim> bcs_m;
};

#endif  // OPALX_SRC_PARTBUNCH_BCHANDLER_HPP