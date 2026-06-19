#ifndef OPAL_SOLVE_REQUEST_H
#define OPAL_SOLVE_REQUEST_H

#include <optional>

#include "Manager/BaseManager.h"

/**
 * @brief Normalized options for one backend solve call.
 *
 * OPALX correction and driver code should add backend features here instead of calling concrete
 * IPPL solver methods directly.
 */
template <typename T, unsigned Dim>
struct SolveRequest {
    /// @brief Optional Green's-function shift used by shifted-Greens Dirichlet correction.
    std::optional<ippl::Vector<T, Dim>> greensShift;

    /// @brief Whether this request needs shifted Green's-function support.
    bool hasShiftedGreens() const { return greensShift.has_value(); }
};

#endif  // OPAL_SOLVE_REQUEST_H
