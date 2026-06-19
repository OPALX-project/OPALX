#ifndef OPAL_POISSON_BACKEND_H
#define OPAL_POISSON_BACKEND_H

#include <string>

#include "PartBunch/Solvers/SolveRequest.hpp"
#include "PartBunch/Solvers/SolverCapabilities.hpp"

/**
 * @brief Solver-independent interface used by OPALX self-field code.
 *
 * Implementations own the concrete IPPL solver access for one backend. Callers pass normalized
 * solve requests and inspect capabilities instead of downcasting to IPPL solver types.
 */
template <typename T, unsigned Dim>
class PoissonBackend {
public:
    virtual ~PoissonBackend() = default;

    /// @brief User-facing solver name handled by this backend.
    virtual std::string name() const = 0;

    /// @brief Feature and policy flags required by orchestration code.
    virtual SolverCapabilities capabilities() const = 0;

    /// @brief Execute one Poisson solve for the current field storage.
    virtual void solve(const SolveRequest<T, Dim>& request) = 0;

    /// @brief Rebind backend field references after a mesh or layout change.
    virtual void refreshAfterLayoutChange() = 0;

    /// @brief Scaling constant applied when preparing charge density for this backend.
    virtual T couplingConstant() const = 0;
};

#endif  // OPAL_POISSON_BACKEND_H
