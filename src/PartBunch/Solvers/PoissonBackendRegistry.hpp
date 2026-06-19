#ifndef OPAL_POISSON_BACKEND_REGISTRY_H
#define OPAL_POISSON_BACKEND_REGISTRY_H

#include <memory>
#include <string>
#include <vector>

#include "Manager/FieldSolverBase.h"
#include "PartBunch/Solvers/PoissonBackend.hpp"

/**
 * @brief Factory and metadata table for supported Poisson backends.
 *
 * This is the only place that maps OPAL solver names to backend wrappers. Adding a new IPPL
 * solver should require a new `PoissonBackend` implementation plus one registry entry here.
 */
template <typename T, unsigned Dim>
class PoissonBackendRegistry {
public:
    /// @brief Construct the backend wrapper for @p solverName.
    static std::unique_ptr<PoissonBackend<T, Dim>> create(
            const std::string& solverName, Solver_t<T, Dim>& solverVariant, Field_t<Dim>* rho,
            VField_t<T, Dim>* E, Field_t<Dim>* phi);

    /// @brief Return backend feature and policy flags without constructing the solver.
    static SolverCapabilities capabilitiesFor(const std::string& solverName);

    /// @brief Return the charge-density scaling constant for @p solverName.
    static T couplingConstantFor(const std::string& solverName);

    /// @brief Whether @p solverName has a registered backend.
    static bool supports(const std::string& solverName);

    /// @brief User-facing solver names currently registered.
    static std::vector<std::string> solverNames();
};

#endif  // OPAL_POISSON_BACKEND_REGISTRY_H
