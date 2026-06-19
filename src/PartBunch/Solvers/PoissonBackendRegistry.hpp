#ifndef OPAL_POISSON_BACKEND_REGISTRY_H
#define OPAL_POISSON_BACKEND_REGISTRY_H

#include <memory>
#include <string>
#include <vector>

#include "Manager/FieldSolverBase.h"
#include "PartBunch/Solvers/PoissonBackend.hpp"

template <typename T, unsigned Dim>
class PoissonBackendRegistry {
public:
    static std::unique_ptr<PoissonBackend<T, Dim>> create(
            const std::string& solverName, Solver_t<T, Dim>& solverVariant, Field_t<Dim>* rho,
            VField_t<T, Dim>* E, Field_t<Dim>* phi);

    static SolverCapabilities capabilitiesFor(const std::string& solverName);
    static T couplingConstantFor(const std::string& solverName);
    static bool supports(const std::string& solverName);
    static std::vector<std::string> solverNames();
};

#endif  // OPAL_POISSON_BACKEND_REGISTRY_H
