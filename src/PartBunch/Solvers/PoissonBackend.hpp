#ifndef OPAL_POISSON_BACKEND_H
#define OPAL_POISSON_BACKEND_H

#include <string>

#include "PartBunch/Solvers/SolveRequest.hpp"
#include "PartBunch/Solvers/SolverCapabilities.hpp"

template <typename T, unsigned Dim>
class PoissonBackend {
public:
    virtual ~PoissonBackend() = default;

    virtual std::string name() const = 0;
    virtual SolverCapabilities capabilities() const = 0;
    virtual void solve(const SolveRequest<T, Dim>& request) = 0;
    virtual void refreshAfterLayoutChange() = 0;
    virtual T couplingConstant() const = 0;
};

#endif  // OPAL_POISSON_BACKEND_H
