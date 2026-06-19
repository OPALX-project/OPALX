#ifndef OPAL_FIELD_DIAGNOSTICS_H
#define OPAL_FIELD_DIAGNOSTICS_H

#include <cstddef>
#include <string>

#include "Manager/BaseManager.h"
#include "Manager/FieldSolverBase.h"

// Owns field-dump diagnostics and solver call-count bookkeeping for FieldSolver.
// Dirichlet-plane diagnostics stay in BinnedFieldSolver for this stage because
// they depend on binning and correction state.
template <typename T, unsigned Dim>
class FieldDiagnostics {
public:
    FieldDiagnostics() = default;

    void incrementCallCounter() { ++callCounter_m; }
    void resetCallCounter() { callCounter_m = 0; }
    std::size_t getCallCounter() const { return callCounter_m; }

    void dumpScalarField(
            const std::string& what, const Field_t<Dim>* rho, const Field_t<Dim>* phi,
            const std::string& solverType) const;
    void dumpVectorField(const std::string& what, const VField_t<T, Dim>* E) const;

private:
    std::size_t callCounter_m = 0;
};

template <>
void FieldDiagnostics<double, 3>::dumpScalarField(
        const std::string& what, const Field_t<3>* rho, const Field_t<3>* phi,
        const std::string& solverType) const;

template <>
void FieldDiagnostics<double, 3>::dumpVectorField(
        const std::string& what, const VField_t<double, 3>* E) const;

extern template class FieldDiagnostics<double, 3>;

#endif  // OPAL_FIELD_DIAGNOSTICS_H
