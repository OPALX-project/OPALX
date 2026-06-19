#ifndef OPAL_FIELD_DIAGNOSTICS_H
#define OPAL_FIELD_DIAGNOSTICS_H

#include <cstddef>
#include <string>

#include "Manager/BaseManager.h"
#include "Manager/FieldSolverBase.h"
#include "SpaceCharge/Solvers/SolverCapabilities.hpp"

/**
 * @brief Owns field-dump diagnostics and solver call-count bookkeeping for `FieldSolver`.
 *
 * Dirichlet-plane diagnostics stay in `BinnedFieldSolver` for this stage because they depend on
 * binning and correction state.
 */
template <typename T, unsigned Dim>
class FieldDiagnostics {
public:
    FieldDiagnostics() = default;

    /// @brief Increment the number of backend solve calls recorded for the current run.
    void incrementCallCounter() { ++callCounter_m; }

    /// @brief Reset backend solve-call accounting, usually after warm-up.
    void resetCallCounter() { callCounter_m = 0; }

    /// @brief Return the number of counted backend solve calls.
    std::size_t getCallCounter() const { return callCounter_m; }

    /**
     * @brief Dump a scalar mesh field for OPALX field-debug diagnostics.
     *
     * `capabilities` decides whether `PHI` is read from a separate potential field or from the
     * in-place scalar field used by the backend.
     */
    void dumpScalarField(
            const std::string& what, const Field_t<Dim>* rho, const Field_t<Dim>* phi,
            const SolverCapabilities& capabilities) const;

    /**
     * @brief Dump a vector mesh field for OPALX field-debug diagnostics.
     */
    void dumpVectorField(const std::string& what, const VField_t<T, Dim>* E) const;

private:
    std::size_t callCounter_m = 0;
};

template <>
void FieldDiagnostics<double, 3>::dumpScalarField(
        const std::string& what, const Field_t<3>* rho, const Field_t<3>* phi,
        const SolverCapabilities& capabilities) const;

template <>
void FieldDiagnostics<double, 3>::dumpVectorField(
        const std::string& what, const VField_t<double, 3>* E) const;

extern template class FieldDiagnostics<double, 3>;

#endif  // OPAL_FIELD_DIAGNOSTICS_H
