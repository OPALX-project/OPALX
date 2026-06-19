#ifndef OPALX_SELF_FIELD_DRIVER_HPP
#define OPALX_SELF_FIELD_DRIVER_HPP

#include "SpaceCharge/Corrections/BoundaryCorrection.hpp"
#include "PartBunch/PartBunchFwd.h"

template <typename T, unsigned Dim>
class BinnedFieldSolver;

/**
 * @brief Strategy interface for self-field algorithm drivers.
 *
 * Drivers own the high-level compute path. The compatibility solver facade owns field storage,
 * backend access, and shared helpers, while the selected driver decides the monolithic or binned
 * algorithm body for the current bunch state.
 */
template <typename T, unsigned Dim>
class SelfFieldDriver {
public:
    using Solver_t    = BinnedFieldSolver<T, Dim>;
    using PartBunch_t = PartBunch<T, Dim>;

    virtual ~SelfFieldDriver() = default;

    virtual const char* name() const = 0;
    virtual void compute(
            Solver_t& solver, PartBunch_t& bunch,
            const BoundaryCorrection<T, Dim>* activeCorrection) = 0;
};

#endif  // OPALX_SELF_FIELD_DRIVER_HPP
