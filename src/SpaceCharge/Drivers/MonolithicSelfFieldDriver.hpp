#ifndef OPALX_MONOLITHIC_SELF_FIELD_DRIVER_HPP
#define OPALX_MONOLITHIC_SELF_FIELD_DRIVER_HPP

#include "SpaceCharge/Drivers/SelfFieldDriver.hpp"

/**
 * @brief Legacy single-solve self-field driver.
 *
 * This driver keeps the non-binned path as one monolithic scatter, solve, and gather pass. It is
 * selected when the bunch has no adaptive binning object.
 */
template <typename T, unsigned Dim>
class MonolithicSelfFieldDriver final : public SelfFieldDriver<T, Dim> {
public:
    using Solver_t    = typename SelfFieldDriver<T, Dim>::Solver_t;
    using PartBunch_t = typename SelfFieldDriver<T, Dim>::PartBunch_t;

    const char* name() const override { return "monolithic"; }
    void compute(
            Solver_t& solver, PartBunch_t& bunch,
            const BoundaryCorrection<T, Dim>* activeCorrection) override;
};

#endif  // OPALX_MONOLITHIC_SELF_FIELD_DRIVER_HPP
