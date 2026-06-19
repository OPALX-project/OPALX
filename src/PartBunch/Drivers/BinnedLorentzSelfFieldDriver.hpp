#ifndef OPALX_BINNED_LORENTZ_SELF_FIELD_DRIVER_HPP
#define OPALX_BINNED_LORENTZ_SELF_FIELD_DRIVER_HPP

#include "PartBunch/Drivers/SelfFieldDriver.hpp"

/**
 * @brief Per-bin quasi-static Lorentz self-field driver.
 *
 * This driver owns the binned algorithm orchestration: rebin, build per-bin context, run primary
 * and optional boundary-correction solves, accumulate lab-frame fields, then gather to particles.
 */
template <typename T, unsigned Dim>
class BinnedLorentzSelfFieldDriver final : public SelfFieldDriver<T, Dim> {
public:
    using Solver_t    = typename SelfFieldDriver<T, Dim>::Solver_t;
    using PartBunch_t = typename SelfFieldDriver<T, Dim>::PartBunch_t;

    const char* name() const override { return "binned-lorentz"; }
    void compute(
            Solver_t& solver, PartBunch_t& bunch,
            const BoundaryCorrection<T, Dim>* activeCorrection) override;
};

#endif  // OPALX_BINNED_LORENTZ_SELF_FIELD_DRIVER_HPP
