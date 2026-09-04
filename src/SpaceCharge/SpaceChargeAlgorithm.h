/**
 * @file SpaceChargeAlgorithm.h
 * @brief Common host-side interface implemented by space-charge algorithms.
 */

#ifndef OPALX_SPACE_CHARGE_ALGORITHM_H
#define OPALX_SPACE_CHARGE_ALGORITHM_H

#include "SpaceCharge/SpaceChargeCapabilities.h"
#include "SpaceCharge/SpaceChargeDiagnostics.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"

namespace opalx::spacecharge {

    /**
     * @brief Run-lifetime space-charge algorithm selected during setup.
     *
     * Implementations may own meshes, Poisson solvers, and persistent device scratch. solve() is
     * a host orchestration boundary selected by SpaceChargeConfig and registered in
     * SpaceChargeSolverFactory. Meshless or tree algorithms implement this interface directly
     * rather than entering the Cartesian PIC Poisson layer.
     */
    class SpaceChargeAlgorithm {
    public:
        virtual ~SpaceChargeAlgorithm() = default;

        [[nodiscard]] virtual SpaceChargeCapabilities capabilities() const = 0;

        /**
         * @brief Compute one configured space-charge update.
         *
         * On entry, selected particle positions and momenta are in the tracker frame and E/B are
         * writable destinations. A successful implementation returns primary positions and all
         * computed fields in that same frame. The context, its particle bindings, transforms, and
         * requests are borrowed only for this call and must not be retained. Layout migration can
         * invalidate native device views, so every later kernel must reacquire them from the stable
         * containers or fields. Record diagnostics only after the corresponding work succeeds.
         *
         * Exceptions propagate to the run-level handler and are terminal for the current run.
         * Temporary particle, mesh, field, backend, and frame state is unspecified after failure;
         * implementations do not provide rollback or retry guarantees.
         */
        virtual void solve(
                SpaceChargeSolveContext& context, SpaceChargeDiagnostics& diagnostics) = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_ALGORITHM_H
