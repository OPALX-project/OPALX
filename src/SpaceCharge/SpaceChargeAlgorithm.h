/**
 * @file SpaceChargeAlgorithm.h
 * @brief Common host-side interface implemented by space-charge algorithms.
 */

#ifndef OPALX_SPACE_CHARGE_ALGORITHM_H
#define OPALX_SPACE_CHARGE_ALGORITHM_H

#include "SpaceCharge/SpaceChargeSolveContext.h"

#include <cstddef>

namespace opalx::spacecharge {

    /** @brief Completed work and tracker-facing state from one solve. */
    struct SpaceChargeSolveResult {
        std::size_t backendSolves   = 0;
        std::size_t redistributions = 0;
        int reportedBins            = 1;
    };

    /**
     * @brief Run-lifetime space-charge algorithm selected during setup.
     *
     * Implementations may own meshes, field solvers, and persistent device scratch. Meshless or
     * tree algorithms implement this interface directly rather than entering Cartesian PIC.
     */
    class SpaceChargeAlgorithm {
    public:
        virtual ~SpaceChargeAlgorithm() = default;

        /**
         * @brief Compute one configured space-charge update.
         *
         * On entry, participating particles are in the tracker frame. A successful implementation
         * returns computed positions and fields in that frame. The context is borrowed only for
         * this call. Layout migration invalidates native device views, which must be reacquired.
         *
         * Exceptions propagate to the run-level handler and are terminal for the current run.
         * Temporary particle, mesh, field, backend, and frame state is unspecified after failure;
         * implementations do not provide rollback or retry guarantees.
         */
        [[nodiscard]] virtual SpaceChargeSolveResult solve(
                const SpaceChargeSolveContext& context) = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_ALGORITHM_H
