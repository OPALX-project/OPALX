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
         * R (metres), normalized momentum P (beta*gamma), E (V/m), and B (tesla) use tracker
         * axes on entry and return. Successful calls replace E/B on participating containers,
         * including zero fields for NONE, and preserve R/P/Q/dt apart from roundoff and particle
         * permutation/migration. Cartesian PIC participates with the primary container; FFT2D5
         * uses tracking-active containers. Other containers' field values are preserved.
         *
         * The context is borrowed only for this call. Spatial frame changes rotate P with R;
         * algorithm-specific Lorentz transformations are separate. Layout migration invalidates
         * native device views, which must be reacquired before later kernels.
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
