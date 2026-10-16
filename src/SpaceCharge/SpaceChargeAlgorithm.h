/**
 * @file SpaceChargeAlgorithm.h
 * @brief Tracker-facing interface for space-charge algorithms.
 */

#ifndef OPALX_SPACE_CHARGE_ALGORITHM_H
#define OPALX_SPACE_CHARGE_ALGORITHM_H

#include "SpaceCharge/SpaceChargeSolveContext.h"

#include <cstddef>
#include <string>

namespace opalx::spacecharge {

    /** @brief Work completed by one space-charge update. */
    struct SpaceChargeSolveResult {
        std::size_t backendSolves   = 0;  ///< Completed Poisson backend calls.
        std::size_t redistributions = 0;  ///< Completed ORB redistributions.
        int reportedBins = 1;  ///< Bin count reported to tracker diagnostics for this step.
    };

    /** @brief Run-lifetime space-charge algorithm selected during setup. */
    class SpaceChargeAlgorithm {
    public:
        virtual ~SpaceChargeAlgorithm() = default;

        /**
         * @brief Compute one configured space-charge update.
         *
         * R (metres), normalized P (beta-gamma), E (V/m), and B (tesla) use tracker axes on entry
         * and return. A successful call replaces E/B on participating containers and preserves
         * R/P/Q/dt, apart from roundoff and particle reordering during migration. TYPE=NONE returns
         * zero self-fields.
         *
         * @note Particle migration invalidates cached device views.
         */
        [[nodiscard]] virtual SpaceChargeSolveResult solve(
                const SpaceChargeSolveContext& context) = 0;

        /** @brief Algorithm-specific counters for the end-of-run log; empty if there are none. */
        [[nodiscard]] virtual std::string runSummary() const { return {}; }
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_ALGORITHM_H
