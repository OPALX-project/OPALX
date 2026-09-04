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
     * a host orchestration boundary; implementations must not retain the borrowed context.
     */
    class SpaceChargeAlgorithm {
    public:
        virtual ~SpaceChargeAlgorithm() = default;

        [[nodiscard]] virtual SpaceChargeCapabilities capabilities() const = 0;
        virtual void solve(
                SpaceChargeSolveContext& context, SpaceChargeDiagnostics& diagnostics) = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_ALGORITHM_H
