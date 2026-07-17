/**
 * @file SelfFieldAlgorithm.h
 * @brief Common host-side interface implemented by self-field algorithms.
 */

#ifndef OPALX_SPACE_CHARGE_SELF_FIELD_ALGORITHM_H
#define OPALX_SPACE_CHARGE_SELF_FIELD_ALGORITHM_H

#include "SpaceCharge/SelfFieldDiagnostics.h"
#include "SpaceCharge/SolveContext.h"
#include "SpaceCharge/SolverCapabilities.h"

namespace opalx::spacecharge {

    /**
     * @brief Run-lifetime self-field algorithm selected during setup.
     *
     * Implementations may own meshes, backends, and persistent device scratch. execute() is a
     * host orchestration boundary; implementations must not retain the borrowed SolveContext.
     */
    class SelfFieldAlgorithm {
    public:
        virtual ~SelfFieldAlgorithm() = default;

        [[nodiscard]] virtual SolverCapabilities capabilities() const                  = 0;
        virtual void execute(SolveContext& context, SelfFieldDiagnostics& diagnostics) = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SELF_FIELD_ALGORITHM_H
