/**
 * @file SelfFieldRequestPolicy.h
 * @brief Resolves immutable run configuration into per-step physics requests.
 */

#ifndef OPALX_SPACE_CHARGE_SELF_FIELD_REQUEST_POLICY_H
#define OPALX_SPACE_CHARGE_SELF_FIELD_REQUEST_POLICY_H

#include "SpaceCharge/SelfFieldConfig.h"
#include "SpaceCharge/SolveContext.h"

#include <cstddef>

namespace opalx::spacecharge {

    /** @brief Tracker-owned snapshot of step-dependent self-field request rules. */
    class SelfFieldRequestPolicy {
    public:
        SelfFieldRequestPolicy() = default;
        explicit SelfFieldRequestPolicy(const SelfFieldConfig& config);

        [[nodiscard]] RequestedPhysics forStep(std::size_t step) const;
        [[nodiscard]] CorrectionRequest configuredCorrection() const { return correction_m; }

    private:
        CorrectionRequest correction_m;
        std::size_t correctionMaximumSteps_m = 0;
        bool useBinning_m                    = false;
        bool writePotential_m                = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SELF_FIELD_REQUEST_POLICY_H
