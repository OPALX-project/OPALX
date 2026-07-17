/**
 * @file SelfFieldSystem.h
 * @brief Owns one configured self-field algorithm and validates solve requests.
 */

#ifndef OPALX_SPACE_CHARGE_SELF_FIELD_SYSTEM_H
#define OPALX_SPACE_CHARGE_SELF_FIELD_SYSTEM_H

#include "SpaceCharge/SelfFieldAlgorithm.h"
#include "SpaceCharge/SelfFieldConfig.h"

#include <memory>

namespace opalx::spacecharge {

    /** @brief Stable tracker-facing entry point for all self-field algorithms. */
    class SelfFieldSystem {
    public:
        SelfFieldSystem(SelfFieldConfig config, std::unique_ptr<SelfFieldAlgorithm> algorithm);

        SelfFieldSystem(const SelfFieldSystem&)            = delete;
        SelfFieldSystem& operator=(const SelfFieldSystem&) = delete;
        SelfFieldSystem(SelfFieldSystem&&)                 = delete;
        SelfFieldSystem& operator=(SelfFieldSystem&&)      = delete;

        /** @brief Validate and dispatch one borrowed solve context. */
        void solve(SolveContext& context);

        [[nodiscard]] const SelfFieldConfig& config() const { return config_m; }
        [[nodiscard]] const SolverCapabilities& capabilities() const { return capabilities_m; }
        [[nodiscard]] SelfFieldDiagnostics& diagnostics() { return diagnostics_m; }
        [[nodiscard]] const SelfFieldDiagnostics& diagnostics() const { return diagnostics_m; }
        void setDiagnosticSink(SelfFieldDiagnosticSink* sink) { diagnostics_m.setSink(sink); }

    private:
        void validateConfiguration() const;
        void validateContext(const SolveContext& context) const;

        SelfFieldConfig config_m;
        std::unique_ptr<SelfFieldAlgorithm> algorithm_m;
        SolverCapabilities capabilities_m;
        SelfFieldDiagnostics diagnostics_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SELF_FIELD_SYSTEM_H
