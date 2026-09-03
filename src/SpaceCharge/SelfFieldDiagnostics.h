/**
 * @file SelfFieldDiagnostics.h
 * @brief Counter and output-cadence state for self-field algorithms.
 */

#ifndef OPALX_SPACE_CHARGE_SELF_FIELD_DIAGNOSTICS_H
#define OPALX_SPACE_CHARGE_SELF_FIELD_DIAGNOSTICS_H

#include <cstddef>

namespace opalx::spacecharge {

    struct SelfFieldDiagnosticSchedule {
        std::size_t binTableFrequency  = 0;
        std::size_t planeDumpFrequency = 0;
    };

    /** @brief Records successful solver work and compatibility output cadence. */
    class SelfFieldDiagnostics {
    public:
        explicit SelfFieldDiagnostics(SelfFieldDiagnosticSchedule schedule = {})
            : schedule_m(schedule) {}

        [[nodiscard]] std::size_t backendSolveCount() const { return backendSolveCount_m; }
        void completeBackendSolve() noexcept { ++backendSolveCount_m; }

        [[nodiscard]] std::size_t redistributionCount() const { return redistributionCount_m; }
        void completeRedistribution() noexcept { ++redistributionCount_m; }

        void recordBinCount(std::size_t count) noexcept {
            currentBinCount_m    = count;
            hasCurrentBinCount_m = true;
        }

        [[nodiscard]] bool hasCurrentBinCount() const { return hasCurrentBinCount_m; }
        [[nodiscard]] std::size_t currentBinCount() const { return currentBinCount_m; }

        [[nodiscard]] bool shouldPrintBinTable(long long step) const noexcept {
            return step >= 0 && schedule_m.binTableFrequency > 0
                   && static_cast<std::size_t>(step) % schedule_m.binTableFrequency == 0;
        }

        [[nodiscard]] bool shouldDumpPlane(long long step) const noexcept {
            return step >= 0 && schedule_m.planeDumpFrequency > 0
                   && static_cast<std::size_t>(step) % schedule_m.planeDumpFrequency == 0;
        }

    private:
        SelfFieldDiagnosticSchedule schedule_m;
        std::size_t backendSolveCount_m   = 0;
        std::size_t redistributionCount_m = 0;
        std::size_t currentBinCount_m     = 0;
        bool hasCurrentBinCount_m         = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SELF_FIELD_DIAGNOSTICS_H
