/** @file DiagnosticOutput.h
 * @brief Internal validation of optional diagnostic output controls.
 */
#ifndef OPALX_DIAGNOSTIC_OUTPUT_H
#define OPALX_DIAGNOSTIC_OUTPUT_H

#include "Utilities/OpalException.h"

#include <charconv>
#include <string_view>
#include <system_error>

namespace opalx::diagnostics {
    /**
     * @brief Select a pre-push space-charge particle/field dump by global step.
     * @param value Value of OPALX_SC_FIELD_H5_STEPS; null or empty disables output.
     * @param step Zero-based global tracking step.
     * @return True through the inclusive nonnegative step limit. A negative limit
     * disables output. Only a complete signed decimal integer is accepted (no
     * whitespace, leading plus or trailing text).
     * @throws OpalException For invalid or overflowing limits; never enables
     * unlimited output as a fallback. This control does not change field solves.
     */
    inline bool spaceChargeFieldDumpRequested(const char* value, long long step) {
        if (value == nullptr || value[0] == '\0') {
            return false;
        }
        const std::string_view text(value);
        long long limit   = 0;
        const auto result = std::from_chars(text.data(), text.data() + text.size(), limit);
        if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
            throw OpalException(
                    "spaceChargeFieldDumpRequested",
                    "OPALX_SC_FIELD_H5_STEPS must be a signed decimal integer step limit.");
        }
        return limit >= 0 && step >= 0 && step <= limit;
    }
}  // namespace opalx::diagnostics

#endif
