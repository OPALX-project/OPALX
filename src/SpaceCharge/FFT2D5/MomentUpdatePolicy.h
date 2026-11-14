/** @file MomentUpdatePolicy.h @brief Eligibility for deferred FFT2D5 tracker statistics. */
#ifndef OPALX_FFT2D5_MOMENT_UPDATE_POLICY_H
#define OPALX_FFT2D5_MOMENT_UPDATE_POLICY_H

#include <cstddef>

namespace opalx::spacecharge::fft2d5_detail {
    /** Conditions under which no consumer needs the three intermediate full moment updates.
     * The tracker converges eligibility across ranks before changing collective scheduling.
     */
    struct MomentUpdateEligibility {
        bool fft2d5                = false;
        bool midpoint              = false;
        bool boundaryControlled    = false;
        bool emissionActive        = false;
        bool globalProcesses       = false;
        bool spinEnabled           = false;
        bool solveActive           = false;
        std::size_t containerCount = 0;

        [[nodiscard]] bool mayDefer() const {
            return fft2d5 && midpoint && !boundaryControlled && !emissionActive && !globalProcesses
                   && !spinEnabled && solveActive && containerCount == 1;
        }
    };
}  // namespace opalx::spacecharge::fft2d5_detail

#endif
