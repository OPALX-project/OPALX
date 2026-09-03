/**
 * @file SolverCapabilities.h
 * @brief Describes self-field algorithm features at the common host-side boundary.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVER_CAPABILITIES_H
#define OPALX_SPACE_CHARGE_SOLVER_CAPABILITIES_H

#include <cstdint>

namespace opalx::spacecharge {

    enum class SelfFieldAlgorithmKind : std::uint8_t { Pic3D, Pic2d5 };

    enum class ParticleSelectionPolicy : std::uint8_t { PrimaryOnly, AllTrackingActive };

    enum class CorrectionKind : std::uint8_t { None, ImageCharge, ShiftedGreen };

    /** @brief Features validated by SelfFieldSystem before algorithm dispatch. */
    struct SolverCapabilities {
        ParticleSelectionPolicy particleSelection = ParticleSelectionPolicy::PrimaryOnly;
        bool supportsBinning                      = false;
        bool supportsImageCharge                  = false;
        bool supportsShiftedGreen                 = false;
        bool supportsRedistribution               = false;
        bool supportsPotentialOutput              = false;

        [[nodiscard]] constexpr bool supports(CorrectionKind correction) const {
            switch (correction) {
                case CorrectionKind::None:
                    return true;
                case CorrectionKind::ImageCharge:
                    return supportsImageCharge;
                case CorrectionKind::ShiftedGreen:
                    return supportsShiftedGreen;
            }
            return false;
        }
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SOLVER_CAPABILITIES_H
