/**
 * @file SpaceChargeCapabilities.h
 * @brief Describes space-charge algorithm features at the common host-side boundary.
 */

#ifndef OPALX_SPACE_CHARGE_CAPABILITIES_H
#define OPALX_SPACE_CHARGE_CAPABILITIES_H

#include <cstdint>

namespace opalx::spacecharge {

    enum class SpaceChargeAlgorithmType : std::uint8_t { CartesianPIC, FFT2D5 };

    enum class ParticleSelectionMode : std::uint8_t { PrimaryOnly, AllTrackingActive };

    enum class SpaceChargeCorrectionType : std::uint8_t { None, ImageCharge, ShiftedGreen };

    /** @brief Features validated by SpaceChargeSolver before algorithm dispatch. */
    struct SpaceChargeCapabilities {
        ParticleSelectionMode particleSelection = ParticleSelectionMode::PrimaryOnly;
        bool supportsBinning                    = false;
        bool supportsImageCharge                = false;
        bool supportsShiftedGreen               = false;
        bool supportsRedistribution             = false;
        bool supportsPotentialOutput            = false;
        bool supportsFixedCartesianDomain       = false;

        [[nodiscard]] constexpr bool supports(SpaceChargeCorrectionType correction) const {
            switch (correction) {
                case SpaceChargeCorrectionType::None:
                    return true;
                case SpaceChargeCorrectionType::ImageCharge:
                    return supportsImageCharge;
                case SpaceChargeCorrectionType::ShiftedGreen:
                    return supportsShiftedGreen;
            }
            return false;
        }
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CAPABILITIES_H
