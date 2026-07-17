/**
 * @file SolverCapabilities.h
 * @brief Describes self-field algorithm features at the common host-side boundary.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVER_CAPABILITIES_H
#define OPALX_SPACE_CHARGE_SOLVER_CAPABILITIES_H

#include <cstdint>

namespace opalx::spacecharge {

    /** @brief Algorithm families selectable by the self-field factory. */
    enum class SelfFieldAlgorithmKind : std::uint8_t { Pic3D };

    /** @brief Per-call particle attributes exposed through ParticleSetView. */
    enum class ParticleAttribute : std::uint8_t {
        Position,
        Momentum,
        Charge,
        Mass,
        TimeStep,
        ElectricField,
        MagneticField,
        InvalidMask,
        Bin,
        Count
    };

    /** @brief Optional correction requested for one self-field solve. */
    enum class CorrectionKind : std::uint8_t { None, ImageCharge, ShiftedGreen };

    /**
     * @brief Compact set of required or available particle attributes.
     *
     * This value type is host-side metadata. It never contains particle data or a device view.
     */
    class ParticleAttributeSet {
    public:
        using storage_type = std::uint32_t;

        constexpr ParticleAttributeSet() = default;
        constexpr explicit ParticleAttributeSet(ParticleAttribute attribute)
            : bits_m(bit(attribute)) {}

        /** @brief Add an attribute to the set. */
        constexpr ParticleAttributeSet& add(ParticleAttribute attribute) {
            bits_m |= bit(attribute);
            return *this;
        }

        /** @brief Return true when the set contains an attribute. */
        [[nodiscard]] constexpr bool contains(ParticleAttribute attribute) const {
            return (bits_m & bit(attribute)) != 0;
        }

        /** @brief Return true when this set contains every attribute in @p other. */
        [[nodiscard]] constexpr bool contains(const ParticleAttributeSet& other) const {
            return (bits_m & other.bits_m) == other.bits_m;
        }

        /** @brief Return the raw bit representation for diagnostics. */
        [[nodiscard]] constexpr storage_type bits() const { return bits_m; }

    private:
        static constexpr storage_type bit(ParticleAttribute attribute) {
            return storage_type{1} << static_cast<storage_type>(attribute);
        }

        storage_type bits_m = 0;
    };

    /** @brief Combine two particle-attribute requirements. */
    [[nodiscard]] constexpr ParticleAttributeSet operator|(
            ParticleAttribute lhs, ParticleAttribute rhs) {
        ParticleAttributeSet result(lhs);
        result.add(rhs);
        return result;
    }

    /** @brief Combine an existing set with one particle-attribute requirement. */
    [[nodiscard]] constexpr ParticleAttributeSet operator|(
            ParticleAttributeSet lhs, ParticleAttribute rhs) {
        lhs.add(rhs);
        return lhs;
    }

    /**
     * @brief Features and data requirements of one concrete self-field algorithm.
     *
     * The system checks these values before dispatch. Concrete mesh, partition, and backend
     * types deliberately do not appear here.
     */
    struct SolverCapabilities {
        SelfFieldAlgorithmKind algorithm = SelfFieldAlgorithmKind::Pic3D;

        bool supportsBinning            = false;
        bool supportsImageCharge        = false;
        bool supportsShiftedGreen       = false;
        bool supportsRedistribution     = false;
        bool supportsMultipleContainers = false;
        bool supportsPotentialOutput    = false;

        ParticleAttributeSet requiredReadableAttributes;
        ParticleAttributeSet requiredWritableAttributes;

        /** @brief Return true when @p correction is implemented by this algorithm. */
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
