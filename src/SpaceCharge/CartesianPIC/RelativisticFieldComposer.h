/**
 * @file RelativisticFieldComposer.h
 * @brief Declares Cartesian PIC field conversion, accumulation, and final gather operations.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC_RELATIVISTIC_FIELD_COMPOSER_H
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC_RELATIVISTIC_FIELD_COMPOSER_H

#include "SpaceCharge/CartesianPIC/CartesianPICFieldStorage.h"
#include "SpaceCharge/CartesianPIC/ParticleMeshFieldTransfer.h"

#include <array>
#include <cstdint>
#include <type_traits>

namespace opalx::spacecharge {

    /** @brief Selects how backend electric-field samples enter field composition. */
    enum class FieldSourceRule : std::uint8_t { Direct, ShiftedGreenImageZ };

    /**
     * @brief Device-safe value policy for one backend-field contribution.
     *
     * Momentum is normalized as beta-gamma. The magnetic sign remains independent of the source
     * rule because the explicit image-charge and shifted-Green paths share the same magnetic sign
     * convention while obtaining their electric samples differently.
     */
    struct FieldCompositionPolicy final {
        std::array<double, 3> meanMomentum{};
        double gamma               = 1.0;
        double magneticSign        = 1.0;
        FieldSourceRule sourceRule = FieldSourceRule::Direct;
    };

    static_assert(std::is_trivially_copyable_v<FieldCompositionPolicy>);

    /**
     * @brief Converts electrostatic backend output to lab-frame fields and gathers final results.
     *
     * The composer owns no fields or particle data. Persistent accumulation and mirror scratch
     * live in @c CartesianPICFieldStorage. Gather operations delegate interpolation to @c
     * ParticleMeshFieldTransfer; the composer only chooses the source fields and replace semantics.
     *
     * All methods enclosing device kernels are public for CUDA builds. Device lambdas capture
     * only current views, scalar values, and small vector values obtained after the latest
     * particle or field-layout change.
     */
    class RelativisticFieldComposer final {
    public:
        using FieldStorage         = CartesianPICFieldStorage<double, 3>;
        using ParticleMeshTransfer = ParticleMeshFieldTransfer;
        using Vector               = typename FieldStorage::Vector;
        using VectorField          = typename FieldStorage::VectorField;
        using PositionAttribute    = typename ParticleMeshTransfer::PositionAttribute;
        using VectorAttribute      = typename ParticleMeshTransfer::VectorAttribute;
        using Policy               = FieldCompositionPolicy;

        /** @brief Clear the persistent lab-frame electric and magnetic accumulators. */
        void clearAccumulation(FieldStorage& fieldStorage) const;

        /**
         * @brief Lorentz-convert and add one backend electric-field contribution.
         *
         * @c Direct reads the backend field at the current cell. @c ShiftedGreenImageZ first
         * builds a global z-mirror in persistent scratch, then applies the image-field component
         * signs before the same Lorentz conversion. Contributions are added in call order.
         */
        void accumulate(FieldStorage& fieldStorage, const Policy& policy) const;

        /**
         * @brief Replace a writable particle electric field with the electrostatic backend field.
         *
         * Interpolation is delegated to @c ParticleMeshFieldTransfer with replace semantics.
         */
        void gatherElectrostatic(
                ParticleMeshTransfer& particleMeshTransfer, VectorAttribute& destination,
                const PositionAttribute& positions, FieldStorage& fieldStorage) const;

        /**
         * @brief Replace writable particle E/B fields with the persistent accumulated fields.
         *
         * Electric interpolation is completed before magnetic interpolation, preserving the
         * existing final-gather order.
         */
        void gatherAccumulated(
                ParticleMeshTransfer& particleMeshTransfer, VectorAttribute& electricDestination,
                VectorAttribute& magneticDestination, const PositionAttribute& positions,
                FieldStorage& fieldStorage) const;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC_RELATIVISTIC_FIELD_COMPOSER_H
