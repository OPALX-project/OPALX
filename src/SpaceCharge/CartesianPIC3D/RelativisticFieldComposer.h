/**
 * @file RelativisticFieldComposer.h
 * @brief CartesianPIC3D field conversion, accumulation, and final gathering.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC3D_RELATIVISTIC_FIELD_COMPOSER_H
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC3D_RELATIVISTIC_FIELD_COMPOSER_H

#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DFieldStorage.h"
#include "SpaceCharge/CartesianPIC3D/ParticleMeshFieldTransfer.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace opalx::spacecharge {

    /** @brief Selects how backend electric-field samples enter field composition. */
    enum class FieldSourceRule : std::uint8_t { Direct, ShiftedGreenImageZ };

    /**
     * @brief Lorentz-conversion settings for one backend-field contribution.
     *
     * Mean momentum uses beta-gamma units in Cartesian solve axes. @c magneticSign is independent
     * of the field reflection selected by @c sourceRule.
     */
    struct FieldCompositionPolicy final {
        std::array<double, 3> meanMomentum{};
        double gamma               = 1.0;
        double magneticSign        = 1.0;
        FieldSourceRule sourceRule = FieldSourceRule::Direct;
    };

    static_assert(std::is_trivially_copyable_v<FieldCompositionPolicy>);

    /**
     * @brief Converts a rest-frame electric field to lab-frame E and B for one boost.
     *
     * The rest frame moves with @c velocity in Cartesian solve axes. The electric field
     * perpendicular to it scales with gamma, the parallel one is unchanged, and
     * B = sign * gamma / c^2 * v x E'. Mesh and particle contributions share this conversion.
     */
    struct LorentzBoost final {
        using Vector = ippl::Vector<double, 3>;

        Vector velocity{0.0};
        Vector velocityUnit{0.0};
        double gamma          = 1.0;
        double gammaMinusOne  = 0.0;
        double magneticFactor = 0.0;

        [[nodiscard]] static LorentzBoost from(const FieldCompositionPolicy& policy);

        KOKKOS_INLINE_FUNCTION void toLab(
                const Vector& electricRest, Vector& electricLab, Vector& magneticLab) const {
            const double electricDotVelocityUnit = electricRest.dot(velocityUnit);
            electricLab =
                    gamma * electricRest - gammaMinusOne * electricDotVelocityUnit * velocityUnit;
            magneticLab = magneticFactor * cross(velocity, electricRest);
        }
    };

    /**
     * @brief Converts and accumulates backend fields in Cartesian solve axes, then gathers them.
     *
     * Fields and particle data are borrowed; persistent scratch lives in
     * CartesianPIC3DFieldStorage.
     */
    class RelativisticFieldComposer final {
    public:
        using FieldStorage         = CartesianPIC3DFieldStorage<double, 3>;
        using ParticleMeshTransfer = ParticleMeshFieldTransfer;
        using Vector               = typename FieldStorage::Vector;
        using VectorField          = typename FieldStorage::VectorField;
        using PositionAttribute    = typename ParticleMeshTransfer::PositionAttribute;
        using VectorAttribute      = typename ParticleMeshTransfer::VectorAttribute;
        using VectorView           = typename VectorAttribute::view_type;
        using Policy               = FieldCompositionPolicy;

        // CUDA requires functions enclosing device lambdas to be public.

        /** @brief Clear the persistent electric and magnetic accumulators in the Cartesian solve
         * axes. */
        void clearAccumulation(FieldStorage& fieldStorage) const;

        /**
         * @brief Lorentz-convert and add one backend electric-field contribution.
         *
         * ShiftedGreenImageZ mirrors the field in z and applies image-field component signs first.
         * Contributions are added in call order.
         * @note Enqueues work on the default execution instance; synchronize before host access
         * or consumption on another execution instance.
         */
        void accumulate(FieldStorage& fieldStorage, const Policy& policy) const;

        /**
         * @brief Lorentz-convert per-particle rest-frame fields and add them to particle E/B.
         *
         * Used for particle-particle contributions, which bypass the mesh accumulators and must
         * be added after the final gather. Only the Direct source rule applies.
         * @param electricRest Rest-frame electric field, valid on [0, @p localCount).
         */
        void addParticleContribution(
                const VectorView& electricRest, VectorAttribute& electricDestination,
                VectorAttribute& magneticDestination, std::size_t localCount,
                const Policy& policy) const;

        /** @brief Replace a particle electric field with gathered backend values. */
        void gatherElectrostatic(
                ParticleMeshTransfer& particleMeshTransfer, VectorAttribute& destination,
                const PositionAttribute& positions, FieldStorage& fieldStorage) const;

        /**
         * @brief Replace particle E/B with one gathered rest-frame field, converted per particle.
         *
         * Equivalent to accumulate() followed by gatherAccumulated() for a single Direct
         * contribution, because CIC gathering and the Lorentz conversion are both linear. It
         * gathers one vector field instead of two and needs no mesh accumulators.
         */
        void gatherConverted(
                ParticleMeshTransfer& particleMeshTransfer, VectorAttribute& electricDestination,
                VectorAttribute& magneticDestination, const PositionAttribute& positions,
                FieldStorage& fieldStorage, std::size_t localCount, const Policy& policy) const;
        /** @brief Replace particle E/B with gathered accumulated fields. */
        void gatherAccumulated(
                ParticleMeshTransfer& particleMeshTransfer, VectorAttribute& electricDestination,
                VectorAttribute& magneticDestination, const PositionAttribute& positions,
                FieldStorage& fieldStorage) const;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC3D_RELATIVISTIC_FIELD_COMPOSER_H
