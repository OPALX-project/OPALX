/**
 * @file FieldComposer.h
 * @brief Declares Cartesian PIC field conversion, accumulation, and final gather operations.
 */

#ifndef OPALX_SPACE_CHARGE_PIC_FIELD_COMPOSER_H
#define OPALX_SPACE_CHARGE_PIC_FIELD_COMPOSER_H

#include "SpaceCharge/Pic/PicScatterGather.h"
#include "SpaceCharge/Pic/PicWorkspace.h"

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
    template <unsigned Dim>
    struct FieldCompositionPolicy final {
        std::array<double, Dim> meanMomentum{};
        double gamma               = 1.0;
        double magneticSign        = 1.0;
        FieldSourceRule sourceRule = FieldSourceRule::Direct;
    };

    static_assert(std::is_trivially_copyable_v<FieldCompositionPolicy<3>>);

    /**
     * @brief Converts electrostatic backend output to lab-frame fields and gathers final results.
     *
     * The composer owns no fields or particle data. Persistent accumulation and mirror scratch
     * live in @c PicWorkspace. Gather operations delegate interpolation to @c PicScatterGather;
     * the composer only chooses the source fields and replace semantics.
     *
     * All methods enclosing device kernels are public for CUDA builds. Device lambdas capture
     * only current views, scalar values, and small vector values obtained after the latest
     * particle or field-layout change.
     */
    template <typename T, unsigned Dim>
    class FieldComposer final {
        static_assert(Dim == 3, "FieldComposer currently supports Dim == 3 only.");

    public:
        using Workspace         = PicWorkspace<T, Dim>;
        using ScatterGather     = PicScatterGather<T, Dim>;
        using Vector            = typename Workspace::Vector;
        using VectorField       = typename Workspace::VectorField;
        using PositionAttribute = typename ScatterGather::PositionAttribute;
        using VectorAttribute   = typename ScatterGather::VectorAttribute;
        using Policy            = FieldCompositionPolicy<Dim>;

        /** @brief Clear the persistent lab-frame electric and magnetic accumulators. */
        void clearAccumulation(Workspace& workspace) const;

        /**
         * @brief Lorentz-convert and add one backend electric-field contribution.
         *
         * @c Direct reads the backend field at the current cell. @c ShiftedGreenImageZ first
         * builds a global z-mirror in persistent scratch, then applies the image-field component
         * signs before the same Lorentz conversion. Contributions are added in call order.
         */
        void accumulate(Workspace& workspace, const Policy& policy) const;

        /**
         * @brief Replace a writable particle electric field with the electrostatic backend field.
         *
         * Interpolation is delegated to @c PicScatterGather with replace semantics.
         */
        void gatherElectrostatic(
                ScatterGather& scatterGather, VectorAttribute& destination,
                const PositionAttribute& positions, Workspace& workspace) const;

        /**
         * @brief Replace writable particle E/B fields with the persistent accumulated fields.
         *
         * Electric interpolation is completed before magnetic interpolation, preserving the
         * existing final-gather order.
         */
        void gatherAccumulated(
                ScatterGather& scatterGather, VectorAttribute& electricDestination,
                VectorAttribute& magneticDestination, const PositionAttribute& positions,
                Workspace& workspace) const;
    };

    extern template class FieldComposer<double, 3>;

}  // namespace opalx::spacecharge

#include "SpaceCharge/Pic/FieldComposer.tpp"

#endif  // OPALX_SPACE_CHARGE_PIC_FIELD_COMPOSER_H
