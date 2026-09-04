/**
 * @file CartesianPICPass.h
 * @brief Defines exhaustive tags for Cartesian PIC solve passes.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC_PASS_H
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC_PASS_H

#include "SpaceCharge/CartesianPIC/RelativisticFieldComposer.h"

#include <cstdint>

namespace opalx::spacecharge {

    /**
     * @brief Physical solve passes supported by Cartesian PIC.
     *
     * Direct image charge combines primary and mirrored charge in one Poisson solve. Binned image
     * and shifted-Green corrections remain separate second passes so their field contributions can
     * be transformed and accumulated with the required component signs.
     */
    enum class CartesianPICPass : std::uint8_t {
        DirectPrimary,
        DirectPrimaryWithImage,
        BinnedPrimary,
        BinnedImage,
        BinnedShiftedGreen
    };

    enum class BackendSolveRule : std::uint8_t { Standard, ShiftedGreen };

    /** @brief Complete behavior derived from one exhaustive pass tag. */
    template <typename T, unsigned Dim>
    struct CartesianPICPassProperties final {
        using ParticleMeshTransfer = ParticleMeshFieldTransfer<T, Dim>;
        using DepositKind          = typename ParticleMeshTransfer::DepositKind;
        using ImagePolicy          = typename ParticleMeshTransfer::ImagePolicy;

        DepositKind depositKind = DepositKind::Primary;
        ImagePolicy imagePolicy;
        BackendSolveRule backendRule = BackendSolveRule::Standard;
        bool binned                  = false;
        bool suppressFieldDump       = false;
        double magneticSign          = 1.0;
        FieldSourceRule sourceRule   = FieldSourceRule::Direct;
        bool dumpDirichletPlaneAfter = false;
    };

    template <typename T, unsigned Dim>
    [[nodiscard]] constexpr CartesianPICPassProperties<T, Dim> cartesianPICPassProperties(
            CartesianPICPass pass, T planeZ) {
        using Policy      = CartesianPICPassProperties<T, Dim>;
        using DepositKind = typename Policy::DepositKind;
        using ImagePolicy = typename Policy::ImagePolicy;

        switch (pass) {
            case CartesianPICPass::DirectPrimary:
                return {};
            case CartesianPICPass::DirectPrimaryWithImage:
                return {DepositKind::PrimaryAndImage,
                        ImagePolicy{true, planeZ},
                        BackendSolveRule::Standard,
                        false,
                        false,
                        1.0,
                        FieldSourceRule::Direct,
                        true};
            case CartesianPICPass::BinnedPrimary:
                return {DepositKind::Primary,    {},   BackendSolveRule::Standard, true, true, 1.0,
                        FieldSourceRule::Direct, false};
            case CartesianPICPass::BinnedImage:
                return {DepositKind::Image,
                        ImagePolicy{true, planeZ},
                        BackendSolveRule::Standard,
                        true,
                        true,
                        -1.0,
                        FieldSourceRule::Direct,
                        true};
            case CartesianPICPass::BinnedShiftedGreen:
                return {DepositKind::Primary,
                        {},
                        BackendSolveRule::ShiftedGreen,
                        true,
                        false,
                        -1.0,
                        FieldSourceRule::ShiftedGreenImageZ,
                        false};
        }
        return {};
    }

    [[nodiscard]] constexpr const char* cartesianPICPassLabel(CartesianPICPass pass) {
        switch (pass) {
            case CartesianPICPass::DirectPrimary:
                return "primary";
            case CartesianPICPass::DirectPrimaryWithImage:
                return "primary-and-image";
            case CartesianPICPass::BinnedPrimary:
                return "primary";
            case CartesianPICPass::BinnedImage:
                return "image";
            case CartesianPICPass::BinnedShiftedGreen:
                return "shifted-green";
        }
        return "unknown";
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC_PASS_H
