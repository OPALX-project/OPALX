/**
 * @file SolvePass.h
 * @brief Defines exhaustive tags for Cartesian PIC solve passes.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_SOLVE_PASS_H
#define OPALX_SPACE_CHARGE_PIC3D_SOLVE_PASS_H

#include "SpaceCharge/Pic3D/FieldComposer.h"

#include <cstdint>

namespace opalx::spacecharge {

    enum class SolvePass : std::uint8_t {
        DirectPrimary,
        DirectPrimaryWithImage,
        BinnedPrimary,
        BinnedImage,
        BinnedShiftedGreen
    };

    enum class BackendSolveRule : std::uint8_t { Standard, ShiftedGreen };

    /** @brief Complete behavior derived from one exhaustive pass tag. */
    template <typename T, unsigned Dim>
    struct SolvePassPolicy final {
        using ScatterGather = PicScatterGather<T, Dim>;
        using DepositKind   = typename ScatterGather::DepositKind;
        using ImagePolicy   = typename ScatterGather::ImagePolicy;

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
    [[nodiscard]] constexpr SolvePassPolicy<T, Dim> solvePassPolicy(SolvePass pass, T planeZ) {
        using Policy      = SolvePassPolicy<T, Dim>;
        using DepositKind = typename Policy::DepositKind;
        using ImagePolicy = typename Policy::ImagePolicy;

        switch (pass) {
            case SolvePass::DirectPrimary:
                return {};
            case SolvePass::DirectPrimaryWithImage:
                return {DepositKind::PrimaryAndImage,
                        ImagePolicy{true, planeZ},
                        BackendSolveRule::Standard,
                        false,
                        false,
                        1.0,
                        FieldSourceRule::Direct,
                        true};
            case SolvePass::BinnedPrimary:
                return {DepositKind::Primary,    {},   BackendSolveRule::Standard, true, true, 1.0,
                        FieldSourceRule::Direct, false};
            case SolvePass::BinnedImage:
                return {DepositKind::Image,
                        ImagePolicy{true, planeZ},
                        BackendSolveRule::Standard,
                        true,
                        true,
                        -1.0,
                        FieldSourceRule::Direct,
                        true};
            case SolvePass::BinnedShiftedGreen:
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

    [[nodiscard]] constexpr const char* solvePassLabel(SolvePass pass) {
        switch (pass) {
            case SolvePass::DirectPrimary:
                return "primary";
            case SolvePass::DirectPrimaryWithImage:
                return "primary-and-image";
            case SolvePass::BinnedPrimary:
                return "primary";
            case SolvePass::BinnedImage:
                return "image";
            case SolvePass::BinnedShiftedGreen:
                return "shifted-green";
        }
        return "unknown";
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_SOLVE_PASS_H
