/**
 * @file SolvePass.h
 * @brief Defines host-side policies for one Cartesian PIC solve pass.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_SOLVE_PASS_H
#define OPALX_SPACE_CHARGE_PIC3D_SOLVE_PASS_H

#include "SpaceCharge/Pic3D/FieldComposer.h"

#include <cstdint>
#include <type_traits>

namespace opalx::spacecharge {

    /** @brief Physical role of one ordered Cartesian PIC solve pass. */
    enum class SolvePassKind : std::uint8_t { Primary, PrimaryAndImage, ImageCharge, ShiftedGreen };

    /** @brief Backend request selected for one solve pass. */
    enum class BackendSolveRule : std::uint8_t { Standard, ShiftedGreen };

    /** @brief How the backend field enters the particle result. */
    enum class FieldOutputRule : std::uint8_t { DirectGather, LorentzAccumulation };

    /**
     * @brief Trivially-copyable policy for one ordered Cartesian PIC solve pass.
     *
     * The record contains no particle selection, field, device view, or owning storage. A solve
     * unit supplies its selection separately, while this policy selects charge deposition,
     * backend behavior, field composition, and diagnostics for that unit.
     */
    template <typename T, unsigned Dim>
    struct SolvePass final {
        using ScatterGather = PicScatterGather<T, Dim>;
        using DepositKind   = typename ScatterGather::DepositKind;
        using ImagePolicy   = typename ScatterGather::ImagePolicy;

        /** @brief Stable solve-pass diagnostic label. */
        [[nodiscard]] static constexpr const char* solveLabel(SolvePassKind passKind) {
            switch (passKind) {
                case SolvePassKind::Primary:
                    return "primary";
                case SolvePassKind::PrimaryAndImage:
                    return "primary-and-image";
                case SolvePassKind::ImageCharge:
                    return "image";
                case SolvePassKind::ShiftedGreen:
                    return "shifted-green";
            }
            return "unknown";
        }

        /** @brief Stable field-composition diagnostic label. */
        [[nodiscard]] static constexpr const char* compositionLabel(SolvePassKind passKind) {
            switch (passKind) {
                case SolvePassKind::Primary:
                    return "primary-accumulation";
                case SolvePassKind::PrimaryAndImage:
                    return "legacy-gather";
                case SolvePassKind::ImageCharge:
                    return "image-accumulation";
                case SolvePassKind::ShiftedGreen:
                    return "shifted-green-accumulation";
            }
            return "unknown";
        }

        [[nodiscard]] constexpr const char* solveLabel() const { return solveLabel(kind); }
        [[nodiscard]] constexpr const char* compositionLabel() const {
            return compositionLabel(kind);
        }

        SolvePassKind kind      = SolvePassKind::Primary;
        DepositKind depositKind = DepositKind::Primary;
        ImagePolicy imagePolicy;
        BackendSolveRule backendRule = BackendSolveRule::Standard;
        bool suppressFieldDump       = false;
        FieldOutputRule outputRule   = FieldOutputRule::DirectGather;
        double magneticSign          = 1.0;
        FieldSourceRule sourceRule   = FieldSourceRule::Direct;
        double planeZ                = 0.0;
        bool dumpDirichletPlaneAfter = false;
    };

    static_assert(
            std::is_trivially_copyable_v<PicScatterGather<double, 3>::ImagePolicy>,
            "ImagePolicy must remain a device-independent value policy.");
    static_assert(
            std::is_trivially_copyable_v<SolvePass<double, 3>>,
            "SolvePass must remain a host value record without owning state.");

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_SOLVE_PASS_H
