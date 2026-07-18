/**
 * @file CorrectionPlan.h
 * @brief Declares ordered Cartesian PIC correction-pass planning.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_CORRECTION_PLAN_H
#define OPALX_SPACE_CHARGE_PIC3D_CORRECTION_PLAN_H

#include "SpaceCharge/Pic3D/IterationPlan.h"
#include "SpaceCharge/Pic3D/SolvePass.h"
#include "SpaceCharge/SelfFieldConfig.h"
#include "SpaceCharge/SolveContext.h"

#include <array>
#include <cstddef>
#include <type_traits>

namespace opalx::spacecharge {

    /** @brief Validated correction policy and ordered passes for one solver call. */
    template <typename T, unsigned Dim>
    struct PreparedCorrection final {
        static constexpr std::size_t maximumPassCount = 2;

        using Pass = SolvePass<T, Dim>;

        [[nodiscard]] const Pass* begin() const { return passes.data(); }
        [[nodiscard]] const Pass* end() const { return passes.data() + passCount; }

        std::array<Pass, maximumPassCount> passes{};
        std::size_t passCount = 0;
        CorrectionRequest configuredCorrection;
        CorrectionRequest activeCorrection;
        std::size_t maximumSteps = 0;
        bool correctionExpired   = false;
        bool shiftedIgnored      = false;
    };

    /**
     * @brief Expands one validated step request into fixed Cartesian PIC solve passes.
     *
     * The plan copies immutable run configuration and retains no particle, mesh, field, parser,
     * or device state. prepare() validates the step-resolved request before an executor mutates
     * particle or field data.
     */
    template <typename T, unsigned Dim>
    class CorrectionPlan final {
        static_assert(Dim == 3, "CorrectionPlan currently supports Dim == 3 only.");

    public:
        using Pass        = SolvePass<T, Dim>;
        using Prepared    = PreparedCorrection<T, Dim>;
        using DepositKind = typename Pass::DepositKind;
        using ImagePolicy = typename Pass::ImagePolicy;

        CorrectionPlan(
                CorrectionConfig config, PoissonBackendKind backend, IterationKind iteration);

        CorrectionPlan(const CorrectionPlan&)            = delete;
        CorrectionPlan& operator=(const CorrectionPlan&) = delete;
        CorrectionPlan(CorrectionPlan&&)                 = delete;
        CorrectionPlan& operator=(CorrectionPlan&&)      = delete;

        /** @brief Validate @p requested for @p step and emit its ordered solve passes. */
        [[nodiscard]] Prepared prepare(const RequestedPhysics& requested, std::size_t step) const;

        [[nodiscard]] const CorrectionConfig& config() const { return config_m; }
        [[nodiscard]] PoissonBackendKind backend() const { return backend_m; }
        [[nodiscard]] IterationKind iteration() const { return iteration_m; }

    private:
        [[nodiscard]] Pass makePass(
                SolvePassKind kind, DepositKind depositKind, ImagePolicy imagePolicy,
                BackendSolveRule backendRule, bool suppressFieldDump, FieldOutputRule outputRule,
                double magneticSign, FieldSourceRule sourceRule,
                bool dumpDirichletPlaneAfter) const;
        void append(Prepared& prepared, const Pass& pass) const;
        void validateRequest(const RequestedPhysics& requested, bool correctionActive) const;

        CorrectionConfig config_m;
        PoissonBackendKind backend_m;
        IterationKind iteration_m;
    };

    static_assert(std::is_trivially_copyable_v<PreparedCorrection<double, 3>>);

    extern template class CorrectionPlan<double, 3>;

}  // namespace opalx::spacecharge

#include "SpaceCharge/Pic3D/CorrectionPlan.tpp"

#endif  // OPALX_SPACE_CHARGE_PIC3D_CORRECTION_PLAN_H
