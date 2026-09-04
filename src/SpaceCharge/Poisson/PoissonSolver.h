/**
 * @file PoissonSolver.h
 * @brief Owns and dispatches the selected concrete IPPL Poisson solver.
 */

#ifndef OPALX_SPACE_CHARGE_POISSON_SOLVER_H
#define OPALX_SPACE_CHARGE_POISSON_SOLVER_H

#include "Ippl.h"
#include "Manager/datatypes.h"
#include "SpaceCharge/SpaceChargeConfig.h"

#include <array>
#include <concepts>
#include <cstddef>
#include <optional>
#include <string_view>
#include <variant>

namespace opalx::spacecharge {

    /** @brief Borrowed grid fields required by every supported Poisson solver. */
    struct PoissonFieldBinding {
        Field_t<3>* chargeDensity          = nullptr;
        VField_t<double, 3>* electricField = nullptr;
    };

    /** @brief Per-call changes to the Green function used by a Poisson solve. */
    struct PoissonSolveRequest {
        std::optional<ippl::Vector<double, 3>> greenFunctionShift;

        [[nodiscard]] bool hasShiftedGreenFunction() const {
            return greenFunctionShift.has_value();
        }
    };

    /** @brief Per-call output controls that do not affect the computed field. */
    struct PoissonSolveOptions {
        bool suppressFieldDump = false;
    };

    /** @brief Storage and normalization behavior of a selected Poisson solver. */
    struct PoissonSolverCapabilities {
        bool isNoOp                         = false;
        bool supportsShiftedGreenFunction   = false;
        bool normalizeChargeByCellVolume    = true;
        bool subtractNeutralizingBackground = true;
        bool debugDumpChargeBeforeSolve     = false;
        bool debugDumpScalarAfterSolve      = false;
        bool debugDumpVectorAfterSolve      = false;
    };

    struct PoissonSolverConfig {
        PoissonSolverType type          = PoissonSolverType::None;
        GreenFunctionType greenFunction = GreenFunctionType::Integrated;
        double p3mCutoff                = 0.0;
        std::array<FieldBoundaryCondition, 3> boundaryConditions{
                FieldBoundaryCondition::Open, FieldBoundaryCondition::Open,
                FieldBoundaryCondition::Open};
    };

    namespace detail {

        /**
         * @brief Compile-time contract for an OPALX adapter around an IPPL Poisson solver.
         *
         * An adapter owns one concrete backend, configures it from the flat solver snapshot, and
         * binds RHS before LHS during construction. Backend-specific temporary behavior belongs in
         * solve(); the facade depends only on this contract.
         */
        template <typename Backend>
        concept PoissonBackend =
                std::constructible_from<Backend, const PoissonSolverConfig&, PoissonFieldBinding>
                && requires(Backend& backend, const PoissonSolveRequest& request) {
                       { Backend::name() } -> std::same_as<std::string_view>;
                       {
                           Backend::capabilities()
                       } -> std::same_as<const PoissonSolverCapabilities&>;
                       { Backend::couplingConstant() } -> std::same_as<double>;
                       { backend.solve(request) } -> std::same_as<void>;
                   };

        class NullPoissonBackendAdapter final {
        public:
            NullPoissonBackendAdapter(
                    const PoissonSolverConfig& config, PoissonFieldBinding fields);

            [[nodiscard]] static std::string_view name() noexcept;
            [[nodiscard]] static const PoissonSolverCapabilities& capabilities() noexcept;
            [[nodiscard]] static double couplingConstant() noexcept;
            void solve(const PoissonSolveRequest& request);

        private:
            NullSolver_t<double, 3> solver_m;
        };

        class PeriodicFFTPoissonBackendAdapter final {
        public:
            PeriodicFFTPoissonBackendAdapter(
                    const PoissonSolverConfig& config, PoissonFieldBinding fields);

            [[nodiscard]] static std::string_view name() noexcept;
            [[nodiscard]] static const PoissonSolverCapabilities& capabilities() noexcept;
            [[nodiscard]] static double couplingConstant() noexcept;
            void solve(const PoissonSolveRequest& request);

        private:
            FFTSolver_t<double, 3> solver_m;
        };

        class OpenPoissonBackendAdapter final {
        public:
            OpenPoissonBackendAdapter(
                    const PoissonSolverConfig& config, PoissonFieldBinding fields);

            [[nodiscard]] static std::string_view name() noexcept;
            [[nodiscard]] static const PoissonSolverCapabilities& capabilities() noexcept;
            [[nodiscard]] static double couplingConstant() noexcept;
            void solve(const PoissonSolveRequest& request);

        private:
            OpenSolver_t<double, 3> solver_m;
        };

        class P3MPoissonBackendAdapter final {
        public:
            P3MPoissonBackendAdapter(const PoissonSolverConfig& config, PoissonFieldBinding fields);

            [[nodiscard]] static std::string_view name() noexcept;
            [[nodiscard]] static const PoissonSolverCapabilities& capabilities() noexcept;
            [[nodiscard]] static double couplingConstant() noexcept;
            void solve(const PoissonSolveRequest& request);

        private:
            FFTTruncatedGreenSolver_t<double, 3> solver_m;
        };

        static_assert(PoissonBackend<NullPoissonBackendAdapter>);
        static_assert(PoissonBackend<PeriodicFFTPoissonBackendAdapter>);
        static_assert(PoissonBackend<OpenPoissonBackendAdapter>);
        static_assert(PoissonBackend<P3MPoissonBackendAdapter>);

    }  // namespace detail

    /**
     * @brief OPALX-owned variant over the supported IPPL Poisson solver implementations.
     *
     * rebuildAfterLayoutChange() reconstructs the typed adapter after a layout change. Each adapter
     * binds the right-hand side before the left-hand side; reconstruction also resizes private FFT
     * fields and heFFTe plans that IPPL cannot refresh through its public interface. Construction
     * and execution are host-only; no concrete backend type crosses this boundary.
     *
     * To integrate another grid-based IPPL solver, add its parser/config enum mapping, implement a
     * PoissonBackend adapter, append that adapter to Backend, and add one construction case. The
     * Cartesian algorithm and SpaceChargeSolveContext remain unchanged. CG remains reserved until
     * its scalar potential and gradient bindings are implemented together.
     */
    class PoissonSolver final {
    public:
        PoissonSolver(PoissonSolverConfig config, PoissonFieldBinding fields);

        PoissonSolver(const PoissonSolver&)            = delete;
        PoissonSolver& operator=(const PoissonSolver&) = delete;
        PoissonSolver(PoissonSolver&&)                 = delete;
        PoissonSolver& operator=(PoissonSolver&&)      = delete;

        void solve(
                const PoissonSolveRequest& request = {}, const PoissonSolveOptions& options = {});
        void warmup();
        void rebuildAfterLayoutChange(PoissonFieldBinding fields);

        [[nodiscard]] std::string_view name() const;
        [[nodiscard]] const PoissonSolverCapabilities& capabilities() const;
        [[nodiscard]] double couplingConstant() const;

    private:
        void constructBackend(PoissonFieldBinding fields);

        using Backend = std::variant<
                std::monostate, detail::NullPoissonBackendAdapter,
                detail::PeriodicFFTPoissonBackendAdapter, detail::OpenPoissonBackendAdapter,
                detail::P3MPoissonBackendAdapter>;

        Backend backend_m;
        const PoissonSolverConfig config_m;
        PoissonFieldBinding fields_m;
        std::size_t runtimeSolveCount_m = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_POISSON_SOLVER_H
