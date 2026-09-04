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
#include <cstddef>
#include <optional>
#include <string_view>
#include <variant>

namespace opalx::spacecharge {

    struct PoissonFieldBinding {
        Field_t<3>* chargeDensity          = nullptr;
        VField_t<double, 3>* electricField = nullptr;
    };

    struct PoissonSolveRequest {
        std::optional<ippl::Vector<double, 3>> greenFunctionShift;

        [[nodiscard]] bool hasShiftedGreenFunction() const {
            return greenFunctionShift.has_value();
        }
    };

    struct PoissonSolveOptions {
        bool suppressFieldDump = false;
    };

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

    /**
     * @brief Host-only variant over concrete IPPL Poisson backends.
     *
     * Native backends are fully emplaced before setRhs() initializes their fields and CUDA/FFT
     * resources. Layout changes reconstruct the backend and preserve RHS-before-LHS binding.
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
        void constructBackend();
        void bindBackendFields(PoissonFieldBinding fields);

        using NullBackend     = NullSolver_t<double, 3>;
        using PeriodicBackend = FFTSolver_t<double, 3>;
        using OpenBackend     = OpenSolver_t<double, 3>;
        using P3MBackend      = FFTTruncatedGreenSolver_t<double, 3>;
        using Backend =
                std::variant<std::monostate, NullBackend, PeriodicBackend, OpenBackend, P3MBackend>;

        Backend backend_m;
        const PoissonSolverConfig config_m;
        PoissonFieldBinding fields_m;
        std::size_t runtimeSolveCount_m = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_POISSON_SOLVER_H
