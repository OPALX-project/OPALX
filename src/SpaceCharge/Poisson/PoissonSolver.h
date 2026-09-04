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

    /**
     * @brief OPALX-owned variant over the supported IPPL Poisson solver implementations.
     *
     * rebuildAfterLayoutChange() reconstructs the typed backend after a layout change, then binds
     * the right-hand side before the left-hand side. Reconstructing also resizes backend-private
     * FFT fields and heFFTe plans that cannot be refreshed through the public IPPL interface.
     * Construction and execution are host-only; no concrete backend type crosses this boundary.
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

        [[nodiscard]] static const PoissonSolverCapabilities& capabilitiesFor(
                PoissonSolverType kind);
        [[nodiscard]] static double couplingConstantFor(PoissonSolverType kind);
        [[nodiscard]] const PoissonSolverCapabilities& capabilities() const {
            return capabilitiesFor(config_m.type);
        }
        [[nodiscard]] double couplingConstant() const { return couplingConstantFor(config_m.type); }

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
