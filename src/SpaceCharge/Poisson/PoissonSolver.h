/**
 * @file PoissonSolver.h
 * @brief Common lifecycle and request interface for 3D IPPL Poisson adapters.
 */

#ifndef OPALX_SPACE_CHARGE_POISSON_SOLVER_H
#define OPALX_SPACE_CHARGE_POISSON_SOLVER_H

#include "Ippl.h"
#include "Manager/datatypes.h"
#include "SpaceCharge/SpaceChargeConfig.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <string_view>

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

    /**
     * @brief Host-side lifecycle shared by concrete 3D IPPL Poisson adapters.
     *
     * Fields are borrowed for the adapter lifetime. Request validation, diagnostics and warmup
     * are shared; native setup and solving belong to the adapters. Exceptions terminate the run
     * and leave transient backend and field state unspecified.
     */
    class PoissonSolver {
    public:
        virtual ~PoissonSolver() = default;

        PoissonSolver(const PoissonSolver&)            = delete;
        PoissonSolver& operator=(const PoissonSolver&) = delete;
        PoissonSolver(PoissonSolver&&)                 = delete;
        PoissonSolver& operator=(PoissonSolver&&)      = delete;

        void solve(
                const PoissonSolveRequest& request = {}, const PoissonSolveOptions& options = {});
        void warmup();
        void rebuildAfterLayoutChange(PoissonFieldBinding fields);

        [[nodiscard]] virtual std::string_view name() const                         = 0;
        [[nodiscard]] virtual const PoissonSolverCapabilities& capabilities() const = 0;
        [[nodiscard]] virtual double couplingConstant() const                       = 0;

    protected:
        PoissonSolver(
                PoissonSolverConfig config, PoissonFieldBinding fields,
                PoissonSolverType expectedType);

        virtual void solveImpl(const PoissonSolveRequest& request) = 0;
        /** @brief Reconstruct native resources in place and bind RHS before LHS.
         *  @note Called from the concrete constructor or after the shared rebuild fence.
         */
        virtual void rebuildImpl(PoissonFieldBinding fields) = 0;

        const PoissonSolverConfig config_m;

    private:
        PoissonFieldBinding fields_m;
        std::size_t runtimeSolveCount_m = 0;
    };

    /** @brief Validate configuration and construct the selected 3D Poisson adapter. */
    [[nodiscard]] std::unique_ptr<PoissonSolver> makePoissonSolver(
            PoissonSolverConfig config, PoissonFieldBinding fields);

    namespace detail {

        inline ippl::ParameterList commonFftParameters() {
            ippl::ParameterList parameters;
            parameters.add("use_heffte_defaults", false);
            parameters.add("use_pencils", true);
            parameters.add("use_reorder", false);
            parameters.add("use_gpu_aware", true);
            parameters.add("comm", ippl::p2p_pl);
            parameters.add("r2c_direction", 0);
            return parameters;
        }

        /** @brief Bind native solver fields in the order required by IPPL initialization. */
        template <typename NativeBackend>
        void bindFields(NativeBackend& backend, PoissonFieldBinding fields) {
            backend.setRhs(*fields.chargeDensity);
            backend.setLhs(*fields.electricField);
        }

    }  // namespace detail

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_POISSON_SOLVER_H
