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

    /** @brief Charge-density input and electric-field output borrowed by a Poisson adapter. */
    struct PoissonFieldBinding {
        Field_t<3>* chargeDensity          = nullptr;
        VField_t<double, 3>* electricField = nullptr;
    };

    /** @brief Per-call options that can change the Poisson kernel. */
    struct PoissonSolveRequest {
        /** @brief Green-function displacement in metres in the current mesh axes. */
        std::optional<ippl::Vector<double, 3>> greenFunctionShift;

        [[nodiscard]] bool hasShiftedGreenFunction() const {
            return greenFunctionShift.has_value();
        }
    };

    /** @brief Per-call diagnostic options. */
    struct PoissonSolveOptions {
        bool suppressFieldDump = false;
    };

    /** @brief Backend properties used by deposition and diagnostics. */
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
     * The adapter borrows its fields. This base class handles validation, diagnostics, warmup, and
     * layout rebuilds; concrete adapters configure and run the native solver.
     */
    class PoissonSolver {
    public:
        virtual ~PoissonSolver() = default;

        PoissonSolver(const PoissonSolver&)            = delete;
        PoissonSolver& operator=(const PoissonSolver&) = delete;
        PoissonSolver(PoissonSolver&&)                 = delete;
        PoissonSolver& operator=(PoissonSolver&&)      = delete;

        /** @brief Run the native solver with the requested kernel and diagnostics. */
        void solve(
                const PoissonSolveRequest& request = {}, const PoissonSolveOptions& options = {});
        /** @brief Run a zero-RHS planning solve without advancing diagnostic numbering. */
        void warmup();
        /** @brief Rebind fields and rebuild native resources after a layout change. */
        void rebuildAfterLayoutChange(PoissonFieldBinding fields);

        [[nodiscard]] virtual std::string_view name() const                         = 0;
        [[nodiscard]] virtual const PoissonSolverCapabilities& capabilities() const = 0;
        /** @brief Charge-normalization factor required by the native backend. */
        [[nodiscard]] virtual double couplingConstant() const = 0;

    protected:
        PoissonSolver(
                PoissonSolverConfig config, PoissonFieldBinding fields,
                PoissonSolverType expectedType);

        virtual void solveImpl(const PoissonSolveRequest& request) = 0;
        /** @brief Reconstruct native resources and bind RHS before LHS. */
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
