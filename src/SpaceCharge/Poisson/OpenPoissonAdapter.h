#ifndef OPALX_OPEN_POISSON_ADAPTER_H
#define OPALX_OPEN_POISSON_ADAPTER_H

#include "SpaceCharge/Poisson/PoissonSolver.h"

#include "Physics/Physics.h"
#include "SpaceCharge/Poisson/PoissonSolverFactory.h"

#include <utility>

namespace opalx::spacecharge {

    /** @brief Adapts the 3D open FFT solver and its per-call shifted Green functions. */
    class OpenPoissonAdapter final : public PoissonSolver {
    public:
        OpenPoissonAdapter(PoissonSolverConfig config, PoissonFieldBinding fields)
            : PoissonSolver(std::move(config), fields, PoissonSolverType::Open) {
            rebuildImpl(fields);
        }

        [[nodiscard]] std::string_view name() const override { return "OPEN"; }
        [[nodiscard]] const PoissonSolverCapabilities& capabilities() const override {
            return capabilities_m;
        }
        [[nodiscard]] double couplingConstant() const override { return 1.0 / Physics::epsilon_0; }

    protected:
        void solveImpl(const PoissonSolveRequest& request) override {
            if (request.hasShiftedGreenFunction()) {
                backend_m->shiftedGreensFunction(*request.greenFunctionShift);
                backend_m->solve();
                backend_m->greensFunction();
            } else {
                backend_m->solve();
            }
        }
        void rebuildImpl(PoissonFieldBinding fields) override {
            auto& backend   = backend_m.emplace();
            auto parameters = detail::commonFftParameters();
            parameters.add("output_type", NativeBackend::SOL_AND_GRAD);
            parameters.add("algorithm", NativeBackend::HOCKNEY);
            parameters.add(
                    "greens_function", config_m.greenFunction == GreenFunctionType::Standard
                                               ? NativeBackend::STANDARD
                                               : NativeBackend::INTEGRATED);
            backend.mergeParameters(parameters);
            detail::bindFields(backend, fields);
        }

    private:
        static constexpr PoissonSolverCapabilities capabilities_m{
                .supportsShiftedGreenFunction   = true,
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = false,
                .debugDumpChargeBeforeSolve     = true,
                .debugDumpScalarAfterSolve      = true,
                .debugDumpVectorAfterSolve      = true};

        using NativeBackend = OpenSolver_t<double, 3>;
        std::optional<NativeBackend> backend_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_OPEN_POISSON_ADAPTER_H
