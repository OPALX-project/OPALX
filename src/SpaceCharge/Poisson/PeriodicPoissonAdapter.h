#ifndef OPALX_PERIODIC_POISSON_ADAPTER_H
#define OPALX_PERIODIC_POISSON_ADAPTER_H

#include "SpaceCharge/Poisson/PoissonSolver.h"

#include "Physics/Physics.h"
#include "SpaceCharge/Poisson/PoissonSolverFactory.h"

#include <utility>

namespace opalx::spacecharge {

    /** @brief Adapts the 3D periodic FFT Poisson solver. */
    class PeriodicPoissonAdapter final : public PoissonSolver {
    public:
        PeriodicPoissonAdapter(PoissonSolverConfig config, PoissonFieldBinding fields)
            : PoissonSolver(std::move(config), fields, PoissonSolverType::PeriodicFFT) {
            rebuildImpl(fields);
        }

        [[nodiscard]] std::string_view name() const override { return "FFT"; }
        [[nodiscard]] const PoissonSolverCapabilities& capabilities() const override {
            return capabilities_m;
        }
        [[nodiscard]] double couplingConstant() const override { return 1.0 / Physics::epsilon_0; }

    protected:
        void solveImpl(const PoissonSolveRequest&) override { backend_m->solve(); }
        void rebuildImpl(PoissonFieldBinding fields) override {
            auto& backend   = backend_m.emplace();
            auto parameters = detail::commonFftParameters();
            parameters.add("output_type", NativeBackend::GRAD);
            backend.mergeParameters(parameters);
            detail::bindFields(backend, fields);
        }

    private:
        static constexpr PoissonSolverCapabilities capabilities_m{
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = true,
                .debugDumpChargeBeforeSolve     = true,
                .debugDumpScalarAfterSolve      = true,
                .debugDumpVectorAfterSolve      = true};

        using NativeBackend = FFTSolver_t<double, 3>;
        std::optional<NativeBackend> backend_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_PERIODIC_POISSON_ADAPTER_H
