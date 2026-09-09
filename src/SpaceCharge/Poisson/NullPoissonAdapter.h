#ifndef OPALX_NULL_POISSON_ADAPTER_H
#define OPALX_NULL_POISSON_ADAPTER_H

#include "SpaceCharge/Poisson/PoissonSolver.h"

#include "Physics/Physics.h"

#include <utility>

namespace opalx::spacecharge {

    /** @brief Adapts the configured no-op Poisson backend. */
    class NullPoissonAdapter final : public PoissonSolver {
    public:
        NullPoissonAdapter(PoissonSolverConfig config, PoissonFieldBinding fields)
            : PoissonSolver(std::move(config), fields, PoissonSolverType::None) {
            rebuildImpl(fields);
        }

        [[nodiscard]] std::string_view name() const override { return "NONE"; }
        [[nodiscard]] const PoissonSolverCapabilities& capabilities() const override {
            return capabilities_m;
        }
        [[nodiscard]] double couplingConstant() const override {
            return 1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
        }

    protected:
        void solveImpl(const PoissonSolveRequest&) override { backend_m->solve(); }
        void rebuildImpl(PoissonFieldBinding fields) override {
            auto& backend = backend_m.emplace();
            backend.mergeParameters(ippl::ParameterList{});
            detail::bindFields(backend, fields);
        }

    private:
        static constexpr PoissonSolverCapabilities capabilities_m{
                .isNoOp                         = true,
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = true};

        using NativeBackend = NullSolver_t<double, 3>;
        std::optional<NativeBackend> backend_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_NULL_POISSON_ADAPTER_H
