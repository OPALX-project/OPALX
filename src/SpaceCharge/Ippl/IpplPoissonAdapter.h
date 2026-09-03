/**
 * @file IpplPoissonAdapter.h
 * @brief Owns and dispatches the selected concrete IPPL Poisson solver.
 */

#ifndef OPALX_SPACE_CHARGE_IPPL_POISSON_ADAPTER_H
#define OPALX_SPACE_CHARGE_IPPL_POISSON_ADAPTER_H

#include "SpaceCharge/Ippl/IpplPoissonTypes.h"
#include "SpaceCharge/SelfFieldConfig.h"

#include <array>
#include <cstddef>
#include <variant>

namespace opalx::spacecharge {

    struct IpplPoissonBackendConfig {
        PoissonBackendKind kind         = PoissonBackendKind::None;
        GreenFunctionKind greenFunction = GreenFunctionKind::Integrated;
        double p3mCutoff                = 0.0;
        std::array<BoundaryConditionKind, 3> boundaryConditions{
                BoundaryConditionKind::Open, BoundaryConditionKind::Open,
                BoundaryConditionKind::Open};
    };

    /**
     * @brief OPALX-owned variant over the supported IPPL Poisson solver implementations.
     *
     * refresh() rebinds the right-hand side before the left-hand side after every layout change.
     * Construction and execution are host-only; no concrete backend type crosses this boundary.
     */
    class IpplPoissonAdapter final {
    public:
        IpplPoissonAdapter(IpplPoissonBackendConfig config, IpplPoissonFields fields);

        IpplPoissonAdapter(const IpplPoissonAdapter&)            = delete;
        IpplPoissonAdapter& operator=(const IpplPoissonAdapter&) = delete;
        IpplPoissonAdapter(IpplPoissonAdapter&&)                 = delete;
        IpplPoissonAdapter& operator=(IpplPoissonAdapter&&)      = delete;

        void solve(
                const IpplPoissonSolveRequest& request = {},
                const IpplPoissonSolveOptions& options = {});
        void warmup();
        void refresh(IpplPoissonFields fields);

        [[nodiscard]] static const IpplPoissonCapabilities& capabilitiesFor(
                PoissonBackendKind kind);
        [[nodiscard]] static double couplingConstantFor(PoissonBackendKind kind);
        [[nodiscard]] const IpplPoissonCapabilities& capabilities() const {
            return capabilitiesFor(kind_m);
        }
        [[nodiscard]] double couplingConstant() const { return couplingConstantFor(kind_m); }

    private:
        using NullBackend     = NullSolver_t<double, 3>;
        using PeriodicBackend = FFTSolver_t<double, 3>;
        using OpenBackend     = OpenSolver_t<double, 3>;
        using P3MBackend      = FFTTruncatedGreenSolver_t<double, 3>;
        using Backend =
                std::variant<std::monostate, NullBackend, PeriodicBackend, OpenBackend, P3MBackend>;

        Backend backend_m;
        IpplPoissonFields fields_m;
        PoissonBackendKind kind_m;
        std::size_t runtimeSolveCount_m = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_IPPL_POISSON_ADAPTER_H
