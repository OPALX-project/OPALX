/**
 * @file IpplPoissonAdapter.h
 * @brief Typed ownership boundary around concrete IPPL Poisson backends.
 */

#ifndef OPALX_SPACE_CHARGE_IPPL_POISSON_ADAPTER_H
#define OPALX_SPACE_CHARGE_IPPL_POISSON_ADAPTER_H

#include "SpaceCharge/Ippl/IpplPoissonTypes.h"
#include "SpaceCharge/SelfFieldConfig.h"

#include <memory>

namespace opalx::spacecharge {
    namespace detail {
        class IpplPoissonBackend;
    }

    /**
     * @brief Owns construction-time backend selection and all concrete IPPL solver access.
     *
     * The adapter borrows the temporary FieldSolverBase variant and stable field objects. It is
     * host-only and never enters a device kernel. refresh() must be called after a field-layout
     * change; each backend rebinds its right-hand side before its left-hand side.
     */
    class IpplPoissonAdapter final {
    public:
        [[nodiscard]] static std::unique_ptr<IpplPoissonAdapter> create(
                PoissonBackendKind kind, GreenFunctionKind greenFunction,
                Solver_t<double, 3>& solverStorage, IpplPoissonFields fields);

        [[nodiscard]] static const IpplPoissonCapabilities& capabilitiesFor(
                PoissonBackendKind kind);
        [[nodiscard]] static double couplingConstantFor(PoissonBackendKind kind);

        ~IpplPoissonAdapter();

        IpplPoissonAdapter(const IpplPoissonAdapter&)            = delete;
        IpplPoissonAdapter& operator=(const IpplPoissonAdapter&) = delete;
        IpplPoissonAdapter(IpplPoissonAdapter&&)                 = delete;
        IpplPoissonAdapter& operator=(IpplPoissonAdapter&&)      = delete;

        void solve(const IpplPoissonSolveRequest& request = {});
        void refresh(IpplPoissonFields fields);

        [[nodiscard]] const IpplPoissonCapabilities& capabilities() const;
        [[nodiscard]] double couplingConstant() const;

    private:
        explicit IpplPoissonAdapter(std::unique_ptr<detail::IpplPoissonBackend> backend);

        std::unique_ptr<detail::IpplPoissonBackend> backend_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_IPPL_POISSON_ADAPTER_H
