/**
 * @file IpplPoissonAdapter.h
 * @brief Typed ownership boundary around concrete IPPL Poisson backends.
 */

#ifndef OPALX_SPACE_CHARGE_IPPL_POISSON_ADAPTER_H
#define OPALX_SPACE_CHARGE_IPPL_POISSON_ADAPTER_H

#include "SpaceCharge/Ippl/IpplPoissonTypes.h"
#include "SpaceCharge/SelfFieldConfig.h"

#include <array>
#include <memory>

namespace opalx::spacecharge {
    namespace detail {
        class IpplPoissonBackend;
    }

    /**
     * @brief Owns construction-time backend selection and all concrete IPPL solver access.
     *
     * The adapter owns the selected IPPL solver variant and borrows stable field objects. It is
     * host-only and never enters a device kernel. refresh() must be called after a field-layout
     * change; each backend rebinds its right-hand side before its left-hand side.
     */
    class IpplPoissonAdapter final {
    public:
        [[nodiscard]] static std::unique_ptr<IpplPoissonAdapter> create(
                PoissonBackendKind kind, GreenFunctionKind greenFunction, IpplPoissonFields fields);

        [[nodiscard]] static const IpplPoissonCapabilities& capabilitiesFor(
                PoissonBackendKind kind);
        [[nodiscard]] static double couplingConstantFor(PoissonBackendKind kind);

        ~IpplPoissonAdapter();

        IpplPoissonAdapter(const IpplPoissonAdapter&)            = delete;
        IpplPoissonAdapter& operator=(const IpplPoissonAdapter&) = delete;
        IpplPoissonAdapter(IpplPoissonAdapter&&)                 = delete;
        IpplPoissonAdapter& operator=(IpplPoissonAdapter&&)      = delete;

        void solve(
                const IpplPoissonSolveRequest& request = {},
                const IpplPoissonSolveOptions& options = {});
        /** @brief Clear deposited charge and perform the construction-time planning solve. */
        void warmup();
        /** @brief Apply IPPL boundary conditions to the backend potential field when required. */
        void setPotentialBoundaryConditions(
                const std::array<BoundaryConditionKind, 3>& boundaryConditions);
        /** @brief Compatibility overload for an already translated IPPL boundary container. */
        void setPotentialBoundaryConditions(ippl::BConds<Field_t<3>, 3> boundaryConditions);
        void refresh(IpplPoissonFields fields);

        [[nodiscard]] const IpplPoissonCapabilities& capabilities() const;
        [[nodiscard]] double couplingConstant() const;

    private:
        IpplPoissonAdapter(
                PoissonBackendKind kind, GreenFunctionKind greenFunction, IpplPoissonFields fields);

        // Backends borrow the active alternative, so storage must outlive backend_m.
        Solver_t<double, 3> solverStorage_m;
        std::unique_ptr<detail::IpplPoissonBackend> backend_m;
        IpplPoissonFields fields_m;
        PoissonBackendKind kind_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_IPPL_POISSON_ADAPTER_H
