/**
 * @file PicDomainManager.h
 * @brief Declares Cartesian PIC domain updates and solver-owned redistribution.
 */

#ifndef OPALX_SPACE_CHARGE_PIC_DOMAIN_MANAGER_H
#define OPALX_SPACE_CHARGE_PIC_DOMAIN_MANAGER_H

#include "SpaceCharge/Pic/PicParticleDomainAdapter.h"
#include "SpaceCharge/Pic/PicWorkspace.h"
#include "SpaceCharge/SelfFieldConfig.h"
#include "SpaceCharge/SelfFieldDiagnostics.h"
#include "SpaceCharge/SolveContext.h"

#include <vector>

namespace opalx::spacecharge {

    class IpplPoissonAdapter;

    /** @brief Physical coordinate frame used to build the current PIC mesh. */
    enum class PicDomainFrame { Beam, Reference };

    /**
     * @brief Owns Cartesian mesh updates, layout changes, and ORB scheduling for 3D PIC.
     *
     * The manager stores immutable run configuration and reusable host communication scratch.
     * It does not own or borrow PartBunch, FieldSolver, or particle views. Every native object is
     * supplied explicitly for the duration of update().
     */
    class PicDomainManager final {
    public:
        using Orb = ippl::OrthogonalRecursiveBisection<Field<double, 3>, double>;

        explicit PicDomainManager(Pic3DConfig config);

        PicDomainManager(const PicDomainManager&)            = delete;
        PicDomainManager& operator=(const PicDomainManager&) = delete;
        PicDomainManager(PicDomainManager&&)                 = delete;
        PicDomainManager& operator=(PicDomainManager&&)      = delete;

        /**
         * @brief Rebuild physical geometry, migrate particles, and optionally redistribute.
         *
         * Beam-frame updates apply emission stretching and the configured ORB cadence. Reference
         * updates restore a mesh around reference-frame particles without either operation.
         * Image-domain extension and longitudinal resizing apply in both frames while active.
         */
        void update(
                PicDomainFrame frame, SolveContext& context, PicWorkspace<double, 3>& workspace,
                PicParticleDomainAdapter& particles, IpplPoissonAdapter& backend,
                SelfFieldDiagnostics& diagnostics);

    private:
        [[nodiscard]] bool imageChargeActive(std::size_t step) const;
        [[nodiscard]] PicWorkspace<double, 3>::Extents targetExtents(std::size_t step) const;
        void extendImageBounds(PicDomainBounds& bounds, std::size_t step) const;
        void expandBounds(
                PicDomainBounds& bounds, bool applyEmissionStretch, double emittedFraction,
                std::size_t longitudinalExtent) const;
        void updateGeometry(
                const PicDomainBounds& bounds, PicWorkspace<double, 3>& workspace) const;
        [[nodiscard]] bool redistributionDue(std::size_t step) const;
        [[nodiscard]] bool redistribute(
                SolveContext& context, PicWorkspace<double, 3>& workspace,
                PicParticleDomainAdapter& particles, SelfFieldDiagnostics& diagnostics);

        Pic3DConfig config_m;
        Orb orb_m;
        std::vector<int> rankFlags_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC_DOMAIN_MANAGER_H
