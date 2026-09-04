/**
 * @file CartesianDomainUpdater.h
 * @brief Declares Cartesian PIC domain updates and solver-owned redistribution.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_DOMAIN_UPDATER_H
#define OPALX_SPACE_CHARGE_CARTESIAN_DOMAIN_UPDATER_H

#include "SpaceCharge/CartesianPIC/CartesianPICFieldStorage.h"
#include "SpaceCharge/CartesianPIC/ParticleDomainOperations.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeDiagnostics.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"

#include <memory>
#include <vector>

class BunchStateHandler;

namespace opalx::spacecharge {

    class PoissonSolver;

    /** @brief Physical coordinate frame used to build the current PIC mesh. */
    enum class DomainCoordinateFrame { Beam, Reference };

    /**
     * @brief Owns Cartesian mesh updates, layout changes, and ORB scheduling for 3D PIC.
     *
     * The updater stores immutable run configuration, a shared read-only bunch-state handle,
     * and reusable host communication scratch. It does not borrow PartBunch, FieldSolver, or
     * particle views. Every native object is supplied explicitly for updateForSolve().
     */
    class CartesianDomainUpdater final {
    public:
        using Orb = ippl::OrthogonalRecursiveBisection<Field<double, 3>, double>;

        CartesianDomainUpdater(
                CartesianPICConfig config, std::shared_ptr<const BunchStateHandler> bunchState);

        CartesianDomainUpdater(const CartesianDomainUpdater&)            = delete;
        CartesianDomainUpdater& operator=(const CartesianDomainUpdater&) = delete;
        CartesianDomainUpdater(CartesianDomainUpdater&&)                 = delete;
        CartesianDomainUpdater& operator=(CartesianDomainUpdater&&)      = delete;

        /**
         * @brief Rebuild physical geometry, migrate particles, and optionally redistribute.
         *
         * Beam-frame updates apply emission stretching and the configured ORB cadence. Reference
         * updates restore a mesh around reference-frame particles without either operation. While
         * fixed state is active, its exact solve-frame bounds override envelope expansion and both
         * emission stretching and ORB are skipped.
         */
        void updateForSolve(
                DomainCoordinateFrame frame, SpaceChargeSolveContext& context,
                CartesianPICFieldStorage<double, 3>& fieldStorage,
                ParticleDomainOperations& particles, PoissonSolver& poissonSolver,
                SpaceChargeDiagnostics& diagnostics);

    private:
        [[nodiscard]] CartesianPICFieldStorage<double, 3>::Extents targetExtents(
                const SpaceChargeCorrectionRequest& correction) const;
        void extendImageBounds(
                CartesianBounds& bounds, const SpaceChargeCorrectionRequest& correction) const;
        void expandBounds(
                CartesianBounds& bounds, bool applyEmissionStretch, double emittedFraction,
                std::size_t longitudinalExtent) const;
        void updateMeshGeometry(
                const CartesianBounds& bounds,
                CartesianPICFieldStorage<double, 3>& fieldStorage) const;
        [[nodiscard]] bool isRedistributionDue(std::size_t step) const;
        [[nodiscard]] bool redistribute(
                SpaceChargeSolveContext& context, CartesianPICFieldStorage<double, 3>& fieldStorage,
                ParticleDomainOperations& particles);

        CartesianPICConfig config_m;
        std::shared_ptr<const BunchStateHandler> bunchState_m;
        Orb orb_m;
        std::vector<int> rankFlags_m;
        /** Tracks an ordinary nonthrowing ORB layout mutation until backend reconstruction. */
        bool poissonRebuildRequired_m = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_DOMAIN_UPDATER_H
