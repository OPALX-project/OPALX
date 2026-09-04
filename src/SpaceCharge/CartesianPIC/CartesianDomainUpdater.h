/**
 * @file CartesianDomainUpdater.h
 * @brief Declares Cartesian PIC geometry, migration, and redistribution updates.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_DOMAIN_UPDATER_H
#define OPALX_SPACE_CHARGE_CARTESIAN_DOMAIN_UPDATER_H

#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/CartesianPIC/CartesianPICFieldStorage.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace opalx::spacecharge {

    class PoissonSolver;

    enum class DomainCoordinateFrame { Beam, Reference };

    struct CartesianBounds {
        ippl::Vector<double, 3> lower{0.0};
        ippl::Vector<double, 3> upper{0.0};
    };

    /** @brief Owns Cartesian mesh geometry, particle migration, and ORB scheduling. */
    class CartesianDomainUpdater final {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;
        using FieldStorage      = CartesianPICFieldStorage<double, 3>;
        using Orb               = ippl::OrthogonalRecursiveBisection<Field<double, 3>, double>;
        using FixedDomain       = BunchStateHandler::FixedCartesianDomainState;

        CartesianDomainUpdater(
                CartesianPICConfig config, std::span<ParticleContainer* const> particles);

        CartesianDomainUpdater(const CartesianDomainUpdater&)            = delete;
        CartesianDomainUpdater& operator=(const CartesianDomainUpdater&) = delete;

        /**
         * @brief Update geometry and ownership for one solve phase.
         *
         * Beam-frame updates use only the primary container. Reference-frame updates restore all
         * containers to one shared tracker-frame layout. Fixed bounds apply only to the beam phase.
         */
        [[nodiscard]] bool updateForSolve(
                DomainCoordinateFrame frame, const SpaceChargeSolveContext& context,
                const CorrectionConfig& correction, const FixedDomain* fixedDomain,
                FieldStorage& fieldStorage, PoissonSolver& poissonSolver);

    private:
        [[nodiscard]] CartesianBounds computeBounds(bool primaryOnly);
        void updateLayoutsAndMigrate(FieldStorage& fieldStorage, bool primaryOnly);
        void updateMoments(bool primaryOnly);
        [[nodiscard]] bool isRedistributionBlocked(
                std::span<const std::uint8_t> trackingActive) const;
        [[nodiscard]] bool loadIsImbalanced(double threshold);
        [[nodiscard]] FieldStorage::Extents targetExtents(const CorrectionConfig& correction) const;
        void extendImageBounds(CartesianBounds& bounds, const CorrectionConfig& correction) const;
        void expandBounds(
                CartesianBounds& bounds, bool applyEmissionStretch, double emittedFraction,
                std::size_t longitudinalExtent) const;
        void updateMeshGeometry(const CartesianBounds& bounds, FieldStorage& fieldStorage) const;
        [[nodiscard]] bool isRedistributionDue(std::size_t step) const;
        [[nodiscard]] bool redistribute(
                const SpaceChargeSolveContext& context, FieldStorage& fieldStorage);
        [[nodiscard]] ParticleContainer& primary() const { return *particles_m.front(); }

        CartesianPICConfig config_m;
        std::vector<ParticleContainer*> particles_m;
        Orb orb_m;
        std::vector<int> rankFlags_m;
        bool poissonRebuildRequired_m = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_DOMAIN_UPDATER_H
