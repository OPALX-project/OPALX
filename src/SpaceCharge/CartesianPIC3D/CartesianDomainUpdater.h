/**
 * @file CartesianDomainUpdater.h
 * @brief CartesianPIC3D geometry, migration, and redistribution updates.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_DOMAIN_UPDATER_H
#define OPALX_SPACE_CHARGE_CARTESIAN_DOMAIN_UPDATER_H

#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DFieldStorage.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace opalx::spacecharge {

    class PoissonSolver;

    /** @brief Coordinate frame used to compute the next particle and mesh domain. */
    enum class DomainCoordinateFrame { Beam, Reference };

    /** @brief Axis-aligned domain bounds in metres in the selected coordinate frame. */
    struct CartesianBounds {
        ippl::Vector<double, 3> lower{0.0};
        ippl::Vector<double, 3> upper{0.0};
    };

    /** @brief Owns Cartesian mesh geometry, particle migration, and ORB scheduling. */
    class CartesianDomainUpdater final {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;
        using FieldStorage      = CartesianPIC3DFieldStorage<double, 3>;
        using Orb               = ippl::OrthogonalRecursiveBisection<Field<double, 3>, double>;
        using FixedDomain       = BunchStateHandler::FixedCartesianDomainState;

        CartesianDomainUpdater(
                CartesianPIC3DConfig config, std::span<ParticleContainer* const> particles);

        CartesianDomainUpdater(const CartesianDomainUpdater&)            = delete;
        CartesianDomainUpdater& operator=(const CartesianDomainUpdater&) = delete;

        /**
         * @brief Update geometry and ownership for one solve phase.
         *
         * Beam-frame updates use only the first (primary) container. Reference-frame updates place
         * all containers on one tracker-frame layout. Fixed bounds apply only in the beam frame.
         *
         * @return Whether ORB changed the particle decomposition.
         */
        [[nodiscard]] bool updateForSolve(
                DomainCoordinateFrame frame, const SpaceChargeSolveContext& context,
                const DirichletPlaneConfig& dirichletPlane, const FixedDomain* fixedDomain,
                FieldStorage& fieldStorage, PoissonSolver& poissonSolver);

    private:
        [[nodiscard]] CartesianBounds computeBounds(bool primaryOnly);
        /**
         * @brief Update the particle layouts to the current mesh and migrate particles.
         *
         * Beam-frame updates handle only the primary container. A mesh-tied P3M cutoff is set
         * from the largest spacing there, since the pair stage runs on these layouts. Reference
         * frame layouts only migrate particles and use the smallest spacing instead. Either value
         * is reduced to the largest cutoff that fits every rank region, see
         * CartesianDomain::largestFittingOverlapCutoff().
         */
        void updateLayoutsAndMigrate(FieldStorage& fieldStorage, DomainCoordinateFrame frame);
        void updateMoments(bool primaryOnly);
        [[nodiscard]] bool isRedistributionBlocked(
                std::span<const std::uint8_t> trackingActive) const;
        [[nodiscard]] bool loadIsImbalanced(double threshold);
        /** @brief Whether the deviation exceeds the threshold on any rank. Collective. */
        [[nodiscard]] bool anyRankExceeds(double deviation, double threshold);
        [[nodiscard]] FieldStorage::Extents targetExtents(
                const DirichletPlaneConfig& dirichletPlane) const;
        void extendImageBounds(
                CartesianBounds& bounds, const DirichletPlaneConfig& dirichletPlane) const;
        void expandBounds(
                CartesianBounds& bounds, bool applyEmissionStretch, double emittedFraction,
                std::size_t longitudinalExtent) const;
        void updateMeshGeometry(const CartesianBounds& bounds, FieldStorage& fieldStorage) const;
        [[nodiscard]] bool isRedistributionDue(std::size_t step) const;
        [[nodiscard]] bool redistributionAllowed(const SpaceChargeSolveContext& context) const;
        /** @brief Reallocate the fields and migrate the primary container after an ORB call. */
        void applyRepartition(FieldStorage& fieldStorage);
        /** @brief ORB on the primary container's particle count. */
        [[nodiscard]] bool redistribute(
                const SpaceChargeSolveContext& context, FieldStorage& fieldStorage);
        /**
         * @brief Deposit the mixed particle-count and pair-count weight into the ORB field.
         * @return This rank's share of the total weight. Collective.
         */
        [[nodiscard]] double depositPairWeight(FieldStorage& fieldStorage);
        /**
         * @brief ORB on a weight that mixes particle count and P3M pair count.
         *
         * Runs only when a rank's weight share deviates from 1/ranks by more than the
         * load-balancing threshold. If the weighted cut leaves a rank region too small for the
         * requested P3M cutoff, while the previous layout fitted it better, the layout is restored
         * and redistribute() cuts by particle count instead.
         */
        [[nodiscard]] bool redistributeByPairWeight(
                const SpaceChargeSolveContext& context, FieldStorage& fieldStorage);
        [[nodiscard]] ParticleContainer& primary() const { return *particles_m.front(); }

        CartesianPIC3DConfig config_m;
        std::vector<ParticleContainer*> particles_m;
        Orb orb_m;
        std::vector<int> rankFlags_m;
        bool poissonRebuildRequired_m     = false;
        bool cutoffReductionReported_m    = false;
        bool pairWeightFallbackReported_m = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_DOMAIN_UPDATER_H
