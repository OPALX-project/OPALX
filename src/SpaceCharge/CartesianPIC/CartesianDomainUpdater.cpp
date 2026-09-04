/**
 * @file CartesianDomainUpdater.cpp
 * @brief Implements Cartesian PIC domain updates and solver-owned redistribution.
 */

#include "SpaceCharge/CartesianPIC/CartesianDomainUpdater.h"

#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/Poisson/PoissonSolver.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace opalx::spacecharge {

    namespace {

        using FieldStorage = CartesianPICFieldStorage<double, 3>;

        PoissonFieldBinding makePoissonFieldBinding(FieldStorage& fieldStorage) {
            return {&fieldStorage.chargeDensity(), &fieldStorage.electricField()};
        }

    }  // namespace

    CartesianDomainUpdater::CartesianDomainUpdater(
            CartesianPICConfig config, std::shared_ptr<const BunchStateHandler> bunchState)
        : config_m(std::move(config)),
          bunchState_m(std::move(bunchState)),
          rankFlags_m(static_cast<std::size_t>(ippl::Comm->size()), 0) {
        if (bunchState_m == nullptr) {
            throw OpalException(
                    "CartesianDomainUpdater::CartesianDomainUpdater",
                    "The bunch state handler is null.");
        }
    }

    void CartesianDomainUpdater::updateForSolve(
            DomainCoordinateFrame frame, SpaceChargeSolveContext& context,
            FieldStorage& fieldStorage, ParticleDomainOperations& particles,
            PoissonSolver& poissonSolver, SpaceChargeDiagnostics& diagnostics) {
        const bool beamFrame    = frame == DomainCoordinateFrame::Beam;
        const auto& fixedDomain = bunchState_m->fixedCartesianDomain();

        const std::size_t step                         = context.stepState().step;
        const SpaceChargeCorrectionRequest& correction = context.request().correction;
        const FieldStorage::Extents extents =
                fixedDomain.has_value() ? config_m.meshSize() : targetExtents(correction);
        const bool repeatFieldLayoutRefresh = poissonRebuildRequired_m;
        if (fieldStorage.layoutExtents() != extents) {
            // Backend-private fields and plans must be rebuilt after changing layout-owned fields.
            poissonRebuildRequired_m = true;
        }
        const std::array<bool, 3>& decomposition =
                fixedDomain.has_value() ? config_m.parallelDimensions()
                                        : config_m.layoutRebuildParallelDimensions();
        const bool layoutChanged =
                fieldStorage.domain().rebuildGlobalLayoutInPlace(extents, decomposition);
        if (layoutChanged || repeatFieldLayoutRefresh) {
            fieldStorage.updateFieldLayoutsAfterLayoutChange();
        }

        CartesianBounds bounds;
        if (fixedDomain.has_value()) {
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                bounds.lower[dimension] = fixedDomain->lower[dimension];
                bounds.upper[dimension] = fixedDomain->upper[dimension];
            }
        } else {
            bounds = particles.computeBounds();
            extendImageBounds(bounds, correction);
            expandBounds(
                    bounds, beamFrame && context.stepState().emissionActive,
                    context.stepState().emittedFraction, fieldStorage.layoutExtents()[2]);
        }
        updateMeshGeometry(bounds, fieldStorage);

        // A mesh change invalidates native particle-layout views even when the global field
        // decomposition is unchanged. Migrate first and reacquire views only in later kernels.
        if (fixedDomain.has_value()) {
            // Secondary containers remain in tracker coordinates and must not be migrated through
            // a mesh whose bounds are expressed in the primary solve frame.
            particles.updatePrimaryLayoutAndMigrate(fieldStorage.layout(), fieldStorage.mesh());
            particles.updatePrimaryMoments();
        } else {
            particles.updateLayoutsAndMigrate(fieldStorage.layout(), fieldStorage.mesh());
            particles.updateMoments();
        }

        bool redistributed = false;
        // Keep the expensive collective imbalance check behind the configured step cadence. This
        // preserves the old short-circuit behavior when no repartition is due.
        if (!fixedDomain.has_value() && beamFrame && isRedistributionDue(step)
            && particles.loadIsImbalanced(config_m.loadBalancingThreshold(), rankFlags_m)) {
            redistributed = redistribute(context, fieldStorage, particles);
        }

        // Layout-owned fields and particle layouts now agree. Rebuild and bind the Poisson RHS
        // before its LHS only after the last possible layout change.
        if (poissonRebuildRequired_m) {
            poissonSolver.rebuildAfterLayoutChange(makePoissonFieldBinding(fieldStorage));
            poissonRebuildRequired_m = false;
        }

        if (redistributed) {
            ippl::Comm->barrier();
            diagnostics.recordRedistribution();
        }
    }

    FieldStorage::Extents CartesianDomainUpdater::targetExtents(
            const SpaceChargeCorrectionRequest& correction) const {
        FieldStorage::Extents extents = config_m.meshSize();
        if (correction.kind == SpaceChargeCorrectionType::ImageCharge) {
            if (extents[2] > std::numeric_limits<std::size_t>::max() / 2) {
                throw OpalException(
                        "CartesianDomainUpdater::targetExtents",
                        "The image-charge longitudinal mesh extent overflows size_t.");
            }
            extents[2] *= 2;
        }
        return extents;
    }

    void CartesianDomainUpdater::extendImageBounds(
            CartesianBounds& bounds, const SpaceChargeCorrectionRequest& correction) const {
        if (correction.kind != SpaceChargeCorrectionType::ImageCharge) {
            return;
        }

        const double planeZ       = correction.planeZ;
        const double mirroredMinZ = 2.0 * planeZ - bounds.upper[2];
        const double mirroredMaxZ = 2.0 * planeZ - bounds.lower[2];
        bounds.lower[2]           = std::min(bounds.lower[2], mirroredMinZ);
        bounds.upper[2]           = std::max(bounds.upper[2], mirroredMaxZ);
    }

    void CartesianDomainUpdater::expandBounds(
            CartesianBounds& bounds, bool applyEmissionStretch, double emittedFraction,
            std::size_t longitudinalExtent) const {
        ippl::Vector<double, 3> span = bounds.upper - bounds.lower;
        for (unsigned dimension = 0; dimension < 3; ++dimension) {
            if (span[dimension] < 1.0e-6) {
                span[dimension] = 1.0e-6;
            }
        }

        const double relativeIncrement = config_m.boundingBoxIncreasePercent() / 100.0;
        bounds.lower                   = bounds.lower - span * relativeIncrement;
        bounds.upper                   = bounds.upper + span * relativeIncrement;

        if (!applyEmissionStretch || longitudinalExtent <= 1) {
            return;
        }

        // Match old OPAL while emission is running by extending the beam-frame domain over the
        // complete source pulse. At emitted fraction f, the physical z length is scaled by 1/f.
        // The one-cell lower bound prevents a nearly empty pulse from requesting an unbounded
        // domain. Without this stretch, early charge is compressed and its self-field is too large.
        const double percent =
                std::max(1.0 / static_cast<double>(longitudinalExtent - 1), emittedFraction);
        const double originalLength =
                std::abs(bounds.upper[2] - bounds.lower[2]) / (1.0 + 2.0 * relativeIncrement);
        if (percent >= 1.0 || percent <= 0.0 || originalLength <= 0.0) {
            return;
        }

        bounds.upper[2] -= relativeIncrement * originalLength;
        bounds.lower[2] = bounds.upper[2] - originalLength / percent;

        const double stretchedLength = originalLength / percent;
        bounds.upper[2] += relativeIncrement * stretchedLength;
        bounds.lower[2] -= relativeIncrement * stretchedLength;
    }

    void CartesianDomainUpdater::updateMeshGeometry(
            const CartesianBounds& bounds, FieldStorage& fieldStorage) const {
        const FieldStorage::Extents extents = fieldStorage.layoutExtents();
        const FieldStorage::Vector span     = bounds.upper - bounds.lower;
        FieldStorage::Vector spacing(0.0);
        FieldStorage::Vector origin = bounds.lower;
        for (unsigned dimension = 0; dimension < 3; ++dimension) {
            spacing[dimension] =
                    extents[dimension] <= 1
                            ? span[dimension]
                            : span[dimension] / static_cast<double>(extents[dimension] - 1);
            origin[dimension] = bounds.lower[dimension] - 0.5 * spacing[dimension];
        }
        fieldStorage.domain().setGeometry(bounds.lower, bounds.upper, spacing, origin);
    }

    bool CartesianDomainUpdater::isRedistributionDue(std::size_t step) const {
        const std::size_t frequency = config_m.repartitionFrequency();
        return frequency > 0 && step % frequency + 1 == frequency;
    }

    bool CartesianDomainUpdater::redistribute(
            SpaceChargeSolveContext& context, FieldStorage& fieldStorage,
            ParticleDomainOperations& particles) {
        Inform m("CartesianDomainUpdater::redistribute");
        const int ranks = context.stepState().mpiSize;
        if (ranks < 2 || particles.primaryTotalCount() == 0) {
            return false;
        }
        if ((ranks & (ranks - 1)) != 0) {
            m << level1 << "Skipping ORB load balancing because ORB requires a power-of-two MPI "
              << "rank count; current rank count is " << ranks << "." << endl;
            return false;
        }
        if (particles.isRedistributionBlocked(context.particles())) {
            // Cartesian PIC solves and constructs ORB weights from the primary container only.
            // Repartitioning it while another active or populated container shares the old layout
            // would leave the bunch with inconsistent spatial ownership, so defer ORB entirely.
            m << level2
              << "Skipping ORB load balancing because a non-primary particle container is "
                 "tracking-active or non-empty."
              << endl;
            return false;
        }

        const std::size_t localBefore = particles.primaryLocalCount();
        orb_m                         = Orb();
        orb_m.initialize(fieldStorage.layout(), fieldStorage.mesh(), fieldStorage.chargeDensity());
        // ORB can mutate the layout even when it reports ordinary nonthrowing failure. Mark the
        // backend dirty before the call so fields and plans are refreshed in either result.
        poissonRebuildRequired_m = true;
        const bool succeeded =
                orb_m.binaryRepartition(particles.primaryPositions(), fieldStorage.layout(), false);
        if (!succeeded) {
            m << level2 << "ORB load balancing failed; keeping the previous layout." << endl;
            return false;
        }

        // ORB mutates the stable FieldLayout object. Refresh every field before migrating every
        // particle container into that decomposition. No pre-migration Kokkos view is retained.
        fieldStorage.updateFieldLayoutsAfterLayoutChange();
        particles.updateLayoutsAndMigrate(fieldStorage.layout(), fieldStorage.mesh());

        m << level2 << "ORB load balancing done. Rank 0: " << localBefore << " -> "
          << particles.primaryLocalCount() << " particles in primary container." << endl;
        return true;
    }

}  // namespace opalx::spacecharge
