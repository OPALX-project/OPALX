/**
 * @file PicDomainManager.cpp
 * @brief Implements Cartesian PIC domain updates and solver-owned redistribution.
 */

#include "SpaceCharge/Pic/PicDomainManager.h"

#include "SpaceCharge/Ippl/IpplPoissonAdapter.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace opalx::spacecharge {

    namespace {

        using Workspace = PicWorkspace<double, 3>;

        IpplPoissonFields backendFields(Workspace& workspace) {
            return {&workspace.chargeDensity(), &workspace.electricField(), &workspace.potential()};
        }

    }  // namespace

    PicDomainManager::PicDomainManager(Pic3DConfig config)
        : config_m(std::move(config)),
          rankFlags_m(static_cast<std::size_t>(ippl::Comm->size()), 0) {}

    void PicDomainManager::update(
            PicDomainFrame frame, SolveContext& context, Workspace& workspace,
            PicParticleDomainAdapter& particles, IpplPoissonAdapter& backend,
            SelfFieldDiagnostics& diagnostics) {
        const bool beamFrame = frame == PicDomainFrame::Beam;
        auto domainEvent     = diagnostics.scopedEvent(
                SelfFieldEventKind::DomainUpdate,
                beamFrame ? "beam-frame-mesh" : "reference-frame-mesh");

        const std::size_t step = context.stepState().step;
        bool refreshBackend    = workspace.rebuildGlobalLayoutInPlace(
                targetExtents(step), config_m.layoutRebuildParallelDimensions());

        PicDomainBounds bounds = particles.computeBounds();
        extendImageBounds(bounds, step);
        expandBounds(
                bounds, beamFrame && context.stepState().emissionActive,
                context.stepState().emittedFraction, workspace.layoutExtents()[2]);
        updateGeometry(bounds, workspace);

        // A mesh change invalidates native particle-layout views even when the global field
        // decomposition is unchanged. Migrate first and reacquire views only in later kernels.
        particles.updateLayoutsAndMigrate(
                workspace.layout(), workspace.mesh(), context.particles());
        particles.updateMoments();

        bool redistributed = false;
        if (beamFrame && redistributionDue(step)
            && particles.loadIsImbalanced(config_m.loadBalancingThreshold(), rankFlags_m)) {
            redistributed  = redistribute(context, workspace, particles, diagnostics);
            refreshBackend = refreshBackend || redistributed;
        }

        // Layout-owned fields and particle layouts now agree. Rebind backend RHS before LHS only
        // after the last possible layout change; IpplPoissonAdapter owns that exact ordering.
        if (refreshBackend) {
            backend.refresh(backendFields(workspace));
        }

        if (redistributed) {
            ippl::Comm->barrier();
            diagnostics.completeRedistribution();
        }
    }

    bool PicDomainManager::imageChargeActive(std::size_t step) const {
        const CorrectionConfig& correction = config_m.correction();
        return correction.kind() == CorrectionKind::ImageCharge
               && (correction.maximumSteps() == 0 || step < correction.maximumSteps());
    }

    Workspace::Extents PicDomainManager::targetExtents(std::size_t step) const {
        Workspace::Extents extents = config_m.meshSize();
        if (imageChargeActive(step)) {
            if (extents[2] > std::numeric_limits<std::size_t>::max() / 2) {
                throw OpalException(
                        "PicDomainManager::targetExtents",
                        "The image-charge longitudinal mesh extent overflows size_t.");
            }
            extents[2] *= 2;
        }
        return extents;
    }

    void PicDomainManager::extendImageBounds(PicDomainBounds& bounds, std::size_t step) const {
        if (!imageChargeActive(step)) {
            return;
        }

        const double planeZ       = config_m.correction().planeZ();
        const double mirroredMinZ = 2.0 * planeZ - bounds.upper[2];
        const double mirroredMaxZ = 2.0 * planeZ - bounds.lower[2];
        bounds.lower[2]           = std::min(bounds.lower[2], mirroredMinZ);
        bounds.upper[2]           = std::max(bounds.upper[2], mirroredMaxZ);
    }

    void PicDomainManager::expandBounds(
            PicDomainBounds& bounds, bool applyEmissionStretch, double emittedFraction,
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

    void PicDomainManager::updateGeometry(
            const PicDomainBounds& bounds, Workspace& workspace) const {
        const Workspace::Extents extents = workspace.layoutExtents();
        const Workspace::Vector span     = bounds.upper - bounds.lower;
        Workspace::Vector spacing(0.0);
        Workspace::Vector origin = bounds.lower;
        for (unsigned dimension = 0; dimension < 3; ++dimension) {
            spacing[dimension] =
                    extents[dimension] <= 1
                            ? span[dimension]
                            : span[dimension] / static_cast<double>(extents[dimension] - 1);
            origin[dimension] = bounds.lower[dimension] - 0.5 * spacing[dimension];
        }
        workspace.setGeometry(bounds.lower, bounds.upper, spacing, origin);
    }

    bool PicDomainManager::redistributionDue(std::size_t step) const {
        const std::size_t frequency = config_m.repartitionFrequency();
        return frequency > 0 && step % frequency + 1 == frequency;
    }

    bool PicDomainManager::redistribute(
            SolveContext& context, Workspace& workspace, PicParticleDomainAdapter& particles,
            SelfFieldDiagnostics& diagnostics) {
        Inform m("PicDomainManager::redistribute");
        const int ranks = context.stepState().communicator.size;
        if (ranks < 2 || particles.primaryTotalCount() == 0) {
            return false;
        }
        if ((ranks & (ranks - 1)) != 0) {
            m << level1 << "Skipping ORB load balancing because ORB requires a power-of-two MPI "
              << "rank count; current rank count is " << ranks << "." << endl;
            return false;
        }
        if (particles.redistributionBlocked(context.particles())) {
            m << level2
              << "Skipping ORB load balancing because a non-primary particle container is "
                 "tracking-active or non-empty."
              << endl;
            return false;
        }

        auto redistributionEvent =
                diagnostics.scopedEvent(SelfFieldEventKind::DomainUpdate, "redistribution");
        const std::size_t localBefore = particles.primaryLocalCount();
        orb_m                         = Orb();
        orb_m.initialize(workspace.layout(), workspace.mesh(), workspace.chargeDensity());
        const bool succeeded =
                orb_m.binaryRepartition(particles.primaryPositions(), workspace.layout(), false);
        if (!succeeded) {
            m << level2 << "ORB load balancing failed; keeping the previous layout." << endl;
            return false;
        }

        // ORB mutates the stable FieldLayout object. Refresh every field before migrating every
        // particle container into that decomposition. No pre-migration Kokkos view is retained.
        workspace.updateFieldLayoutsAfterLayoutChange();
        particles.updateLayoutsAndMigrate(
                workspace.layout(), workspace.mesh(), context.particles());

        m << level2 << "ORB load balancing done. Rank 0: " << localBefore << " -> "
          << particles.primaryLocalCount() << " particles in primary container." << endl;
        return true;
    }

}  // namespace opalx::spacecharge
