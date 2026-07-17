/**
 * @file LegacyPic3DAlgorithm.cpp
 * @brief Implements the temporary PartBunch self-field bridge.
 */

#include "SpaceCharge/LegacyPic3DAlgorithm.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "PartBunch/BinnedFieldSolver.h"
#include "PartBunch/PartBunch.h"

#include <exception>
#include <utility>

namespace opalx::spacecharge {

    LegacyPic3DAlgorithm::LegacyPic3DAlgorithm(
            PartBunch_t& bunch, Pic3DConfig config,
            std::shared_ptr<PicWorkspace<double, 3>> workspace)
        : workspace_m(std::move(workspace)),
          particleDomain_m(bunch.getParticleContainers()),
          domainManager_m(std::move(config)),
          bunch_m(&bunch) {
        if (workspace_m == nullptr) {
            throw OpalException(
                    "LegacyPic3DAlgorithm::LegacyPic3DAlgorithm", "PIC workspace is null.");
        }
    }

    SolverCapabilities LegacyPic3DAlgorithm::capabilities() const {
        SolverCapabilities result;
        result.algorithm                  = SelfFieldAlgorithmKind::Pic3D;
        result.supportsBinning            = true;
        result.supportsImageCharge        = true;
        result.supportsShiftedGreen       = true;
        result.supportsRedistribution     = true;
        result.supportsMultipleContainers = false;
        result.supportsPotentialOutput    = true;

        result.requiredReadableAttributes =
                ParticleAttribute::Position | ParticleAttribute::Momentum
                | ParticleAttribute::Charge | ParticleAttribute::Mass | ParticleAttribute::TimeStep
                | ParticleAttribute::InvalidMask | ParticleAttribute::Bin;
        result.requiredWritableAttributes =
                ParticleAttribute::Position | ParticleAttribute::TimeStep
                | ParticleAttribute::ElectricField | ParticleAttribute::MagneticField
                | ParticleAttribute::Bin;
        return result;
    }

    void LegacyPic3DAlgorithm::execute(SolveContext& context, SelfFieldDiagnostics& diagnostics) {
        const FrameState& frames = context.stepState().frames;
        if (!frames.trackerToSolve.has_value() || !frames.solveToTracker.has_value()) {
            throw OpalException(
                    "LegacyPic3DAlgorithm::execute",
                    "The 3D PIC algorithm requires both per-call coordinate transforms.");
        }

        const auto& referenceToBeam = frames.trackerToSolve->native<CoordinateSystemTrafo>();
        const auto& beamToReference = frames.solveToTracker->native<CoordinateSystemTrafo>();
        auto primary                = bunch_m->getParticleContainer();
        auto& backend               = bunch_m->getFieldSolver()->getBackendAdapter();
        context.particles().updateGeneration(particleDomain_m.generation());

        bool positionsInBeam = false;
        bool electricInBeam  = false;
        bool magneticInBeam  = false;
        try {
            positionsInBeam = true;
            referenceToBeam.transformBunchTo(primary->R.getView(), primary->getLocalNum());

            domainManager_m.update(
                    PicDomainFrame::Beam, context, *workspace_m, particleDomain_m, backend,
                    diagnostics);

            // Domain migration may replace particle storage. The legacy solver reacquires every
            // native Kokkos view internally after this point.
            electricInBeam = true;
            magneticInBeam = true;
            bunch_m->computeSelfFields(diagnostics);

            const std::size_t localCount = primary->getLocalNum();
            beamToReference.transformBunchTo(primary->R.getView(), localCount);
            positionsInBeam = false;
            beamToReference.rotateBunchTo(primary->E.getView(), localCount);
            electricInBeam = false;
            beamToReference.rotateBunchTo(primary->B.getView(), localCount);
            magneticInBeam = false;

            domainManager_m.update(
                    PicDomainFrame::Reference, context, *workspace_m, particleDomain_m, backend,
                    diagnostics);
        } catch (...) {
            const std::exception_ptr originalException = std::current_exception();
            // Keep borrowed particle state in the tracker frame even when a host-side domain or
            // backend operation fails. Reacquire all views because migration may have occurred;
            // cleanup failure must not replace the original exception.
            try {
                const std::size_t localCount = primary->getLocalNum();
                if (positionsInBeam) {
                    beamToReference.transformBunchTo(primary->R.getView(), localCount);
                    positionsInBeam = false;
                }
                if (electricInBeam) {
                    beamToReference.rotateBunchTo(primary->E.getView(), localCount);
                }
                if (magneticInBeam) {
                    beamToReference.rotateBunchTo(primary->B.getView(), localCount);
                }
                if (!positionsInBeam) {
                    domainManager_m.update(
                            PicDomainFrame::Reference, context, *workspace_m, particleDomain_m,
                            backend, diagnostics);
                }
            } catch (...) {
            }
            std::rethrow_exception(originalException);
        }
    }

}  // namespace opalx::spacecharge
