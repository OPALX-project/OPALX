/**
 * @file LegacyPic3DAlgorithm.cpp
 * @brief Implements the temporary PartBunch self-field bridge.
 */

#include "SpaceCharge/LegacyPic3DAlgorithm.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "PartBunch/BinnedFieldSolver.h"
#include "PartBunch/PartBunch.h"
#include "SpaceCharge/Pic/IterationPlanFactory.h"
#include "Structure/DataSink.h"

#include <exception>
#include <optional>
#include <utility>
#include <vector>

namespace opalx::spacecharge {

    namespace {

        using ParticleContainer = ::ParticleContainer<double, 3>;

        ParticleContainer& requirePrimaryParticles(PartBunch_t& bunch) {
            const std::shared_ptr<ParticleContainer> primary = bunch.getParticleContainer();
            if (primary == nullptr) {
                throw OpalException(
                        "LegacyPic3DAlgorithm::LegacyPic3DAlgorithm",
                        "The primary particle container is null.");
            }
            return *primary;
        }

        /** @brief Synchronously forwards requested host bin snapshots to the run DataSink. */
        class DataSinkBinConfigurationObserver final : public BinConfigurationObserver {
        public:
            DataSinkBinConfigurationObserver(
                    DataSink* dataSink, long long step, double time, const BinningConfig& config)
                : dataSink_m(dataSink), step_m(step), time_m(time), config_m(config) {}

            [[nodiscard]] bool wants(BinConfigurationPoint) const override {
                return dataSink_m != nullptr && !config_m.dumpFile().empty()
                       && config_m.dumpFrequency() > 0
                       && step_m % static_cast<long long>(config_m.dumpFrequency()) == 0;
            }

            void record(BinConfigurationPoint point, const BinConfigurationSnapshot& snapshot)
                    override {
                if (!wants(point)) {
                    return;
                }

                const std::vector<std::size_t> particleCounts(
                        snapshot.particleCounts.begin(), snapshot.particleCounts.end());
                dataSink_m->dumpBinConfig(
                        step_m, time_m, point == BinConfigurationPoint::BeforeMerge, particleCounts,
                        snapshot.widths, snapshot.lowerBound, config_m.dumpFile());
            }

        private:
            DataSink* dataSink_m = nullptr;
            long long step_m     = 0;
            double time_m        = 0.0;
            const BinningConfig& config_m;
        };

    }  // namespace

    LegacyPic3DAlgorithm::LegacyPic3DAlgorithm(
            PartBunch_t& bunch, Pic3DConfig config,
            std::shared_ptr<PicWorkspace<double, 3>> workspace)
        : workspace_m(std::move(workspace)),
          binningConfig_m(config.binning()),
          iterationPlan_m(
                  IterationPlanFactory<double, 3>::create(
                          binningConfig_m, requirePrimaryParticles(bunch))),
          particleDomain_m(bunch.getParticleContainers()),
          domainManager_m(std::move(config)),
          bunch_m(&bunch) {
        if (workspace_m == nullptr) {
            throw OpalException(
                    "LegacyPic3DAlgorithm::LegacyPic3DAlgorithm", "PIC workspace is null.");
        }

        BinnedFieldSolver<double, 3>* fieldSolver = bunch_m->getFieldSolver();
        if (fieldSolver == nullptr) {
            throw OpalException(
                    "LegacyPic3DAlgorithm::LegacyPic3DAlgorithm", "The PIC field solver is null.");
        }
        fieldSolver->configureIterationMetadata(
                iterationPlan_m->kind() == IterationKind::Binning,
                iterationPlan_m->maximumBinCount());
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
            std::optional<DataSinkBinConfigurationObserver> binObserver;
            if (binningConfig_m.has_value()) {
                binObserver.emplace(
                        bunch_m->getDataSink(), bunch_m->getGlobalTrackStep(), bunch_m->getT(),
                        *binningConfig_m);
            }

            electricInBeam = true;
            magneticInBeam = true;
            bunch_m->computeSelfFields(
                    *iterationPlan_m, context.particles().generation(),
                    binObserver.has_value() ? &*binObserver : nullptr, diagnostics);

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
