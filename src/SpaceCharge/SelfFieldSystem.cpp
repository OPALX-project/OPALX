/**
 * @file SelfFieldSystem.cpp
 * @brief Implements common self-field request validation and dispatch.
 */

#include "SpaceCharge/SelfFieldSystem.h"

#include "Utilities/OpalException.h"

#include <string>
#include <utility>

namespace opalx::spacecharge {

    namespace {

        SelfFieldDiagnosticSchedule makeDiagnosticSchedule(const SelfFieldConfig& config) {
            const Pic3DConfig& picConfig = config.get<Pic3DConfig>();
            SelfFieldDiagnosticSchedule schedule;
            schedule.planeDumpFrequency = picConfig.correction().planeDumpFrequency();
            if (picConfig.binning().has_value()) {
                schedule.binTableFrequency = picConfig.binning()->tablePrintFrequency();
            }
            return schedule;
        }

    }  // namespace

    SelfFieldSystem::SelfFieldSystem(
            SelfFieldConfig config, std::unique_ptr<SelfFieldAlgorithm> algorithm)
        : config_m(std::move(config)),
          algorithm_m(std::move(algorithm)),
          diagnostics_m(makeDiagnosticSchedule(config_m)) {
        if (algorithm_m == nullptr) {
            throw OpalException(
                    "SelfFieldSystem::SelfFieldSystem", "The self-field algorithm is null.");
        }
        capabilities_m = algorithm_m->capabilities();
        validateConfiguration();
    }

    void SelfFieldSystem::solve(SolveContext& context) {
        validateContext(context);
        auto solveEvent = diagnostics_m.scopedEvent(SelfFieldEventKind::Solve, "self-field");
        algorithm_m->execute(context, diagnostics_m);
    }

    RequestedPhysics SelfFieldSystem::requestedPhysicsForStep(std::size_t step) const {
        const Pic3DConfig& picConfig       = config_m.get<Pic3DConfig>();
        const CorrectionConfig& correction = picConfig.correction();
        const bool correctionActive =
                correction.enabled()
                && (correction.maximumSteps() == 0 || step < correction.maximumSteps());

        RequestedPhysics requested;
        requested.useBinning = picConfig.binning().has_value();
        if (correctionActive) {
            requested.correction     = {correction.kind(), correction.planeZ()};
            requested.writePotential = correction.kind() == CorrectionKind::ImageCharge
                                       && correction.planeDumpFrequency() != 0;
        }
        return requested;
    }

    CorrectionRequest SelfFieldSystem::configuredCorrection() const {
        const CorrectionConfig& correction = config_m.get<Pic3DConfig>().correction();
        return {correction.kind(), correction.planeZ()};
    }

    int SelfFieldSystem::reportedBinCount() const {
        const auto& binning = config_m.get<Pic3DConfig>().binning();
        if (!binning.has_value() || !diagnostics_m.hasCurrentBinCount()
            || diagnostics_m.currentBinCount() == binning->maximumBins()) {
            return 1;
        }
        return static_cast<int>(diagnostics_m.currentBinCount());
    }

    void SelfFieldSystem::validateConfiguration() const {
        if (capabilities_m.algorithm != config_m.algorithmKind()) {
            throw OpalException(
                    "SelfFieldSystem::SelfFieldSystem",
                    "The selected algorithm does not match the self-field configuration.");
        }

        const Pic3DConfig& picConfig = config_m.get<Pic3DConfig>();
        if (picConfig.binning().has_value() && !capabilities_m.supportsBinning) {
            throw OpalException(
                    "SelfFieldSystem::SelfFieldSystem",
                    "The selected self-field algorithm does not support binning.");
        }
        if (!capabilities_m.supports(picConfig.correction().kind())) {
            throw OpalException(
                    "SelfFieldSystem::SelfFieldSystem",
                    "The selected self-field algorithm does not support the configured "
                    "source-plane correction.");
        }
        if (picConfig.correction().planeDumpFrequency() != 0
            && !capabilities_m.supportsPotentialOutput) {
            throw OpalException(
                    "SelfFieldSystem::SelfFieldSystem",
                    "The selected self-field algorithm cannot write the configured potential "
                    "diagnostic.");
        }
        if (picConfig.repartitionFrequency() != 0 && !capabilities_m.supportsRedistribution) {
            throw OpalException(
                    "SelfFieldSystem::SelfFieldSystem",
                    "The selected self-field algorithm does not support particle "
                    "redistribution.");
        }
    }

    void SelfFieldSystem::validateContext(const SolveContext& context) const {
        const ParticleSetView& particles = context.particles();
        if (!capabilities_m.supportsMultipleContainers && particles.activeContainerCount() > 1) {
            throw OpalException(
                    "SelfFieldSystem::solve",
                    "The selected self-field algorithm supports only one active particle "
                    "container.");
        }
        if (!particles.activeContainersProvide(capabilities_m.requiredReadableAttributes)) {
            throw OpalException(
                    "SelfFieldSystem::solve",
                    "An active particle container is missing a readable attribute required by "
                    "the self-field algorithm (required bits: "
                            + std::to_string(capabilities_m.requiredReadableAttributes.bits())
                            + ").");
        }
        if (!particles.activeContainersProvideWritable(capabilities_m.requiredWritableAttributes)) {
            throw OpalException(
                    "SelfFieldSystem::solve",
                    "An active particle container is missing a writable attribute required by "
                    "the self-field algorithm (required bits: "
                            + std::to_string(capabilities_m.requiredWritableAttributes.bits())
                            + ").");
        }

        const RequestedPhysics& requested = context.requestedPhysics();
        if (requested.useBinning && !capabilities_m.supportsBinning) {
            throw OpalException(
                    "SelfFieldSystem::solve",
                    "The solve requests binning, but the selected algorithm does not support it.");
        }
        if (requested.writePotential && !capabilities_m.supportsPotentialOutput) {
            throw OpalException(
                    "SelfFieldSystem::solve",
                    "The solve requests potential output, but the selected algorithm does not "
                    "support it.");
        }
        if (!capabilities_m.supports(requested.correction.kind)) {
            throw OpalException(
                    "SelfFieldSystem::solve",
                    "The solve requests a correction that the selected algorithm does not "
                    "support.");
        }

        const CorrectionKind configuredCorrection = config_m.get<Pic3DConfig>().correction().kind();
        if (requested.correction.kind != CorrectionKind::None
            && requested.correction.kind != configuredCorrection) {
            throw OpalException(
                    "SelfFieldSystem::solve",
                    "The solve requests a correction different from the immutable run "
                    "configuration.");
        }
    }

}  // namespace opalx::spacecharge
