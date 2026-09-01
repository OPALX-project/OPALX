/**
 * @file SelfFieldSystem.cpp
 * @brief Implements common self-field request validation and dispatch.
 */

#include "SpaceCharge/SelfFieldSystem.h"

#include "Utilities/OpalException.h"

#include <string>
#include <type_traits>
#include <utility>

namespace opalx::spacecharge {

    namespace {

        SelfFieldDiagnosticSchedule makeDiagnosticSchedule(const SelfFieldConfig& config) {
            return std::visit(
                    [](const auto& algorithmConfig) {
                        using Config = std::decay_t<decltype(algorithmConfig)>;
                        SelfFieldDiagnosticSchedule schedule;
                        if constexpr (std::is_same_v<Config, Pic3DConfig>) {
                            schedule.planeDumpFrequency =
                                    algorithmConfig.correction().planeDumpFrequency();
                            if (algorithmConfig.binning().has_value()) {
                                schedule.binTableFrequency =
                                        algorithmConfig.binning()->tablePrintFrequency();
                            }
                        }
                        return schedule;
                    },
                    config.algorithmConfig());
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
        context.particles().applySelection(capabilities_m.particleSelection);
        validateContext(context);
        auto solveEvent = diagnostics_m.scopedEvent(SelfFieldEventKind::Solve, "self-field");
        algorithm_m->execute(context, diagnostics_m);
    }

    RequestedPhysics SelfFieldSystem::requestedPhysicsForStep(std::size_t step) const {
        return std::visit(
                [step](const auto& algorithmConfig) {
                    using Config = std::decay_t<decltype(algorithmConfig)>;
                    RequestedPhysics requested;
                    if constexpr (std::is_same_v<Config, Pic3DConfig>) {
                        const CorrectionConfig& correction = algorithmConfig.correction();
                        const bool correctionActive        = correction.enabled()
                                                      && (correction.maximumSteps() == 0
                                                          || step < correction.maximumSteps());
                        requested.useBinning = algorithmConfig.binning().has_value();
                        if (correctionActive) {
                            requested.correction = {correction.kind(), correction.planeZ()};
                            requested.writePotential =
                                    correction.kind() == CorrectionKind::ImageCharge
                                    && correction.planeDumpFrequency() != 0;
                        }
                    }
                    return requested;
                },
                config_m.algorithmConfig());
    }

    CorrectionRequest SelfFieldSystem::configuredCorrection() const {
        return std::visit(
                [](const auto& algorithmConfig) {
                    using Config = std::decay_t<decltype(algorithmConfig)>;
                    if constexpr (std::is_same_v<Config, Pic3DConfig>) {
                        const CorrectionConfig& correction = algorithmConfig.correction();
                        return CorrectionRequest{correction.kind(), correction.planeZ()};
                    }
                    return CorrectionRequest{};
                },
                config_m.algorithmConfig());
    }

    int SelfFieldSystem::reportedBinCount() const {
        return std::visit(
                [this](const auto& algorithmConfig) {
                    using Config = std::decay_t<decltype(algorithmConfig)>;
                    if constexpr (std::is_same_v<Config, Pic3DConfig>) {
                        const auto& binning = algorithmConfig.binning();
                        if (binning.has_value() && diagnostics_m.hasCurrentBinCount()
                            && diagnostics_m.currentBinCount() != binning->maximumBins()) {
                            return static_cast<int>(diagnostics_m.currentBinCount());
                        }
                    }
                    return 1;
                },
                config_m.algorithmConfig());
    }

    void SelfFieldSystem::validateConfiguration() const {
        if (capabilities_m.algorithm != config_m.algorithmKind()) {
            throw OpalException(
                    "SelfFieldSystem::SelfFieldSystem",
                    "The selected algorithm does not match the self-field configuration.");
        }

        std::visit(
                [this](const auto& algorithmConfig) {
                    using Config = std::decay_t<decltype(algorithmConfig)>;
                    if constexpr (std::is_same_v<Config, Pic3DConfig>) {
                        if (algorithmConfig.binning().has_value()
                            && !capabilities_m.supportsBinning) {
                            throw OpalException(
                                    "SelfFieldSystem::SelfFieldSystem",
                                    "The selected self-field algorithm does not support binning.");
                        }
                        if (!capabilities_m.supports(algorithmConfig.correction().kind())) {
                            throw OpalException(
                                    "SelfFieldSystem::SelfFieldSystem",
                                    "The selected self-field algorithm does not support the "
                                    "configured source-plane correction.");
                        }
                        if (algorithmConfig.correction().planeDumpFrequency() != 0
                            && !capabilities_m.supportsPotentialOutput) {
                            throw OpalException(
                                    "SelfFieldSystem::SelfFieldSystem",
                                    "The selected self-field algorithm cannot write the "
                                    "configured potential diagnostic.");
                        }
                        if (algorithmConfig.repartitionFrequency() != 0
                            && !capabilities_m.supportsRedistribution) {
                            throw OpalException(
                                    "SelfFieldSystem::SelfFieldSystem",
                                    "The selected self-field algorithm does not support particle "
                                    "redistribution.");
                        }
                    }
                },
                config_m.algorithmConfig());
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

        const CorrectionKind configuredKind = this->configuredCorrection().kind;
        if (requested.correction.kind != CorrectionKind::None
            && requested.correction.kind != configuredKind) {
            throw OpalException(
                    "SelfFieldSystem::solve",
                    "The solve requests a correction different from the immutable run "
                    "configuration.");
        }
    }

}  // namespace opalx::spacecharge
