/**
 * @file SelfFieldSystem.cpp
 * @brief Implements common self-field binding and request validation.
 */

#include "SpaceCharge/SelfFieldSystem.h"

#include "SpaceCharge/SelfFieldRequestPolicy.h"
#include "Utilities/OpalException.h"

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
            SelfFieldConfig config, std::unique_ptr<SelfFieldAlgorithm> algorithm,
            std::vector<ParticleFieldBinding3d> bindings, std::size_t primaryIndex)
        : config_m(std::move(config)),
          algorithm_m(std::move(algorithm)),
          bindings_m(std::move(bindings)),
          primaryIndex_m(primaryIndex),
          diagnostics_m(makeDiagnosticSchedule(config_m)) {
        if (algorithm_m == nullptr) {
            throw OpalException(
                    "SelfFieldSystem::SelfFieldSystem", "The self-field algorithm is null.");
        }
        if (bindings_m.empty() || primaryIndex_m >= bindings_m.size()) {
            throw OpalException(
                    "SelfFieldSystem::SelfFieldSystem", "The particle binding set is invalid.");
        }
        for (const ParticleFieldBinding3d& binding : bindings_m) {
            if (!binding.hasCompleteIdentity()) {
                throw OpalException(
                        "SelfFieldSystem::SelfFieldSystem",
                        "The particle binding set contains an incomplete binding.");
            }
        }
        capabilities_m = algorithm_m->capabilities();
        validateConfiguration();
    }

    void SelfFieldSystem::solve(SolveContext& context) {
        validateBindings(context.particles());
        context.particles().applySelection(capabilities_m.particleSelection);
        validateRequest(context);
        algorithm_m->execute(context, diagnostics_m);
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

    void SelfFieldSystem::validateBindings(const ParticleSetView& particles) const {
        const auto bindings = particles.bindings();
        if (bindings.size() != bindings_m.size() || particles.primaryIndex() != primaryIndex_m) {
            throw OpalException(
                    "SelfFieldSystem::solve",
                    "The solve context particle binding set does not match construction.");
        }
        for (std::size_t index = 0; index < bindings.size(); ++index) {
            if (!bindings[index].sameIdentity(bindings_m[index])) {
                throw OpalException(
                        "SelfFieldSystem::solve",
                        "The solve context particle or R/P/E/B binding does not match "
                        "construction at index "
                                + std::to_string(index) + ".");
            }
        }
    }

    void SelfFieldSystem::validateRequest(const SolveContext& context) const {
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
        const RequestedPhysics expected =
                SelfFieldRequestPolicy(config_m).forStep(context.stepState().step);
        const bool wrongPlane = requested.correction.kind != CorrectionKind::None
                                && requested.correction.planeZ != expected.correction.planeZ;
        if (requested.useBinning != expected.useBinning
            || requested.writePotential != expected.writePotential
            || requested.correction.kind != expected.correction.kind || wrongPlane) {
            throw OpalException(
                    "SelfFieldSystem::solve",
                    "The requested physics differs from the immutable step-resolved "
                    "configuration.");
        }
    }

}  // namespace opalx::spacecharge
