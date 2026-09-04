/**
 * @file SpaceChargeSolver.cpp
 * @brief Implements common space-charge binding and request validation.
 */

#include "SpaceCharge/SpaceChargeSolver.h"

#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/SpaceChargeRequestSchedule.h"
#include "Utilities/OpalException.h"

#include <type_traits>
#include <utility>

namespace opalx::spacecharge {
    namespace {

        SpaceChargeDiagnosticSchedule makeDiagnosticSchedule(const SpaceChargeConfig& config) {
            return std::visit(
                    [](const auto& algorithmConfig) {
                        using Config = std::decay_t<decltype(algorithmConfig)>;
                        SpaceChargeDiagnosticSchedule schedule;
                        if constexpr (std::is_same_v<Config, CartesianPICConfig>) {
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

    SpaceChargeSolver::SpaceChargeSolver(
            SpaceChargeConfig config, std::unique_ptr<SpaceChargeAlgorithm> algorithm,
            std::vector<ParticleFieldBinding3D> bindings,
            std::shared_ptr<const BunchStateHandler> bunchState, std::size_t primaryIndex)
        : config_m(std::move(config)),
          algorithm_m(std::move(algorithm)),
          bindings_m(std::move(bindings)),
          primaryIndex_m(primaryIndex),
          bunchState_m(std::move(bunchState)),
          diagnostics_m(makeDiagnosticSchedule(config_m)) {
        if (algorithm_m == nullptr) {
            throw OpalException(
                    "SpaceChargeSolver::SpaceChargeSolver", "The space-charge algorithm is null.");
        }
        if (bindings_m.empty() || primaryIndex_m >= bindings_m.size()) {
            throw OpalException(
                    "SpaceChargeSolver::SpaceChargeSolver", "The particle binding set is invalid.");
        }
        if (bunchState_m == nullptr) {
            throw OpalException(
                    "SpaceChargeSolver::SpaceChargeSolver", "The bunch state handler is null.");
        }
        for (const ParticleFieldBinding3D& binding : bindings_m) {
            if (!binding.hasCompleteIdentity()) {
                throw OpalException(
                        "SpaceChargeSolver::SpaceChargeSolver",
                        "The particle binding set contains an incomplete binding.");
            }
        }
        capabilities_m = algorithm_m->capabilities();
        validateConfiguration();
    }

    void SpaceChargeSolver::solve(SpaceChargeSolveContext& context) {
        validateBindings(context.particles());
        context.particles().selectForSolve(capabilities_m.particleSelection);
        validateRequest(context);
        algorithm_m->solve(context, diagnostics_m);
    }

    int SpaceChargeSolver::reportedBinCount() const {
        return std::visit(
                [this](const auto& algorithmConfig) {
                    using Config = std::decay_t<decltype(algorithmConfig)>;
                    if constexpr (std::is_same_v<Config, CartesianPICConfig>) {
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

    void SpaceChargeSolver::validateConfiguration() const {
        std::visit(
                [this](const auto& algorithmConfig) {
                    using Config = std::decay_t<decltype(algorithmConfig)>;
                    if constexpr (std::is_same_v<Config, CartesianPICConfig>) {
                        if (algorithmConfig.binning().has_value()
                            && !capabilities_m.supportsBinning) {
                            throw OpalException(
                                    "SpaceChargeSolver::SpaceChargeSolver",
                                    "The selected space-charge algorithm does not support "
                                    "binning.");
                        }
                        if (!capabilities_m.supports(algorithmConfig.correction().kind())) {
                            throw OpalException(
                                    "SpaceChargeSolver::SpaceChargeSolver",
                                    "The selected space-charge algorithm does not support the "
                                    "configured source-plane correction.");
                        }
                        if (algorithmConfig.correction().planeDumpFrequency() != 0
                            && !capabilities_m.supportsPotentialOutput) {
                            throw OpalException(
                                    "SpaceChargeSolver::SpaceChargeSolver",
                                    "The selected space-charge algorithm cannot write the "
                                    "configured potential diagnostic.");
                        }
                        if (algorithmConfig.repartitionFrequency() != 0
                            && !capabilities_m.supportsRedistribution) {
                            throw OpalException(
                                    "SpaceChargeSolver::SpaceChargeSolver",
                                    "The selected space-charge algorithm does not support particle "
                                    "redistribution.");
                        }
                    }
                },
                config_m.algorithmConfig());
    }

    void SpaceChargeSolver::validateBindings(const ParticleFieldSet& particles) const {
        const auto bindings = particles.bindings();
        if (bindings.size() != bindings_m.size() || particles.primaryIndex() != primaryIndex_m) {
            throw OpalException(
                    "SpaceChargeSolver::solve",
                    "The solve context particle binding set does not match construction.");
        }
        for (std::size_t index = 0; index < bindings.size(); ++index) {
            if (!bindings[index].sameIdentity(bindings_m[index])) {
                throw OpalException(
                        "SpaceChargeSolver::solve",
                        "The solve context particle or R/P/E/B binding does not match "
                        "construction at index "
                                + std::to_string(index) + ".");
            }
        }
    }

    void SpaceChargeSolver::validateRequest(const SpaceChargeSolveContext& context) const {
        const SpaceChargeRequest& requested = context.request();
        if (bunchState_m->fixedCartesianDomain().has_value()) {
            if (!capabilities_m.supportsFixedCartesianDomain
                || config_m.algorithmType() != SpaceChargeAlgorithmType::CartesianPIC) {
                throw OpalException(
                        "SpaceChargeSolver::solve",
                        "The selected space-charge algorithm does not support a fixed Cartesian "
                        "domain.");
            }
            const CartesianPICConfig& cartesian = config_m.get<CartesianPICConfig>();
            if (cartesian.backend() != PoissonSolverType::Open) {
                throw OpalException(
                        "SpaceChargeSolver::solve",
                        "A fixed Cartesian domain currently requires the OPEN Poisson backend.");
            }
            if (cartesian.correction().enabled()
                || requested.correction.kind != SpaceChargeCorrectionType::None) {
                throw OpalException(
                        "SpaceChargeSolver::solve",
                        "A fixed Cartesian domain does not support source-plane corrections.");
            }
            if (cartesian.repartitionFrequency() != 0) {
                throw OpalException(
                        "SpaceChargeSolver::solve",
                        "ORB redistribution must be disabled while a fixed Cartesian domain is "
                        "active.");
            }
        }
        if (requested.useBinning && !capabilities_m.supportsBinning) {
            throw OpalException(
                    "SpaceChargeSolver::solve",
                    "The solve requests binning, but the selected algorithm does not support it.");
        }
        if (requested.writePotential && !capabilities_m.supportsPotentialOutput) {
            throw OpalException(
                    "SpaceChargeSolver::solve",
                    "The solve requests potential output, but the selected algorithm does not "
                    "support it.");
        }
        if (!capabilities_m.supports(requested.correction.kind)) {
            throw OpalException(
                    "SpaceChargeSolver::solve",
                    "The solve requests a correction that the selected algorithm does not "
                    "support.");
        }
        const SpaceChargeRequest expected =
                SpaceChargeRequestSchedule(config_m).requestForStep(context.stepState().step);
        const bool wrongPlane = requested.correction.kind != SpaceChargeCorrectionType::None
                                && requested.correction.planeZ != expected.correction.planeZ;
        if (requested.useBinning != expected.useBinning
            || requested.writePotential != expected.writePotential
            || requested.correction.kind != expected.correction.kind || wrongPlane) {
            throw OpalException(
                    "SpaceChargeSolver::solve",
                    "The requested physics differs from the immutable step-resolved "
                    "configuration.");
        }
    }

}  // namespace opalx::spacecharge
