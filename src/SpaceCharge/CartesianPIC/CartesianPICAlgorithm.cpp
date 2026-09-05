/**
 * @file CartesianPICAlgorithm.cpp
 * @brief Implements the complete Cartesian 3D PIC space-charge algorithm.
 */

#include "SpaceCharge/CartesianPIC/CartesianPICAlgorithm.h"

#include "PartBunch/BunchStateHandler.h"
#include "Structure/DataSink.h"
#include "Utilities/OpalException.h"

#include <iomanip>
#include <optional>
#include <string>
#include <utility>

namespace opalx::spacecharge {
    namespace {

        using ParticleContainer = CartesianPICAlgorithm::ParticleContainer;

        PoissonFieldBinding makePoissonFieldBinding(
                CartesianPICAlgorithm::FieldStorage& fieldStorage) {
            return {&fieldStorage.chargeDensity(), &fieldStorage.electricField()};
        }

        ParticleContainer& requirePrimaryParticles(std::span<ParticleContainer* const> particles) {
            if (particles.empty() || particles.front() == nullptr) {
                throw OpalException(
                        "CartesianPICAlgorithm::CartesianPICAlgorithm",
                        "The primary particle container is not available.");
            }
            return *particles.front();
        }

    }  // namespace

    CartesianPICAlgorithm::CartesianPICAlgorithm(
            CartesianPICConfig config, std::span<ParticleContainer* const> particles,
            std::unique_ptr<FieldStorage> fieldStorage, DataSink* dataSink,
            std::shared_ptr<const BunchStateHandler> bunchState)
        : config_m(config),
          primary_m(&requirePrimaryParticles(particles)),
          bunchState_m(std::move(bunchState)),
          fieldStorage_m(std::move(fieldStorage)),
          dataSink_m(dataSink),
          domainUpdater_m(config, particles) {
        if (fieldStorage_m == nullptr) {
            throw OpalException(
                    "CartesianPICAlgorithm::CartesianPICAlgorithm",
                    "The Cartesian PIC field storage is null.");
        }
        if (dataSink_m == nullptr) {
            throw OpalException(
                    "CartesianPICAlgorithm::CartesianPICAlgorithm", "The data sink is null.");
        }
        if (bunchState_m == nullptr) {
            throw OpalException(
                    "CartesianPICAlgorithm::CartesianPICAlgorithm",
                    "The bunch state handler is null.");
        }
        validateSpaceChargeConfig(SpaceChargeConfig(config_m));
        if (config_m.binning.has_value()) {
            particleBinTraversal_m =
                    std::make_unique<ParticleBinTraversalType>(*primary_m, *config_m.binning);
            binStats_m.reserve(particleBinTraversal_m->maximumBinCount());
        }

        fieldStorage_m->initializeFields(config_m.backend);
        poissonSolver_m = std::make_unique<PoissonSolver>(
                makePoissonSolverConfig(config_m), makePoissonFieldBinding(*fieldStorage_m));
        if (config_m.backend == PoissonSolverType::P3M) {
            shortRangeInteraction_m.emplace(config.p3mCutoff);
        }
        // Preserve the construction-time planning solve, but reset the Poisson solver's runtime
        // count so diagnostics and debug-file numbering begin with the first physical solve.
        poissonSolver_m->warmup();
    }

    CartesianPICAlgorithm::SolvePlan CartesianPICAlgorithm::makeSolvePlan(std::size_t step) const {
        SolvePlan plan;
        const CorrectionConfig& configured = config_m.correction;
        plan.correctionExpired             = configured.enabled() && configured.maximumSteps != 0
                                 && step >= configured.maximumSteps;
        plan.activeCorrection = plan.correctionExpired ? CorrectionConfig() : configured;

        if (config_m.binning.has_value()
            || plan.activeCorrection.kind == SpaceChargeCorrectionType::ShiftedGreen) {
            plan.passes[plan.passCount++] = PassKind::Primary;
            if (plan.activeCorrection.kind == SpaceChargeCorrectionType::ImageCharge) {
                plan.passes[plan.passCount++] = PassKind::Image;
            } else if (plan.activeCorrection.kind == SpaceChargeCorrectionType::ShiftedGreen) {
                plan.passes[plan.passCount++] = PassKind::ShiftedImage;
            }
        } else {
            plan.passes[plan.passCount++] =
                    plan.activeCorrection.kind == SpaceChargeCorrectionType::ImageCharge
                            ? PassKind::PrimaryAndImage
                            : PassKind::Primary;
        }
        return plan;
    }

    CartesianPICAlgorithm::PassProperties CartesianPICAlgorithm::passProperties(
            PassKind pass, double planeZ, bool binned) {
        using DepositKind = ParticleMeshTransfer::DepositKind;
        switch (pass) {
            case PassKind::Primary:
                return {.depositKind = DepositKind::Primary, .suppressFieldDump = binned};
            case PassKind::PrimaryAndImage:
                return {.depositKind             = DepositKind::PrimaryAndImage,
                        .imagePolicy             = {true, planeZ},
                        .dumpDirichletPlaneAfter = true,
                        .label                   = "primary-and-image"};
            case PassKind::Image:
                return {.depositKind       = DepositKind::Image,
                        .imagePolicy       = {true, planeZ},
                        .suppressFieldDump = true,
                        .magneticSign      = -1.0,
                        .label             = "image"};
            case PassKind::ShiftedImage:
                return {.shiftedGreen = true,
                        .magneticSign = -1.0,
                        .sourceRule   = FieldSourceRule::ShiftedGreenImageZ,
                        .label        = "shifted-green"};
        }
        throw OpalException("CartesianPICAlgorithm::passProperties", "Unknown solve pass.");
    }

    SpaceChargeSolveResult CartesianPICAlgorithm::solve(const SpaceChargeSolveContext& context) {
        SpaceChargeSolveResult result;
        const SolvePlan plan = makeSolvePlan(context.stepState().step);

        const auto& fixedState = bunchState_m->fixedCartesianDomain();
        const bool fixedDomain = fixedState.has_value();
        if (fixedDomain) {
            if (config_m.backend != PoissonSolverType::Open) {
                throw OpalException(
                        "CartesianPICAlgorithm::solve",
                        "A fixed Cartesian domain currently requires the OPEN Poisson backend.");
            }
            if (config_m.correction.enabled()) {
                throw OpalException(
                        "CartesianPICAlgorithm::solve",
                        "A fixed Cartesian domain does not support source-plane corrections.");
            }
            if (config_m.repartitionFrequency != 0) {
                throw OpalException(
                        "CartesianPICAlgorithm::solve",
                        "ORB redistribution must be disabled while a fixed Cartesian domain is "
                        "active.");
            }
        }
        clearSelfFields(*primary_m);
        enterSolveFrame(context.stepState().frames, *primary_m);
        result.redistributions += domainUpdater_m.updateForSolve(
                DomainCoordinateFrame::Beam, context, plan.activeCorrection,
                fixedDomain ? &*fixedState : nullptr, *fieldStorage_m, *poissonSolver_m);

        solveInBeamFrame(context, plan, result);

        leaveSolveFrame(context.stepState().frames, *primary_m);
        if (fixedDomain) {
            // Keep the fixed beam-frame mesh and decomposition for the next interaction. Only the
            // restored primary coordinates changed, so refresh moments without migrating them.
            primary_m->markMomentsDirty();
            primary_m->updateMoments();
        } else {
            result.redistributions += domainUpdater_m.updateForSolve(
                    DomainCoordinateFrame::Reference, context, plan.activeCorrection, nullptr,
                    *fieldStorage_m, *poissonSolver_m);
        }
        return result;
    }

    void CartesianPICAlgorithm::solveInBeamFrame(
            const SpaceChargeSolveContext& context, const SolvePlan& plan,
            SpaceChargeSolveResult& result) {
        Inform m("CartesianPICAlgorithm::solveInBeamFrame");
        const auto& poissonCapabilities = poissonSolver_m->capabilities();
        if (poissonCapabilities.isNoOp) {
            m << level5 << "Skipping scatter/gather and space-charge computation for NONE solver."
              << endl;
            return;
        }

        // Use the global count because early emission can leave this rank empty while another rank
        // owns the only emitted particle. A single global charge has no space-charge field here.
        if (primary_m->getTotalNum() <= 1) {
            return;
        }
        // Deposition temporarily stores dt*Q in the particle time-step attribute. Reject Q=0
        // before that operation because restoring dt would otherwise evaluate 0/0 and leave NaNs
        // for later passes. In practice this nearly always indicates a missing BEAM BCHARGE.
        if (primary_m->getChargePerParticle() == 0.0) {
            throw OpalException(
                    "CartesianPICAlgorithm::solveInBeamFrame",
                    "Per-particle charge is zero but a space-charge solver is active (type="
                            + std::string(poissonSolver_m->name())
                            + "). This almost always means the BEAM command is missing BCHARGE. "
                              "Set BCHARGE on the BEAM definition, or use TYPE=NONE when no "
                              "space charge is intended.");
        }
        if (plan.passCount == 0) {
            throw OpalException(
                    "CartesianPICAlgorithm::solveInBeamFrame",
                    "The prepared correction contains no primary solve pass.");
        }

        const bool binned = particleBinTraversal_m != nullptr;
        m << level4 << "Entry: rank=" << ippl::Comm->rank()
          << ", localParticles=" << primary_m->getLocalNum()
          << ", totalParticles=" << primary_m->getTotalNum() << ", hasBins=" << (binned ? 1 : 0)
          << ", stype=" << poissonSolver_m->name() << endl;

        if (plan.correctionExpired) {
            m << level3 << "ZEROFACE_MAXSTEPS reached (step=" << context.stepState().step
              << ", maxSteps=" << config_m.correction.maximumSteps << "); disabling ";
            if (config_m.correction.kind == SpaceChargeCorrectionType::ShiftedGreen) {
                m << "SHIFTED_GREENS_FUNCTION correction for this step." << endl;
            } else {
                m << "image charges for this step." << endl;
            }
        }

        if (binned) {
            solveBinned(context, plan, result);
        } else {
            solveWholeBunch(context, plan, result);
        }
    }

    void CartesianPICAlgorithm::solveWholeBunch(
            const SpaceChargeSolveContext& context, const SolvePlan& plan,
            SpaceChargeSolveResult& result) {
        result.reportedBins = 1;

        const bool accumulatePasses = plan.passCount > 1;
        if (accumulatePasses) {
            relativisticFieldComposer_m.clearAccumulation(*fieldStorage_m);
        }
        for (std::size_t index = 0; index < plan.passCount; ++index) {
            solvePass(context, plan, nullptr, plan.passes[index], result);
        }
        if (accumulatePasses) {
            relativisticFieldComposer_m.gatherAccumulated(
                    particleMeshTransfer_m, primary_m->E, primary_m->B, primary_m->R,
                    *fieldStorage_m);
        }
    }

    void CartesianPICAlgorithm::solveBinned(
            const SpaceChargeSolveContext& context, const SolvePlan& plan,
            SpaceChargeSolveResult& result) {
        const long long step = static_cast<long long>(context.stepState().step);
        const bool captureSnapshots =
                dataSink_m != nullptr && config_m.binning.has_value()
                && !config_m.binning->dumpFile.empty() && config_m.binning->dumpFrequency > 0
                && step % static_cast<long long>(config_m.binning->dumpFrequency) == 0;

        // A pre-merge snapshot records the freshly rebuilt histogram. Adaptive binning supplies a
        // second snapshot after merging; fixed binning deliberately has no post-merge snapshot.
        const BinPreparationResult prepared = particleBinTraversal_m->prepareBins(captureSnapshots);
        if (prepared.beforeMerge.has_value()) {
            dumpBinSnapshot(context, *prepared.beforeMerge, true);
        }
        if (prepared.afterMerge.has_value()) {
            dumpBinSnapshot(context, *prepared.afterMerge, false);
        }
        result.reportedBins = prepared.mergedBinCount == config_m.binning->maximumBins
                                      ? 1
                                      : static_cast<int>(prepared.mergedBinCount);
        relativisticFieldComposer_m.clearAccumulation(*fieldStorage_m);

        Inform m("CartesianPICAlgorithm::solveBinned");
        m << level4 << "Iteration mode=binned, nBins=" << static_cast<int>(prepared.mergedBinCount)
          << ", stype=" << poissonSolver_m->name() << endl;

        binStats_m.clear();
        while (const std::optional<ParticleBinType> unit =
                       particleBinTraversal_m->nextNonemptyBin()) {
            m << level4 << "binIndex=" << static_cast<int>(unit->ordinal)
              << " nPartGlobal=" << static_cast<unsigned long long>(unit->globalParticleCount)
              << " gammaBin=" << std::setprecision(10) << unit->gamma << endl;
            binStats_m.push_back(
                    BinStatsRow{
                            static_cast<long long>(unit->ordinal),
                            static_cast<unsigned long long>(unit->globalParticleCount),
                            unit->gamma});

            for (std::size_t index = 0; index < plan.passCount; ++index) {
                solvePass(context, plan, &*unit, plan.passes[index], result);
            }
        }

        relativisticFieldComposer_m.gatherAccumulated(
                particleMeshTransfer_m, primary_m->E, primary_m->B, primary_m->R, *fieldStorage_m);
        if (config_m.binning->tablePrintFrequency > 0
            && step % static_cast<long long>(config_m.binning->tablePrintFrequency) == 0) {
            printBinStatsTable();
        }
    }

    void CartesianPICAlgorithm::solvePass(
            const SpaceChargeSolveContext& context, const SolvePlan& plan,
            const ParticleBinType* unit, PassKind pass, SpaceChargeSolveResult& result) {
        const double planeZ         = plan.activeCorrection.planeZ;
        const PassProperties policy = passProperties(pass, planeZ, unit != nullptr);

        const auto& capabilities = poissonSolver_m->capabilities();
        if (unit == nullptr) {
            ParticleMeshTransfer::ChargeNormalization normalization;
            normalization.timeStep              = context.stepState().timeStep;
            normalization.gamma                 = 1.0;
            normalization.selectedCharge        = primary_m->getTotalCharge();
            normalization.couplingConstant      = poissonSolver_m->couplingConstant();
            normalization.normalizeByCellVolume = capabilities.normalizeChargeByCellVolume;
            normalization.subtractNeutralizingBackground =
                    capabilities.subtractNeutralizingBackground;
            particleMeshTransfer_m.depositCharge(
                    *primary_m, *fieldStorage_m, policy.depositKind,
                    ParticleMeshTransfer::Selection::direct(0, primary_m->getLocalNum()),
                    normalization, policy.imagePolicy);
        } else {
            depositChargeForBin(context, *unit, policy);
        }

        fieldStorage_m->electricField() = 0.0;
        auto& mesh                      = fieldStorage_m->mesh();
        const auto originalSpacing      = mesh.getMeshSpacing();
        auto stretchedSpacing           = originalSpacing;
        const double gamma              = unit == nullptr ? 1.0 : unit->gamma;

        // Each binned Poisson solve is performed in the bin rest frame. Stretching the
        // longitudinal mesh spacing by gamma applies the Lorentz contraction convention used by
        // the legacy solver; the original spacing is restored after successful composition.
        stretchedSpacing[2] *= gamma;

        PoissonSolveRequest poissonRequest;
        if (policy.shiftedGreen) {
            const auto origin = mesh.getOrigin();
            const int longitudinalExtent =
                    static_cast<int>(fieldStorage_m->layout().getDomain()[2].length());
            const double zCenter =
                    origin[2] + 0.5 * static_cast<double>(longitudinalExtent) * originalSpacing[2];
            ippl::Vector<double, 3> shift(0.0);

            // Old OPAL expresses this displacement with the opposite sign. IPPL evaluates
            // G(r-shift), so the translated kernel requires +2*gamma*(zCenter-planeZ).
            shift[2]                          = 2.0 * gamma * (zCenter - planeZ);
            poissonRequest.greenFunctionShift = shift;
        }

        if (unit != nullptr) {
            mesh.setMeshSpacing(stretchedSpacing);
        }

        Inform m("CartesianPICAlgorithm::solvePass");
        m << level4 << "pass=" << policy.label
          << ", binIndex=" << static_cast<int>(unit == nullptr ? 0 : unit->ordinal)
          << ", suppressFieldDump=" << (policy.suppressFieldDump ? 1 : 0);
        if (poissonRequest.hasShiftedGreenFunction()) {
            m << ", plane=" << planeZ << ", shift_z=" << (*poissonRequest.greenFunctionShift)[2];
        }
        m << endl;

        poissonSolver_m->solve(poissonRequest, {.suppressFieldDump = policy.suppressFieldDump});
        ++result.backendSolves;

        if (unit == nullptr && policy.dumpDirichletPlaneAfter) {
            dumpDirichletPlaneDiagnosticsIfRequested(context, "legacy", planeZ);
        }

        if (unit == nullptr && plan.passCount == 1) {
            relativisticFieldComposer_m.gatherElectrostatic(
                    particleMeshTransfer_m, primary_m->E, primary_m->R, *fieldStorage_m);
            if (shortRangeInteraction_m.has_value()) {
                shortRangeInteraction_m->apply(*primary_m);
            }
        } else {
            RelativisticFieldComposerType::Policy compositionPolicy;
            if (unit != nullptr) {
                compositionPolicy.meanMomentum = unit->meanMomentum;
                compositionPolicy.gamma        = unit->gamma;
            }
            compositionPolicy.magneticSign = policy.magneticSign;
            compositionPolicy.sourceRule   = policy.sourceRule;

            // Shifted Green reuses the real-charge RHS. Its image sign enters through field
            // reflection and component signs here; negating rho would invert the image twice.
            relativisticFieldComposer_m.accumulate(*fieldStorage_m, compositionPolicy);
        }

        if (unit != nullptr) {
            mesh.setMeshSpacing(originalSpacing);
        }
    }

    void CartesianPICAlgorithm::depositChargeForBin(
            const SpaceChargeSolveContext& context, const ParticleBinType& unit,
            const PassProperties& pass) {
        Inform m("CartesianPICAlgorithm::depositChargeForBin");
        const auto& indexedSelection = unit.indexedSelection;
        const auto& selectionPolicy  = indexedSelection.policy();
        const auto& hash             = indexedSelection.hash();
        const std::size_t begin      = static_cast<std::size_t>(selectionPolicy.begin());
        const std::size_t end        = static_cast<std::size_t>(selectionPolicy.end());
        const std::size_t hashExtent = static_cast<std::size_t>(hash.extent(0));
        const std::size_t localCount = primary_m->getLocalNum();

        if (end > hashExtent) {
            throw OpalException(
                    "CartesianPICAlgorithm::depositChargeForBin",
                    "Bin scatter policy exceeds its hash extent.");
        }
        if (unit.coversAllLocalParticles != (begin == 0 && end == localCount)) {
            throw OpalException(
                    "CartesianPICAlgorithm::depositChargeForBin",
                    "The solve-unit all-local flag does not match its indexed selection.");
        }

        // A merged bin often covers every local particle during early emission. Use the direct
        // scatter path in that case to avoid dereferencing the bin hash in the hot kernel.
        const char* depositName = "PrimaryAndImage";
        if (pass.depositKind == ParticleMeshTransfer::DepositKind::Primary) {
            depositName = "PrimaryOnly";
        } else if (pass.depositKind == ParticleMeshTransfer::DepositKind::Image) {
            depositName = "ImageOnly";
        }
        m << level5 << "prepareRho: scatter mode=" << depositName
          << ", path=" << (unit.coversAllLocalParticles ? "all-local" : "subset")
          << ", localP=" << static_cast<unsigned long long>(localCount) << ", policy=[" << begin
          << "," << end << ")"
          << ", hashExtent=" << static_cast<unsigned long long>(hashExtent) << endl;

        const auto& capabilities = poissonSolver_m->capabilities();
        ParticleMeshTransfer::ChargeNormalization normalization;
        normalization.timeStep = context.stepState().timeStep;
        normalization.gamma    = unit.gamma;
        normalization.selectedCharge =
                primary_m->getChargePerParticle() * static_cast<double>(unit.globalParticleCount);
        normalization.couplingConstant               = poissonSolver_m->couplingConstant();
        normalization.normalizeByCellVolume          = capabilities.normalizeChargeByCellVolume;
        normalization.subtractNeutralizingBackground = capabilities.subtractNeutralizingBackground;

        particleMeshTransfer_m.depositCharge(
                *primary_m, *fieldStorage_m, pass.depositKind, unit.depositSelection(),
                normalization, pass.imagePolicy);
    }

    void CartesianPICAlgorithm::dumpBinSnapshot(
            const SpaceChargeSolveContext& context, const BinConfigurationSnapshot& snapshot,
            bool beforeMerge) const {
        const std::vector<std::size_t> particleCounts(
                snapshot.particleCounts.begin(), snapshot.particleCounts.end());
        dataSink_m->dumpBinConfig(
                static_cast<long long>(context.stepState().step), context.stepState().time,
                beforeMerge, particleCounts, snapshot.widths, snapshot.lowerBound,
                config_m.binning->dumpFile);
    }

    void CartesianPICAlgorithm::dumpDirichletPlaneDiagnosticsIfRequested(
            const SpaceChargeSolveContext& context, const std::string& solveTag, double planeZ) {
        const long long step        = static_cast<long long>(context.stepState().step);
        const std::size_t frequency = config_m.correction.planeDumpFrequency;
        if (frequency == 0 || step < 0 || step % static_cast<long long>(frequency) != 0) {
            return;
        }
        Inform m("CartesianPICAlgorithm::dumpDirichletPlaneDiagnosticsIfRequested");
        if (ippl::Comm->size() != 1) {
            if (!warnedPlaneDumpParallelUnsupported_m) {
                warnedPlaneDumpParallelUnsupported_m = true;
                m << level3
                  << "Dirichlet-plane diagnostics currently support only single-rank runs. "
                     "Skipping dump and statistics output."
                  << endl;
            }
            return;
        }
        if (dataSink_m == nullptr) {
            return;
        }

        const auto planeDiagnostics = dataSink_m->dumpDirichletPlane(
                step, context.stepState().time, planeZ, fieldStorage_m->chargeDensity(), solveTag);
        if (planeDiagnostics.sampleCount == 0) {
            return;
        }
        m << level2 << "Dirichlet-plane potential diagnostics (" << solveTag << ") at step " << step
          << ": z=" << planeZ << " m, mean(phi)=" << planeDiagnostics.mean
          << " V, var(phi)=" << planeDiagnostics.variance << " V^2" << endl;
    }

    void CartesianPICAlgorithm::printBinStatsTable() const {
        const std::string& diagnosticName = particleBinTraversal_m->diagnosticName();
        const std::string informName =
                diagnosticName.empty()
                        ? "CartesianPICAlgorithm::printBinStatsTable"
                        : "CartesianPICAlgorithm::printBinStatsTable[" + diagnosticName + "]";
        Inform m(informName.c_str());
        m << level2 << std::setw(9) << "bin"
          << " | " << std::setw(13) << "nParticles"
          << " | " << std::setw(15) << "gammaBin" << endl;
        for (const BinStatsRow& row : binStats_m) {
            m << level2 << std::setw(9) << row.binNumber << " | " << std::setw(13)
              << row.particleCount << " | " << std::setw(15) << std::setprecision(10) << row.gamma
              << endl;
        }
    }

}  // namespace opalx::spacecharge
