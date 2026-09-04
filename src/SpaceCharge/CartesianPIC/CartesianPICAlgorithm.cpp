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

        std::string_view poissonSolverName(PoissonSolverType kind) {
            switch (kind) {
                case PoissonSolverType::None:
                    return "NONE";
                case PoissonSolverType::PeriodicFFT:
                    return "FFT";
                case PoissonSolverType::Open:
                    return "OPEN";
                case PoissonSolverType::ConjugateGradient:
                    return "CG";
                case PoissonSolverType::P3M:
                    return "P3M";
            }
            throw OpalException(
                    "CartesianPICAlgorithm::CartesianPICAlgorithm", "Unknown Poisson backend.");
        }

        PoissonFieldBinding makePoissonFieldBinding(
                CartesianPICAlgorithm::FieldStorage& fieldStorage) {
            return {&fieldStorage.chargeDensity(), &fieldStorage.electricField()};
        }

        ParticleContainer& requirePrimaryParticles(
                std::span<const ParticleFieldBinding3D> bindings) {
            if (bindings.empty() || bindings.front().container == nullptr) {
                throw OpalException(
                        "CartesianPICAlgorithm::CartesianPICAlgorithm",
                        "The primary particle container is not available.");
            }
            return *bindings.front().container;
        }

    }  // namespace

    CartesianPICAlgorithm::CartesianPICAlgorithm(
            CartesianPICConfig config, std::span<const ParticleFieldBinding3D> particleBindings,
            std::unique_ptr<FieldStorage> fieldStorage, DataSink* dataSink,
            std::shared_ptr<const BunchStateHandler> bunchState)
        : primary_m(&requirePrimaryParticles(particleBindings)),
          bunchState_m(std::move(bunchState)),
          fieldStorage_m(std::move(fieldStorage)),
          poissonSolverType_m(config.backend()),
          dataSink_m(dataSink),
          binningConfig_m(config.binning()),
          correctionPassSchedule_m(config.correction(), config.binning().has_value()),
          particleDomainOperations_m(particleBindings),
          domainUpdater_m(config, bunchState_m) {
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
        if (poissonSolverType_m == PoissonSolverType::ConjugateGradient) {
            throw OpalException(
                    "CartesianPICAlgorithm::CartesianPICAlgorithm",
                    "The CG Poisson backend is recognized but not implemented.");
        }
        if (binningConfig_m.has_value()) {
            particleBinTraversal_m =
                    std::make_unique<ParticleBinTraversalType>(*primary_m, *binningConfig_m);
            binStats_m.reserve(particleBinTraversal_m->maximumBinCount());
        }

        fieldStorage_m->initializeFields(poissonSolverName(poissonSolverType_m));
        PoissonSolverConfig poissonConfig;
        poissonConfig.type               = poissonSolverType_m;
        poissonConfig.greenFunction      = config.greenFunction();
        poissonConfig.p3mCutoff          = config.p3mCutoff();
        poissonConfig.boundaryConditions = config.boundaryConditions();
        poissonSolver_m                  = std::make_unique<PoissonSolver>(
                poissonConfig, makePoissonFieldBinding(*fieldStorage_m));
        if (poissonSolverType_m == PoissonSolverType::P3M) {
            shortRangeInteraction_m.emplace(config.p3mCutoff());
        }
        // Preserve the construction-time planning solve, but reset the Poisson solver's runtime
        // count so diagnostics and debug-file numbering begin with the first physical solve.
        poissonSolver_m->warmup();
    }

    SpaceChargeCapabilities CartesianPICAlgorithm::capabilities() const {
        SpaceChargeCapabilities result;
        result.particleSelection            = ParticleSelectionMode::PrimaryOnly;
        result.supportsBinning              = true;
        result.supportsImageCharge          = true;
        result.supportsShiftedGreen         = true;
        result.supportsRedistribution       = true;
        result.supportsPotentialOutput      = true;
        result.supportsFixedCartesianDomain = true;
        return result;
    }

    void CartesianPICAlgorithm::solve(
            SpaceChargeSolveContext& context, SpaceChargeDiagnostics& diagnostics) {
        const CorrectionPassSequence correction =
                correctionPassSchedule_m.passesForStep(context.request(), context.stepState().step);

        const bool fixedDomain = bunchState_m->fixedCartesianDomain().has_value();
        SpaceChargeFrameTransform<double, 3> frameTransform(context.stepState().frames, *primary_m);
        frameTransform.enter();
        domainUpdater_m.updateForSolve(
                DomainCoordinateFrame::Beam, context, *fieldStorage_m, particleDomainOperations_m,
                *poissonSolver_m, diagnostics);

        frameTransform.markComputedFields();
        solveInBeamFrame(context, correction, diagnostics);

        frameTransform.leave();
        if (fixedDomain) {
            // Keep the fixed beam-frame mesh and decomposition for the next interaction. Only the
            // restored primary coordinates changed, so refresh moments without migrating them.
            primary_m->markMomentsDirty();
            particleDomainOperations_m.updatePrimaryMoments();
        } else {
            domainUpdater_m.updateForSolve(
                    DomainCoordinateFrame::Reference, context, *fieldStorage_m,
                    particleDomainOperations_m, *poissonSolver_m, diagnostics);
        }
    }

    void CartesianPICAlgorithm::solveInBeamFrame(
            SpaceChargeSolveContext& context, const CorrectionPassSequence& correction,
            SpaceChargeDiagnostics& diagnostics) {
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
            Kokkos::deep_copy(primary_m->E.getView(), Vector_t<double, 3>(0.0));
            Kokkos::deep_copy(primary_m->B.getView(), Vector_t<double, 3>(0.0));
            return;
        }
        // Deposition temporarily stores dt*Q in the particle time-step attribute. Reject Q=0
        // before that operation because restoring dt would otherwise evaluate 0/0 and leave NaNs
        // for later passes. In practice this nearly always indicates a missing BEAM BCHARGE.
        if (primary_m->getChargePerParticle() == 0.0) {
            throw OpalException(
                    "CartesianPICAlgorithm::solveInBeamFrame",
                    "Per-particle charge is zero but a space-charge solver is active (type="
                            + std::string(poissonSolverName(poissonSolverType_m))
                            + "). This almost always means the BEAM command is missing BCHARGE. "
                              "Set BCHARGE on the BEAM definition, or use TYPE=NONE when no "
                              "space charge is intended.");
        }
        if (correction.passCount == 0) {
            throw OpalException(
                    "CartesianPICAlgorithm::solveInBeamFrame",
                    "The prepared correction contains no primary solve pass.");
        }

        const bool binned = particleBinTraversal_m != nullptr;
        m << level4 << "Entry: rank=" << ippl::Comm->rank()
          << ", localParticles=" << primary_m->getLocalNum()
          << ", totalParticles=" << primary_m->getTotalNum() << ", hasBins=" << (binned ? 1 : 0)
          << ", stype=" << poissonSolverName(poissonSolverType_m) << endl;

        if (correction.correctionExpired) {
            m << level3 << "ZEROFACE_MAXSTEPS reached (step=" << context.stepState().step
              << ", maxSteps=" << correction.maximumSteps << "); disabling ";
            if (correction.configuredCorrection.kind == SpaceChargeCorrectionType::ShiftedGreen) {
                m << "SHIFTED_GREENS_FUNCTION correction for this step." << endl;
            } else {
                m << "image charges for this step." << endl;
            }
        }

        if (!binned && correction.shiftedIgnored) {
            m << level3 << "SHIFTED_GREENS_FUNCTION is set but no binning is active; "
              << "the whole-bunch path does not apply the Dirichlet correction." << endl;
        }

        if (binned) {
            solveBinned(context, correction, diagnostics);
        } else {
            solveWholeBunch(context, correction, diagnostics);
        }
    }

    void CartesianPICAlgorithm::solveWholeBunch(
            SpaceChargeSolveContext& context, const CorrectionPassSequence& correction,
            SpaceChargeDiagnostics& diagnostics) {
        diagnostics.recordBinCount(0);
        bool dumpedDirichletPlaneThisStep = false;

        // The direct image-charge mode deposits primary and mirrored charge into one RHS and uses
        // one standard solve. Shifted Green remains a binned-only correction for compatibility.
        for (const CartesianPICPass pass : correction) {
            solvePass(context, nullptr, pass, dumpedDirichletPlaneThisStep, diagnostics);
        }
    }

    void CartesianPICAlgorithm::solveBinned(
            SpaceChargeSolveContext& context, const CorrectionPassSequence& correction,
            SpaceChargeDiagnostics& diagnostics) {
        const long long step = static_cast<long long>(context.stepState().step);
        const bool captureSnapshots =
                dataSink_m != nullptr && binningConfig_m.has_value()
                && !binningConfig_m->dumpFile().empty() && binningConfig_m->dumpFrequency() > 0
                && step % static_cast<long long>(binningConfig_m->dumpFrequency()) == 0;

        // A pre-merge snapshot records the freshly rebuilt histogram. Adaptive binning supplies a
        // second snapshot after merging; fixed binning deliberately has no post-merge snapshot.
        const BinPreparationResult prepared = particleBinTraversal_m->prepareBins(captureSnapshots);
        if (prepared.beforeMerge.has_value()) {
            dumpBinSnapshot(context, *prepared.beforeMerge, true);
        }
        if (prepared.afterMerge.has_value()) {
            dumpBinSnapshot(context, *prepared.afterMerge, false);
        }
        diagnostics.recordBinCount(prepared.mergedBinCount);
        relativisticFieldComposer_m.clearAccumulation(*fieldStorage_m);

        Inform m("CartesianPICAlgorithm::solveBinned");
        m << level4 << "Iteration mode=binned, nBins=" << static_cast<int>(prepared.mergedBinCount)
          << ", stype=" << poissonSolverName(poissonSolverType_m) << endl;

        binStats_m.clear();
        bool dumpedDirichletPlaneThisStep = false;
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

            for (const CartesianPICPass pass : correction) {
                solvePass(context, &*unit, pass, dumpedDirichletPlaneThisStep, diagnostics);
            }
        }

        relativisticFieldComposer_m.gatherAccumulated(
                particleMeshTransfer_m, primary_m->E, primary_m->B, primary_m->R, *fieldStorage_m);
        if (diagnostics.shouldPrintBinTable(step)) {
            printBinStatsTable();
        }
    }

    void CartesianPICAlgorithm::solvePass(
            SpaceChargeSolveContext& context, const ParticleBinType* unit, CartesianPICPass pass,
            bool& dumpedDirichletPlaneThisStep, SpaceChargeDiagnostics& diagnostics) {
        const double planeZ         = context.request().correction.planeZ;
        const PassProperties policy = cartesianPICPassProperties<double, 3>(pass, planeZ);
        if (policy.binned != (unit != nullptr)) {
            throw OpalException(
                    "CartesianPICAlgorithm::solvePass",
                    "The solve-pass tag does not match the particle traversal.");
        }

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
        if (policy.backendRule == BackendSolveRule::ShiftedGreen) {
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
        m << level4 << "pass=" << cartesianPICPassLabel(pass)
          << ", binIndex=" << static_cast<int>(unit == nullptr ? 0 : unit->ordinal)
          << ", suppressFieldDump=" << (policy.suppressFieldDump ? 1 : 0);
        if (poissonRequest.hasShiftedGreenFunction()) {
            m << ", plane=" << planeZ << ", shift_z=" << (*poissonRequest.greenFunctionShift)[2];
        }
        m << endl;

        poissonSolver_m->solve(poissonRequest, {.suppressFieldDump = policy.suppressFieldDump});
        diagnostics.recordBackendSolve();

        if (unit == nullptr && policy.dumpDirichletPlaneAfter) {
            dumpDirichletPlaneDiagnosticsIfRequested(context, "legacy", planeZ, diagnostics);
            dumpedDirichletPlaneThisStep = true;
        }

        if (unit == nullptr) {
            relativisticFieldComposer_m.gatherElectrostatic(
                    particleMeshTransfer_m, primary_m->E, primary_m->R, *fieldStorage_m);
            if (shortRangeInteraction_m.has_value()) {
                shortRangeInteraction_m->apply(*primary_m);
            }
        } else {
            RelativisticFieldComposerType::Policy compositionPolicy;
            compositionPolicy.meanMomentum = unit->meanMomentum;
            compositionPolicy.gamma        = unit->gamma;
            compositionPolicy.magneticSign = policy.magneticSign;
            compositionPolicy.sourceRule   = policy.sourceRule;

            // Shifted Green reuses the real-charge RHS. Its image sign enters through field
            // reflection and component signs here; negating rho would invert the image twice.
            relativisticFieldComposer_m.accumulate(*fieldStorage_m, compositionPolicy);
        }

        if (unit != nullptr) {
            mesh.setMeshSpacing(originalSpacing);
        }

        // The explicit-image mesh contains the source plane, so evaluating the Dirichlet residual
        // is meaningful. A shifted-Green correction may solve on a domain far from that plane and
        // therefore never requests this diagnostic.
        if (unit != nullptr && policy.dumpDirichletPlaneAfter && !dumpedDirichletPlaneThisStep) {
            dumpDirichletPlaneDiagnosticsIfRequested(context, "binned", planeZ, diagnostics);
            dumpedDirichletPlaneThisStep = true;
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
                binningConfig_m->dumpFile());
    }

    void CartesianPICAlgorithm::dumpDirichletPlaneDiagnosticsIfRequested(
            const SpaceChargeSolveContext& context, const std::string& solveTag, double planeZ,
            SpaceChargeDiagnostics& diagnostics) {
        const long long step = static_cast<long long>(context.stepState().step);
        if (!diagnostics.shouldDumpPlane(step)) {
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
