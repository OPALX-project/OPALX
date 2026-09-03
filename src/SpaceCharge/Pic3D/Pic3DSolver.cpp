/**
 * @file Pic3DSolver.cpp
 * @brief Implements the complete Cartesian 3D PIC self-field algorithm.
 */

#include "SpaceCharge/Pic3D/Pic3DSolver.h"

#include "Structure/DataSink.h"
#include "Utilities/OpalException.h"

#include <exception>
#include <iomanip>
#include <optional>
#include <string>
#include <utility>

namespace opalx::spacecharge {
    namespace {

        using ParticleContainer = Pic3DSolver::ParticleContainer;

        std::string_view backendName(PoissonBackendKind kind) {
            switch (kind) {
                case PoissonBackendKind::None:
                    return "NONE";
                case PoissonBackendKind::FftPeriodic:
                    return "FFT";
                case PoissonBackendKind::Open:
                    return "OPEN";
                case PoissonBackendKind::ConjugateGradient:
                    return "CG";
                case PoissonBackendKind::P3M:
                    return "P3M";
            }
            throw OpalException("Pic3DSolver::Pic3DSolver", "Unknown Poisson backend.");
        }

        IpplPoissonFields bindFields(Pic3DSolver::Workspace& workspace) {
            return {&workspace.chargeDensity(), &workspace.electricField()};
        }

        ParticleContainer& requirePrimaryParticles(
                std::span<const ParticleFieldBinding3d> bindings) {
            if (bindings.empty() || bindings.front().container == nullptr) {
                throw OpalException(
                        "Pic3DSolver::Pic3DSolver",
                        "The primary particle container is not available.");
            }
            return *bindings.front().container;
        }

    }  // namespace

    Pic3DSolver::Pic3DSolver(
            Pic3DConfig config, std::span<const ParticleFieldBinding3d> particleBindings,
            std::unique_ptr<Workspace> workspace, DataSink* dataSink)
        : primary_m(&requirePrimaryParticles(particleBindings)),
          workspace_m(std::move(workspace)),
          backendKind_m(config.backend()),
          dataSink_m(dataSink),
          binningConfig_m(config.binning()),
          correctionPlan_m(config.correction(), config.binning().has_value()),
          particleDomain_m(particleBindings),
          domainManager_m(config) {
        if (workspace_m == nullptr) {
            throw OpalException("Pic3DSolver::Pic3DSolver", "The PIC workspace is null.");
        }
        if (dataSink_m == nullptr) {
            throw OpalException("Pic3DSolver::Pic3DSolver", "The data sink is null.");
        }
        if (backendKind_m == PoissonBackendKind::ConjugateGradient) {
            throw OpalException(
                    "Pic3DSolver::Pic3DSolver",
                    "The CG Poisson backend is recognized but not implemented.");
        }
        if (binningConfig_m.has_value()) {
            binningPlan_m = std::make_unique<BinningPlan_t>(*primary_m, *binningConfig_m);
            binStats_m.reserve(binningPlan_m->maximumBinCount());
        }

        workspace_m->initializeFields(backendName(backendKind_m));
        IpplPoissonBackendConfig backendConfig;
        backendConfig.kind               = backendKind_m;
        backendConfig.greenFunction      = config.greenFunction();
        backendConfig.p3mCutoff          = config.p3mCutoff();
        backendConfig.boundaryConditions = config.boundaryConditions();
        backend_m = std::make_unique<IpplPoissonAdapter>(backendConfig, bindFields(*workspace_m));
        if (backendKind_m == PoissonBackendKind::P3M) {
            shortRangeInteraction_m.emplace(config.p3mCutoff());
        }
        backend_m->warmup();
    }

    SolverCapabilities Pic3DSolver::capabilities() const {
        SolverCapabilities result;
        result.particleSelection       = ParticleSelectionPolicy::PrimaryOnly;
        result.supportsBinning         = true;
        result.supportsImageCharge     = true;
        result.supportsShiftedGreen    = true;
        result.supportsRedistribution  = true;
        result.supportsPotentialOutput = true;
        return result;
    }

    void Pic3DSolver::execute(SolveContext& context, SelfFieldDiagnostics& diagnostics) {
        const PreparedCorrection correction =
                correctionPlan_m.prepare(context.requestedPhysics(), context.stepState().step);

        ParticleFrameGuard<double, 3> frameGuard(context.stepState().frames, *primary_m);
        try {
            frameGuard.enter();
            domainManager_m.update(
                    PicDomainFrame::Beam, context, *workspace_m, particleDomain_m, *backend_m,
                    diagnostics);

            frameGuard.markComputedFields();
            solveInBeamFrame(context, correction, diagnostics);

            frameGuard.leave();
            domainManager_m.update(
                    PicDomainFrame::Reference, context, *workspace_m, particleDomain_m, *backend_m,
                    diagnostics);
        } catch (...) {
            const std::exception_ptr originalException = std::current_exception();
            frameGuard.restoreNoThrow();
            if (!frameGuard.positionsInSolveFrame()) {
                try {
                    domainManager_m.update(
                            PicDomainFrame::Reference, context, *workspace_m, particleDomain_m,
                            *backend_m, diagnostics);
                } catch (...) {
                }
            }
            std::rethrow_exception(originalException);
        }
    }

    void Pic3DSolver::solveInBeamFrame(
            SolveContext& context, const PreparedCorrection& correction,
            SelfFieldDiagnostics& diagnostics) {
        Inform m("Pic3DSolver::solveInBeamFrame");
        const auto& backendCapabilities = backend_m->capabilities();
        if (backendCapabilities.isNoOp) {
            m << level5 << "Skipping scatter/gather and self-field computation for NONE solver."
              << endl;
            return;
        }

        if (primary_m->getTotalNum() <= 1) {
            Kokkos::deep_copy(primary_m->E.getView(), Vector_t<double, 3>(0.0));
            Kokkos::deep_copy(primary_m->B.getView(), Vector_t<double, 3>(0.0));
            return;
        }
        if (primary_m->getChargePerParticle() == 0.0) {
            throw OpalException(
                    "Pic3DSolver::solveInBeamFrame",
                    "Per-particle charge is zero but a self-field solver is active (type="
                            + std::string(backendName(backendKind_m))
                            + "). This almost always means the BEAM command is missing BCHARGE. "
                              "Set BCHARGE on the BEAM definition, or use TYPE=NONE when no "
                              "space charge is intended.");
        }
        if (correction.passCount == 0) {
            throw OpalException(
                    "Pic3DSolver::solveInBeamFrame",
                    "The prepared correction contains no primary solve pass.");
        }

        const bool binned = binningPlan_m != nullptr;
        m << level4 << "Entry: rank=" << ippl::Comm->rank()
          << ", localParticles=" << primary_m->getLocalNum()
          << ", totalParticles=" << primary_m->getTotalNum() << ", hasBins=" << (binned ? 1 : 0)
          << ", stype=" << backendName(backendKind_m) << endl;

        if (correction.correctionExpired) {
            m << level3 << "ZEROFACE_MAXSTEPS reached (step=" << context.stepState().step
              << ", maxSteps=" << correction.maximumSteps << "); disabling ";
            if (correction.configuredCorrection.kind == CorrectionKind::ShiftedGreen) {
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
            executeBinned(context, correction, diagnostics);
        } else {
            executeWholeBunch(context, correction, diagnostics);
        }
    }

    void Pic3DSolver::executeWholeBunch(
            SolveContext& context, const PreparedCorrection& correction,
            SelfFieldDiagnostics& diagnostics) {
        diagnostics.recordBinCount(0);
        bool dumpedDirichletPlaneThisStep = false;
        for (const SolvePass pass : correction) {
            executeSolvePass(context, nullptr, pass, dumpedDirichletPlaneThisStep, diagnostics);
        }
    }

    void Pic3DSolver::executeBinned(
            SolveContext& context, const PreparedCorrection& correction,
            SelfFieldDiagnostics& diagnostics) {
        const long long step = static_cast<long long>(context.stepState().step);
        const bool captureSnapshots =
                dataSink_m != nullptr && binningConfig_m.has_value()
                && !binningConfig_m->dumpFile().empty() && binningConfig_m->dumpFrequency() > 0
                && step % static_cast<long long>(binningConfig_m->dumpFrequency()) == 0;
        const BinningPreparation prepared = binningPlan_m->prepare(captureSnapshots);
        if (prepared.beforeMerge.has_value()) {
            dumpBinSnapshot(context, *prepared.beforeMerge, true);
        }
        if (prepared.afterMerge.has_value()) {
            dumpBinSnapshot(context, *prepared.afterMerge, false);
        }
        diagnostics.recordBinCount(prepared.mergedBinCount);
        fieldComposer_m.clearAccumulation(*workspace_m);

        Inform m("Pic3DSolver::executeBinned");
        m << level4 << "Iteration mode=binned, nBins=" << static_cast<int>(prepared.mergedBinCount)
          << ", stype=" << backendName(backendKind_m) << endl;

        binStats_m.clear();
        bool dumpedDirichletPlaneThisStep = false;
        while (const std::optional<BinnedSolveUnit_t> unit = binningPlan_m->next()) {
            m << level4 << "binIndex=" << static_cast<int>(unit->ordinal)
              << " nPartGlobal=" << static_cast<unsigned long long>(unit->globalParticleCount)
              << " gammaBin=" << std::setprecision(10) << unit->gamma << endl;
            binStats_m.push_back(
                    BinStatsRow{
                            static_cast<long long>(unit->ordinal),
                            static_cast<unsigned long long>(unit->globalParticleCount),
                            unit->gamma});

            for (const SolvePass pass : correction) {
                executeSolvePass(context, &*unit, pass, dumpedDirichletPlaneThisStep, diagnostics);
            }
        }

        fieldComposer_m.gatherAccumulated(
                scatterGather_m, primary_m->E, primary_m->B, primary_m->R, *workspace_m);
        if (diagnostics.shouldPrintBinTable(step)) {
            printBinStatsTable();
        }
    }

    void Pic3DSolver::executeSolvePass(
            SolveContext& context, const BinnedSolveUnit_t* unit, SolvePass pass,
            bool& dumpedDirichletPlaneThisStep, SelfFieldDiagnostics& diagnostics) {
        const double planeZ     = context.requestedPhysics().correction.planeZ;
        const PassPolicy policy = solvePassPolicy<double, 3>(pass, planeZ);
        if (policy.binned != (unit != nullptr)) {
            throw OpalException(
                    "Pic3DSolver::executeSolvePass",
                    "The solve-pass tag does not match the particle traversal.");
        }

        const auto& capabilities = backend_m->capabilities();
        if (unit == nullptr) {
            ScatterGather::ChargeNormalization normalization;
            normalization.timeStep              = context.stepState().timeStep;
            normalization.gamma                 = 1.0;
            normalization.selectedCharge        = primary_m->getTotalCharge();
            normalization.couplingConstant      = backend_m->couplingConstant();
            normalization.normalizeByCellVolume = capabilities.normalizeChargeByCellVolume;
            normalization.subtractNeutralizingBackground =
                    capabilities.subtractNeutralizingBackground;
            scatterGather_m.depositCharge(
                    *primary_m, *workspace_m, policy.depositKind,
                    ScatterGather::Selection::direct(0, primary_m->getLocalNum()), normalization,
                    policy.imagePolicy);
        } else {
            prepareRhoForUnit(context, *unit, policy);
        }

        workspace_m->electricField() = 0.0;
        auto& mesh                   = workspace_m->mesh();
        const auto originalSpacing   = mesh.getMeshSpacing();
        auto stretchedSpacing        = originalSpacing;
        const double gamma           = unit == nullptr ? 1.0 : unit->gamma;
        stretchedSpacing[2] *= gamma;

        IpplPoissonSolveRequest backendRequest;
        if (policy.backendRule == BackendSolveRule::ShiftedGreen) {
            const auto origin = mesh.getOrigin();
            const int longitudinalExtent =
                    static_cast<int>(workspace_m->layout().getDomain()[2].length());
            const double zCenter =
                    origin[2] + 0.5 * static_cast<double>(longitudinalExtent) * originalSpacing[2];
            ippl::Vector<double, 3> shift(0.0);
            shift[2]                          = 2.0 * gamma * (zCenter - planeZ);
            backendRequest.greenFunctionShift = shift;
        }

        bool spacingStretched = false;
        try {
            if (unit != nullptr) {
                mesh.setMeshSpacing(stretchedSpacing);
                spacingStretched = true;
            }

            Inform m("Pic3DSolver::executeSolvePass");
            m << level4 << "pass=" << solvePassLabel(pass)
              << ", binIndex=" << static_cast<int>(unit == nullptr ? 0 : unit->ordinal)
              << ", suppressFieldDump=" << (policy.suppressFieldDump ? 1 : 0);
            if (backendRequest.hasShiftedGreenFunction()) {
                m << ", plane=" << planeZ
                  << ", shift_z=" << (*backendRequest.greenFunctionShift)[2];
            }
            m << endl;

            backend_m->solve(backendRequest, {.suppressFieldDump = policy.suppressFieldDump});
            diagnostics.completeBackendSolve();

            if (unit == nullptr && policy.dumpDirichletPlaneAfter) {
                dumpDirichletPlaneDiagnosticsIfRequested(context, "legacy", planeZ, diagnostics);
                dumpedDirichletPlaneThisStep = true;
            }

            if (unit == nullptr) {
                fieldComposer_m.gatherElectrostatic(
                        scatterGather_m, primary_m->E, primary_m->R, *workspace_m);
                if (shortRangeInteraction_m.has_value()) {
                    shortRangeInteraction_m->apply(*primary_m);
                }
            } else {
                FieldComposer_t::Policy compositionPolicy;
                compositionPolicy.meanMomentum = unit->meanMomentum;
                compositionPolicy.gamma        = unit->gamma;
                compositionPolicy.magneticSign = policy.magneticSign;
                compositionPolicy.sourceRule   = policy.sourceRule;
                fieldComposer_m.accumulate(*workspace_m, compositionPolicy);
            }

            if (spacingStretched) {
                mesh.setMeshSpacing(originalSpacing);
                spacingStretched = false;
            }
        } catch (...) {
            if (spacingStretched) {
                mesh.setMeshSpacing(originalSpacing);
            }
            throw;
        }

        if (unit != nullptr && policy.dumpDirichletPlaneAfter && !dumpedDirichletPlaneThisStep) {
            dumpDirichletPlaneDiagnosticsIfRequested(context, "binned", planeZ, diagnostics);
            dumpedDirichletPlaneThisStep = true;
        }
    }

    void Pic3DSolver::prepareRhoForUnit(
            const SolveContext& context, const BinnedSolveUnit_t& unit, const PassPolicy& pass) {
        Inform m("Pic3DSolver::prepareRhoForUnit");
        const auto& indexedSelection = unit.indexedSelection;
        const auto& selectionPolicy  = indexedSelection.policy();
        const auto& hash             = indexedSelection.hash();
        const std::size_t begin      = static_cast<std::size_t>(selectionPolicy.begin());
        const std::size_t end        = static_cast<std::size_t>(selectionPolicy.end());
        const std::size_t hashExtent = static_cast<std::size_t>(hash.extent(0));
        const std::size_t localCount = primary_m->getLocalNum();

        if (end > hashExtent) {
            throw OpalException(
                    "Pic3DSolver::prepareRhoForUnit",
                    "Bin scatter policy exceeds its hash extent.");
        }
        if (unit.coversAllLocalParticles != (begin == 0 && end == localCount)) {
            throw OpalException(
                    "Pic3DSolver::prepareRhoForUnit",
                    "The solve-unit all-local flag does not match its indexed selection.");
        }

        const char* depositName = "PrimaryAndImage";
        if (pass.depositKind == ScatterGather::DepositKind::Primary) {
            depositName = "PrimaryOnly";
        } else if (pass.depositKind == ScatterGather::DepositKind::Image) {
            depositName = "ImageOnly";
        }
        m << level5 << "prepareRho: scatter mode=" << depositName
          << ", path=" << (unit.coversAllLocalParticles ? "all-local" : "subset")
          << ", localP=" << static_cast<unsigned long long>(localCount) << ", policy=[" << begin
          << "," << end << ")"
          << ", hashExtent=" << static_cast<unsigned long long>(hashExtent) << endl;

        const auto& capabilities = backend_m->capabilities();
        ScatterGather::ChargeNormalization normalization;
        normalization.timeStep = context.stepState().timeStep;
        normalization.gamma    = unit.gamma;
        normalization.selectedCharge =
                primary_m->getChargePerParticle() * static_cast<double>(unit.globalParticleCount);
        normalization.couplingConstant               = backend_m->couplingConstant();
        normalization.normalizeByCellVolume          = capabilities.normalizeChargeByCellVolume;
        normalization.subtractNeutralizingBackground = capabilities.subtractNeutralizingBackground;

        scatterGather_m.depositCharge(
                *primary_m, *workspace_m, pass.depositKind, unit.depositSelection(), normalization,
                pass.imagePolicy);
    }

    void Pic3DSolver::dumpBinSnapshot(
            const SolveContext& context, const BinConfigurationSnapshot& snapshot,
            bool beforeMerge) const {
        const std::vector<std::size_t> particleCounts(
                snapshot.particleCounts.begin(), snapshot.particleCounts.end());
        dataSink_m->dumpBinConfig(
                static_cast<long long>(context.stepState().step), context.stepState().time,
                beforeMerge, particleCounts, snapshot.widths, snapshot.lowerBound,
                binningConfig_m->dumpFile());
    }

    void Pic3DSolver::dumpDirichletPlaneDiagnosticsIfRequested(
            const SolveContext& context, const std::string& solveTag, double planeZ,
            SelfFieldDiagnostics& diagnostics) {
        const long long step = static_cast<long long>(context.stepState().step);
        if (!diagnostics.shouldDumpPlane(step)) {
            return;
        }
        Inform m("Pic3DSolver::dumpDirichletPlaneDiagnosticsIfRequested");
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
                step, context.stepState().time, planeZ, workspace_m->chargeDensity(), solveTag);
        if (planeDiagnostics.sampleCount == 0) {
            return;
        }
        m << level2 << "Dirichlet-plane potential diagnostics (" << solveTag << ") at step " << step
          << ": z=" << planeZ << " m, mean(phi)=" << planeDiagnostics.mean
          << " V, var(phi)=" << planeDiagnostics.variance << " V^2" << endl;
    }

    void Pic3DSolver::printBinStatsTable() const {
        const std::string& diagnosticName = binningPlan_m->diagnosticName();
        const std::string informName =
                diagnosticName.empty() ? "Pic3DSolver::printBinStatsTable"
                                       : "Pic3DSolver::printBinStatsTable[" + diagnosticName + "]";
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
