/**
 * @file Pic3DSolver.cpp
 * @brief Implements the complete Cartesian 3D PIC self-field algorithm.
 */

#include "SpaceCharge/Pic3D/Pic3DSolver.h"

#include "SpaceCharge/Pic3D/IterationPlanFactory.h"
#include "Structure/DataSink.h"
#include "Utilities/OpalException.h"

#include <exception>
#include <iomanip>
#include <optional>
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
            }
            throw OpalException("Pic3DSolver::Pic3DSolver", "Unknown Poisson backend.");
        }

        IpplPoissonFields bindFields(Pic3DSolver::Workspace& workspace) {
            return {&workspace.chargeDensity(), &workspace.electricField(), &workspace.potential()};
        }

        ParticleContainer& requirePrimaryParticles(
                std::span<const std::shared_ptr<ParticleContainer>> particles) {
            if (particles.empty() || particles.front() == nullptr) {
                throw OpalException(
                        "Pic3DSolver::Pic3DSolver",
                        "The primary particle container is not available.");
            }
            return *particles.front();
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

    Pic3DSolver::Pic3DSolver(
            Pic3DConfig config,
            std::span<const std::shared_ptr<ParticleContainer>> particleContainers,
            std::shared_ptr<Workspace> workspace, DataSink* dataSink)
        : primary_m(&requirePrimaryParticles(particleContainers)),
          workspace_m(std::move(workspace)),
          backendKind_m(config.backend()),
          dataSink_m(dataSink),
          binningConfig_m(config.binning()),
          iterationPlan_m(IterationPlanFactory<double, 3>::create(binningConfig_m, *primary_m)),
          correctionPlan_m(config.correction(), config.backend(), iterationPlan_m->kind()),
          particleDomain_m(particleContainers),
          domainManager_m(config) {
        if (workspace_m == nullptr) {
            throw OpalException("Pic3DSolver::Pic3DSolver", "The PIC workspace is null.");
        }
        if (dataSink_m == nullptr) {
            throw OpalException("Pic3DSolver::Pic3DSolver", "The data sink is null.");
        }

        workspace_m->initializeFields(backendName(backendKind_m));
        backend_m = IpplPoissonAdapter::create(
                backendKind_m, config.greenFunction(), bindFields(*workspace_m));
        backend_m->setPotentialBoundaryConditions(config.boundaryConditions());
        backend_m->warmup();

        if (iterationPlan_m->kind() == IterationKind::Binning) {
            binStats_m.reserve(iterationPlan_m->maximumBinCount());
        }
    }

    SolverCapabilities Pic3DSolver::capabilities() const {
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

    void Pic3DSolver::execute(SolveContext& context, SelfFieldDiagnostics& diagnostics) {
        framePolicy_m.validate(context);
        const PreparedCorrection_t correction =
                correctionPlan_m.prepare(context.requestedPhysics(), context.stepState().step);
        context.particles().updateGeneration(particleDomain_m.generation());

        FieldFramePolicy::State frameState;
        try {
            framePolicy_m.enter(context, *primary_m, frameState);
            domainManager_m.update(
                    PicDomainFrame::Beam, context, *workspace_m, particleDomain_m, *backend_m,
                    diagnostics);

            // Mark both results before work begins so exception cleanup rotates any partially
            // written fields back with freshly acquired particle views.
            framePolicy_m.markComputedFieldsInBeam(frameState);
            solveInBeamFrame(context, correction, diagnostics);

            framePolicy_m.leave(context, *primary_m, frameState);
            domainManager_m.update(
                    PicDomainFrame::Reference, context, *workspace_m, particleDomain_m, *backend_m,
                    diagnostics);
        } catch (...) {
            const std::exception_ptr originalException = std::current_exception();
            framePolicy_m.restore(context, *primary_m, frameState);
            if (!frameState.positionsInBeam) {
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
            SolveContext& context, const PreparedCorrection_t& correction,
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

        const bool binned = iterationPlan_m->kind() == IterationKind::Binning;
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

        std::optional<DataSinkBinConfigurationObserver> binObserver;
        if (binningConfig_m.has_value()) {
            binObserver.emplace(
                    dataSink_m, static_cast<long long>(context.stepState().step),
                    context.stepState().time, *binningConfig_m);
        }
        const PreparedIteration prepared = iterationPlan_m->prepare(
                context.particles().generation(),
                binObserver.has_value() ? &*binObserver : nullptr);
        if (prepared.kind != iterationPlan_m->kind()) {
            throw OpalException(
                    "Pic3DSolver::solveInBeamFrame",
                    "The prepared iteration kind does not match its plan.");
        }
        diagnostics.recordBinCount(prepared.mergedBinCount);

        if (!binned && correction.shiftedIgnored) {
            m << level3 << "SHIFTED_GREENS_FUNCTION is set but no binning is active; "
              << "the whole-bunch path does not apply the Dirichlet correction." << endl;
        }

        executeIterationPlan(context, prepared, correction, diagnostics);
    }

    void Pic3DSolver::executeIterationPlan(
            SolveContext& context, const PreparedIteration& prepared,
            const PreparedCorrection_t& correction, SelfFieldDiagnostics& diagnostics) {
        const bool binned = prepared.kind == IterationKind::Binning;
        if (binned) {
            fieldComposer_m.clearAccumulation(*workspace_m);
        }

        Inform m("Pic3DSolver::executeIterationPlan");
        m << level4 << "Iteration mode=" << (binned ? "binned" : "whole-bunch")
          << ", nBins=" << static_cast<int>(prepared.mergedBinCount)
          << ", stype=" << backendName(backendKind_m) << endl;

        binStats_m.clear();
        bool dumpedDirichletPlaneThisStep = false;
        std::size_t emittedUnitCount      = 0;
        const auto generation             = context.particles().generation();

        while (iterationPlan_m->hasNext(prepared, generation)) {
            auto solveUnitEvent = diagnostics.scopedEvent(
                    SelfFieldEventKind::SolveUnit, binned ? "bin" : "whole-bunch");
            const std::optional<SolveUnit_t> nextUnit = iterationPlan_m->next(prepared, generation);
            if (!nextUnit.has_value()) {
                throw OpalException(
                        "Pic3DSolver::executeIterationPlan",
                        "An iteration plan reported a solve unit but did not emit it.");
            }
            const SolveUnit_t& unit = *nextUnit;
            ++emittedUnitCount;

            if (unit.kind != prepared.kind) {
                throw OpalException(
                        "Pic3DSolver::executeIterationPlan",
                        "The solve-unit kind does not match its prepared iteration.");
            }
            if (unit.fieldMode == SolveUnitFieldMode::Direct) {
                if (binned || emittedUnitCount != 1) {
                    throw OpalException(
                            "Pic3DSolver::executeIterationPlan",
                            "A direct solve unit is invalid for this prepared iteration.");
                }
            } else if (!binned || unit.fieldMode != SolveUnitFieldMode::LorentzTransformed) {
                throw OpalException(
                        "Pic3DSolver::executeIterationPlan",
                        "A Lorentz-transformed solve unit requires a binning plan.");
            } else if (unit.ordinal >= prepared.mergedBinCount) {
                throw OpalException(
                        "Pic3DSolver::executeIterationPlan",
                        "A binning plan emitted a solve-unit ordinal outside its prepared range.");
            }

            if (binned) {
                if (unit.gamma <= 0.0) {
                    throw OpalException(
                            "Pic3DSolver::executeIterationPlan",
                            "Computed non-positive gamma for bin.");
                }
                m << level4 << "binIndex=" << static_cast<int>(unit.ordinal)
                  << " nPartGlobal=" << static_cast<unsigned long long>(unit.globalParticleCount)
                  << " gammaBin=" << std::setprecision(10) << unit.gamma << endl;
                binStats_m.push_back(
                        BinStatsRow{
                                static_cast<long long>(unit.ordinal),
                                static_cast<unsigned long long>(unit.globalParticleCount),
                                unit.gamma});
            }

            for (const SolvePass_t& pass : correction) {
                executeSolvePass(context, unit, pass, dumpedDirichletPlaneThisStep, diagnostics);
            }
        }

        if (!binned) {
            if (emittedUnitCount != 1) {
                throw OpalException(
                        "Pic3DSolver::executeIterationPlan",
                        "A whole-bunch iteration must emit exactly one direct solve unit.");
            }
            return;
        }

        {
            auto compositionEvent =
                    diagnostics.scopedEvent(SelfFieldEventKind::FieldComposition, "final-gather");
            fieldComposer_m.gatherAccumulated(
                    scatterGather_m, primary_m->E, primary_m->B, primary_m->R, *workspace_m);
        }

        if (diagnostics.shouldPrintBinTable(static_cast<long long>(context.stepState().step))) {
            printBinStatsTable();
        }
    }

    void Pic3DSolver::executeSolvePass(
            SolveContext& context, const SolveUnit_t& unit, const SolvePass_t& pass,
            bool& dumpedDirichletPlaneThisStep, SelfFieldDiagnostics& diagnostics) {
        const bool binned = unit.fieldMode == SolveUnitFieldMode::LorentzTransformed;
        const bool direct = unit.fieldMode == SolveUnitFieldMode::Direct;
        if ((direct && pass.outputRule != FieldOutputRule::DirectGather)
            || (binned && pass.outputRule != FieldOutputRule::LorentzAccumulation)
            || (!direct && !binned)) {
            throw OpalException(
                    "Pic3DSolver::executeSolvePass",
                    "The solve-pass output rule does not match its solve unit.");
        }

        auto solvePassEvent =
                diagnostics.scopedEvent(SelfFieldEventKind::SolvePass, pass.solveLabel());
        const auto& capabilities = backend_m->capabilities();
        if (direct) {
            ScatterGather::ChargeNormalization normalization;
            normalization.timeStep              = context.stepState().timeStep;
            normalization.gamma                 = 1.0;
            normalization.selectedCharge        = primary_m->getTotalCharge();
            normalization.couplingConstant      = backend_m->couplingConstant();
            normalization.normalizeByCellVolume = capabilities.normalizeChargeByCellVolume;
            normalization.subtractNeutralizingBackground =
                    capabilities.subtractNeutralizingBackground;
            scatterGather_m.depositCharge(
                    *primary_m, *workspace_m, pass.depositKind, unit.depositSelection(),
                    normalization, pass.imagePolicy);
        } else {
            prepareRhoForUnit(context, unit, pass);
        }

        workspace_m->electricField() = 0.0;
        auto& mesh                   = workspace_m->mesh();
        const auto originalSpacing   = mesh.getMeshSpacing();
        auto stretchedSpacing        = originalSpacing;
        stretchedSpacing[2] *= unit.gamma;

        IpplPoissonSolveRequest backendRequest;
        if (pass.backendRule == BackendSolveRule::ShiftedGreen) {
            if (!binned) {
                throw OpalException(
                        "Pic3DSolver::executeSolvePass",
                        "A shifted Green solve requires a binned solve unit.");
            }
            const auto origin = mesh.getOrigin();
            const int longitudinalExtent =
                    static_cast<int>(workspace_m->layout().getDomain()[2].length());
            const double zCenter =
                    origin[2] + 0.5 * static_cast<double>(longitudinalExtent) * originalSpacing[2];
            ippl::Vector<double, 3> shift(0.0);
            shift[2]                          = 2.0 * unit.gamma * (zCenter - pass.planeZ);
            backendRequest.greenFunctionShift = shift;
        }

        bool spacingStretched = false;
        try {
            if (binned) {
                mesh.setMeshSpacing(stretchedSpacing);
                spacingStretched = true;
            }

            Inform m("Pic3DSolver::executeSolvePass");
            m << level4 << "pass=" << pass.solveLabel()
              << ", binIndex=" << static_cast<int>(unit.ordinal)
              << ", suppressFieldDump=" << (pass.suppressFieldDump ? 1 : 0);
            if (backendRequest.hasShiftedGreenFunction()) {
                m << ", plane=" << pass.planeZ
                  << ", shift_z=" << (*backendRequest.greenFunctionShift)[2];
            }
            m << endl;

            backend_m->solve(
                    backendRequest,
                    {.suppressFieldDump = pass.suppressFieldDump, .diagnostics = &diagnostics});

            if (direct && pass.dumpDirichletPlaneAfter) {
                dumpDirichletPlaneDiagnosticsIfRequested(
                        context, "legacy", pass.planeZ, diagnostics);
                dumpedDirichletPlaneThisStep = true;
            }

            {
                auto compositionEvent = diagnostics.scopedEvent(
                        SelfFieldEventKind::FieldComposition, pass.compositionLabel());
                if (direct) {
                    fieldComposer_m.gatherElectrostatic(
                            scatterGather_m, primary_m->E, primary_m->R, *workspace_m);
                } else {
                    FieldComposer_t::Policy compositionPolicy;
                    compositionPolicy.meanMomentum = unit.meanMomentum;
                    compositionPolicy.gamma        = unit.gamma;
                    compositionPolicy.magneticSign = pass.magneticSign;
                    compositionPolicy.sourceRule   = pass.sourceRule;
                    fieldComposer_m.accumulate(*workspace_m, compositionPolicy);
                }
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

        if (binned && pass.dumpDirichletPlaneAfter && !dumpedDirichletPlaneThisStep) {
            dumpDirichletPlaneDiagnosticsIfRequested(context, "binned", pass.planeZ, diagnostics);
            dumpedDirichletPlaneThisStep = true;
        }
    }

    void Pic3DSolver::prepareRhoForUnit(
            const SolveContext& context, const SolveUnit_t& unit, const SolvePass_t& pass) {
        Inform m("Pic3DSolver::prepareRhoForUnit");
        const auto& indexedSelection = unit.indexedSelection;
        const auto& policy           = indexedSelection.policy();
        const auto& hash             = indexedSelection.hash();
        const std::size_t begin      = static_cast<std::size_t>(policy.begin());
        const std::size_t end        = static_cast<std::size_t>(policy.end());
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

        const auto& capabilities               = backend_m->capabilities();
        Workspace::ScalarField* potentialField = capabilities.usesSeparatePotentialField
                                                         ? &workspace_m->potential()
                                                         : &workspace_m->chargeDensity();
        const auto planeDiagnostics            = dataSink_m->dumpDirichletPlane(
                step, context.stepState().time, planeZ, *potentialField, solveTag);
        if (planeDiagnostics.sampleCount == 0) {
            return;
        }
        m << level2 << "Dirichlet-plane potential diagnostics (" << solveTag << ") at step " << step
          << ": z=" << planeZ << " m, mean(phi)=" << planeDiagnostics.mean
          << " V, var(phi)=" << planeDiagnostics.variance << " V^2" << endl;
    }

    void Pic3DSolver::printBinStatsTable() const {
        const std::string& diagnosticName = iterationPlan_m->diagnosticName();
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
