#include "Structure/DataSink.h"

#include <utility>
#include <vector>

template <typename T, unsigned Dim>
BinnedFieldSolver<T, Dim>::BinnedFieldSolver(
        std::string solver, std::shared_ptr<Workspace_t> workspace,
        std::shared_ptr<BCHandler_t> bcHandler, std::string greensFunction)
    : FieldSolver<T, Dim>(solver, std::move(workspace), bcHandler, std::move(greensFunction)) {
    scatterAttribute_m = ScatterAttribute::ChargeQ;
    gatherAttribute_m  = GatherAttribute::ElectricFieldE;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::configureIterationMetadata(
        bool binningConfigured, std::size_t maximumBinCount) {
    if (binningConfigured && maximumBinCount == 0) {
        throw OpalException(
                "BinnedFieldSolver::configureIterationMetadata",
                "Configured binning requires a positive maximum bin count.");
    }

    binningConfigured_m = binningConfigured;
    maximumBinCount_m   = binningConfigured ? maximumBinCount : 0;
    currentBinCount_m   = binningConfigured ? maximumBinCount : 1;
    if (binningConfigured) {
        binStats_m.reserve(maximumBinCount);
    }
}

template <typename T, unsigned Dim>
int BinnedFieldSolver<T, Dim>::legacyReportedBinCount() const {
    if (!binningConfigured_m) {
        return 1;
    }

    if (currentBinCount_m == maximumBinCount_m) {
        Inform m("PartBunch::getCurrentNBins");
        m << level4
          << "WARNING: Number of bins is the same as the maximum number of bins, we haven't "
             "merged bins yet (likely because the simulation is too empty). Returning 1. If "
             "that is not the case, check e.g. binning parameters."
          << endl;
        return 1;
    }
    return static_cast<int>(currentBinCount_m);
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::computeSelfFields(
        PartBunch_t& bunch, IterationPlan_t& iterationPlan, std::uint64_t particleGeneration,
        const PreparedCorrection_t& correction,
        opalx::spacecharge::BinConfigurationObserver* binObserver,
        opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
    // Validate inputs and decide between binned vs legacy solver.
    std::shared_ptr<ParticleCtr_t> pc = bunch.getParticleContainer();
    if (!pc) {
        throw OpalException(
                "BinnedFieldSolver::computeSelfFields",
                "Bunch particle container is not available.");
    }

    Inform m("BinnedFieldSolver::computeSelfFields");
    // TYPE=NONE is a true no-op: skip all binning/scatter/solve/gather work.
    if (this->getBackendCapabilities().isNoOp) {
        // Already called in ParallelTracker::resetFields()
        // pc->E = 0.0;
        // pc->B = 0.0;
        m << level5 << "Skipping scatter/gather and self-field computation for NONE solver."
          << endl;
        return;
    }

    // Trivial global case where self-field has no physical effect. This must be
    // based on the global particle count, because early emission can leave some
    // MPI ranks empty while another rank owns the single emitted particle.
    if (pc->getTotalNum() <= 1) {
        Kokkos::deep_copy(pc->E.getView(), Vector_t<T, Dim>(0.0));
        Kokkos::deep_copy(pc->B.getView(), Vector_t<T, Dim>(0.0));
        return;
    }

    // Fail fast on a zero per-particle charge. prepareRhoForUnit scatters
    // dt*Q via scaleDtByCharge / unscaleDtByCharge, which computes 0 / 0
    // when Q == 0 and silently poisons the per-particle dt attribute with
    // NaN. The first scatter then returns rho = 0 but leaves dt = NaN,
    // and any subsequent scatter in the same timestep (e.g. the shifted-
    // Green's correction pass) propagates NaN into rho -> E -> particles.
    // Almost always this means BCHARGE was omitted from the BEAM
    // definition in the input file.
    if (pc->getTotalNum() > 0 && pc->getChargePerParticle() == 0.0) {
        throw OpalException(
                "BinnedFieldSolver::computeSelfFields",
                "Per-particle charge is zero but a self-field solver is active "
                "(type=" + this->getStype()
                        + "). This almost always means the BEAM command in the "
                          "input file is missing BCHARGE (bunch charge, in [C]). "
                          "Set e.g. 'BCHARGE = 1e-9' on the BEAM definition, or "
                          "switch the field solver to TYPE=NONE if no space "
                          "charge is intended.");
    }

    const bool hasBins = iterationPlan.kind() == opalx::spacecharge::IterationKind::Binning;

    m << level4 << "Entry: rank=" << ippl::Comm->rank() << ", localParticles=" << pc->getLocalNum()
      << ", totalParticles=" << pc->getTotalNum() << ", hasBins=" << (hasBins ? 1 : 0)
      << ", stype=" << this->getStype() << endl;

    if (correction.passCount == 0) {
        throw OpalException(
                "BinnedFieldSolver::computeSelfFields",
                "The prepared correction contains no primary solve pass.");
    }
    if (correction.correctionExpired) {
        m << level3 << "ZEROFACE_MAXSTEPS reached (step=" << bunch.getGlobalTrackStep()
          << ", maxSteps=" << correction.maximumSteps << "); disabling ";
        if (correction.configuredCorrection.kind
            == opalx::spacecharge::CorrectionKind::ShiftedGreen) {
            m << "SHIFTED_GREENS_FUNCTION correction for this step." << endl;
        } else {
            m << "image charges for this step." << endl;
        }
    }

    const PreparedIteration_t prepared = iterationPlan.prepare(particleGeneration, binObserver);
    if (prepared.kind != iterationPlan.kind()) {
        throw OpalException(
                "BinnedFieldSolver::computeSelfFields",
                "The prepared iteration kind does not match its plan.");
    }
    diagnostics.recordBinCount(prepared.mergedBinCount);

    if (hasBins) {
        if (!binningConfigured_m || iterationPlan.maximumBinCount() != maximumBinCount_m) {
            throw OpalException(
                    "BinnedFieldSolver::computeSelfFields",
                    "Binning-plan metadata does not match the configured solver facade.");
        }
        currentBinCount_m = prepared.mergedBinCount;
    } else {
        if (binningConfigured_m || iterationPlan.maximumBinCount() != 0) {
            throw OpalException(
                    "BinnedFieldSolver::computeSelfFields",
                    "Whole-bunch plan metadata does not match the configured solver facade.");
        }
        if (correction.shiftedIgnored) {
            m << level3 << "SHIFTED_GREENS_FUNCTION is set but no binning is active; "
              << "the legacy path does not apply the Dirichlet correction." << endl;
        }
    }

    m << level4 << "Dispatching prepared iteration plan." << endl;
    executeIterationPlan(
            bunch, iterationPlan, prepared, particleGeneration, correction, diagnostics);
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::setScatterAttribute(const ScatterAttribute attr) {
    // store the scatter attribute selection.
    scatterAttribute_m = attr;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::setGatherAttribute(const GatherAttribute attr) {
    // store the gather attribute selection.
    gatherAttribute_m = attr;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::dumpDirichletPlaneDiagnosticsIfRequested(
        PartBunch_t& bunch, const std::string& solveTag, double planeZ,
        opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
    const long long step = bunch.getGlobalTrackStep();
    if (!diagnostics.shouldDumpPlane(step)) {
        return;
    }

    Inform m("BinnedFieldSolver::dumpDirichletPlaneDiagnosticsIfRequested");

    if (ippl::Comm->size() != 1) {
        if (!warnedPlaneDumpParallelUnsupported_m) {
            warnedPlaneDumpParallelUnsupported_m = true;
            m << level3 << "Dirichlet-plane diagnostics currently support only single-rank runs. "
              << "Skipping dump and statistics output." << endl;
        }
        return;
    }

    Field_t<Dim>* potentialField = this->getBackendCapabilities().usesSeparatePotentialField
                                           ? this->getPhi()
                                           : this->getRho();
    if (!potentialField) {
        return;
    }
    DataSink* dataSink = bunch.getDataSink();
    if (!dataSink) {
        return;
    }

    const auto planeDiagnostics =
            dataSink->dumpDirichletPlane(step, bunch.getT(), planeZ, *potentialField, solveTag);
    if (planeDiagnostics.sampleCount == 0) {
        return;
    }

    m << level2 << "Dirichlet-plane potential diagnostics (" << solveTag << ") at step " << step
      << ": z=" << planeZ << " m, mean(phi)=" << planeDiagnostics.mean
      << " V, var(phi)=" << planeDiagnostics.variance << " V^2" << endl;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::printBinStatsTable(
        const std::string& binningCmdName, const std::vector<BinStatsRow>& rows) {
    // print the table header (metadata + column names).
    const std::string informName =
            binningCmdName.empty()
                    ? "BinnedFieldSolver::printBinStatsTable"
                    : ("BinnedFieldSolver::printBinStatsTable[" + binningCmdName + "]");
    Inform m(informName.c_str());
    // m << level2 << tableName << " | nBins=" << static_cast<int>(nBinsOrZero)
    //   << " | stype=" << this->getStype() << endl;
    m << level2 << std::setw(9) << "bin"
      << " | " << std::setw(13) << "nParticles"
      << " | " << std::setw(15) << "gammaBin" << endl;

    for (const auto& r : rows) {
        m << level2 << std::setw(9) << r.binNumber << " | " << std::setw(13) << r.nParticles
          << " | " << std::setw(15) << std::setprecision(10) << r.gammaBin << endl;
    }
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::executeIterationPlan(
        PartBunch_t& bunch, IterationPlan_t& iterationPlan, const PreparedIteration_t& prepared,
        std::uint64_t particleGeneration, const PreparedCorrection_t& correction,
        opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
    Workspace_t& workspace = this->getWorkspace();

    const bool binned = prepared.kind == opalx::spacecharge::IterationKind::Binning;
    if (binned) {
        fieldComposer_m.clearAccumulation(workspace);
    }

    const std::size_t nBins = prepared.mergedBinCount;
    Inform m("BinnedFieldSolver::executeIterationPlan");
    m << level4 << "Iteration mode=" << (binned ? "binned" : "whole-bunch")
      << ", nBins=" << static_cast<int>(nBins) << ", stype=" << this->getStype() << endl;

    binStats_m.clear();

    bool dumpedDirichletPlaneThisStep = false;
    std::size_t emittedUnitCount      = 0;

    while (iterationPlan.hasNext(prepared, particleGeneration)) {
        auto solveUnitEvent = diagnostics.scopedEvent(
                opalx::spacecharge::SelfFieldEventKind::SolveUnit, binned ? "bin" : "whole-bunch");
        const std::optional<SolveUnit_t> nextUnit =
                iterationPlan.next(prepared, particleGeneration);
        if (!nextUnit.has_value()) {
            throw OpalException(
                    "BinnedFieldSolver::executeIterationPlan",
                    "An iteration plan reported a solve unit but did not emit it.");
        }
        const SolveUnit_t& unit = *nextUnit;
        ++emittedUnitCount;

        if (unit.kind != prepared.kind) {
            throw OpalException(
                    "BinnedFieldSolver::executeIterationPlan",
                    "The solve-unit kind does not match its prepared iteration.");
        }
        if (unit.fieldMode == opalx::spacecharge::SolveUnitFieldMode::Direct) {
            if (binned || emittedUnitCount != 1) {
                throw OpalException(
                        "BinnedFieldSolver::executeIterationPlan",
                        "A direct solve unit is invalid for this prepared iteration.");
            }
        } else if (
                !binned
                || unit.fieldMode != opalx::spacecharge::SolveUnitFieldMode::LorentzTransformed) {
            throw OpalException(
                    "BinnedFieldSolver::executeIterationPlan",
                    "A Lorentz-transformed solve unit requires a binning plan.");
        } else if (unit.ordinal >= nBins) {
            throw OpalException(
                    "BinnedFieldSolver::executeIterationPlan",
                    "A binning plan emitted a solve-unit ordinal outside its prepared range.");
        }

        if (binned) {
            if (unit.gamma <= 0.0) {
                throw OpalException(
                        "BinnedFieldSolver::executeIterationPlan",
                        "Computed non-positive gamma for bin.");
            }
            m << level4 << "binIndex=" << static_cast<int>(unit.ordinal)
              << " nPartGlobal=" << static_cast<unsigned long long>(unit.globalParticleCount)
              << " gammaBin=" << std::setprecision(10) << unit.gamma << endl;
            binStats_m.push_back(
                    BinStatsRow{
                            static_cast<long long>(unit.ordinal),
                            static_cast<unsigned long long>(unit.globalParticleCount), unit.gamma});
        }

        for (const SolvePass_t& pass : correction) {
            executeSolvePass(bunch, unit, pass, dumpedDirichletPlaneThisStep, diagnostics);
        }
    }

    if (!binned) {
        if (emittedUnitCount != 1) {
            throw OpalException(
                    "BinnedFieldSolver::executeIterationPlan",
                    "A whole-bunch iteration must emit exactly one direct solve unit.");
        }
        return;
    }

    // after all bins, gather the accumulated lab-frame field back to particles.
    {
        auto compositionEvent = diagnostics.scopedEvent(
                opalx::spacecharge::SelfFieldEventKind::FieldComposition, "final-gather");
        if (gatherAttribute_m != GatherAttribute::ElectricFieldE) {
            throw OpalException(
                    "BinnedFieldSolver::executeIterationPlan",
                    "Unsupported gather attribute in binned solver.");
        }

        std::shared_ptr<ParticleCtr_t> pc = bunch.getParticleContainer();
        fieldComposer_m.gatherAccumulated(scatterGather_m, pc->E, pc->B, pc->R, workspace);
    }

    // per-call table: gammaBin / nParticles / binNumber.
    if (diagnostics.shouldPrintBinTable(bunch.getGlobalTrackStep())) {
        printBinStatsTable(iterationPlan.diagnosticName(), binStats_m);
    }
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::executeSolvePass(
        PartBunch_t& bunch, const SolveUnit_t& unit, const SolvePass_t& pass,
        bool& dumpedDirichletPlaneThisStep, opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
    const bool binned =
            unit.fieldMode == opalx::spacecharge::SolveUnitFieldMode::LorentzTransformed;
    const bool direct = unit.fieldMode == opalx::spacecharge::SolveUnitFieldMode::Direct;
    if ((direct && pass.outputRule != opalx::spacecharge::FieldOutputRule::DirectGather)
        || (binned && pass.outputRule != opalx::spacecharge::FieldOutputRule::LorentzAccumulation)
        || (!direct && !binned)) {
        throw OpalException(
                "BinnedFieldSolver::executeSolvePass",
                "The solve-pass output rule does not match its solve unit.");
    }

    auto solvePassEvent = diagnostics.scopedEvent(
            opalx::spacecharge::SelfFieldEventKind::SolvePass, pass.solveLabel());

    std::shared_ptr<ParticleCtr_t> pc = bunch.getParticleContainer();
    if (!pc) {
        throw OpalException(
                "BinnedFieldSolver::executeSolvePass",
                "Bunch particle container is not available.");
    }
    if (scatterAttribute_m != ScatterAttribute::ChargeQ) {
        throw OpalException(
                "BinnedFieldSolver::executeSolvePass",
                "Unsupported scatter attribute in PIC solver.");
    }

    Workspace_t& workspace          = this->getWorkspace();
    const auto& backendCapabilities = this->getBackendCapabilities();
    if (direct) {
        typename ScatterGather_t::ChargeNormalization normalization;
        normalization.timeStep              = bunch.getdT();
        normalization.gamma                 = 1.0;
        normalization.selectedCharge        = pc->getTotalCharge();
        normalization.couplingConstant      = this->getCouplingConstant();
        normalization.normalizeByCellVolume = backendCapabilities.normalizeChargeByCellVolume;
        normalization.subtractNeutralizingBackground =
                backendCapabilities.subtractNeutralizingBackground;
        scatterGather_m.depositCharge(
                *pc, workspace, pass.depositKind, unit.depositSelection(), normalization,
                pass.imagePolicy);
    } else {
        prepareRhoForUnit(bunch, unit, pass);
    }

    // Ensure deterministic output even for solver types that do not update `E`.
    *(this->getE()) = 0.0;

    auto& mesh        = this->getRho()->get_mesh();
    const auto hrOrig = mesh.getMeshSpacing();
    auto hrStretched  = hrOrig;
    hrStretched[Dim - 1] *= unit.gamma;

    opalx::spacecharge::IpplPoissonSolveRequest backendRequest;
    if (pass.backendRule == opalx::spacecharge::BackendSolveRule::ShiftedGreen) {
        if (!binned) {
            throw OpalException(
                    "BinnedFieldSolver::executeSolvePass",
                    "A shifted Green solve requires a binned solve unit.");
        }
        const auto origin = mesh.getOrigin();
        const int longitudinalExtent =
                static_cast<int>(this->getRho()->getLayout().getDomain()[Dim - 1].length());
        const double zCenter =
                origin[Dim - 1] + 0.5 * static_cast<double>(longitudinalExtent) * hrOrig[Dim - 1];
        ippl::Vector<double, Dim> shift(0.0);
        shift[Dim - 1]                    = 2.0 * unit.gamma * (zCenter - pass.planeZ);
        backendRequest.greenFunctionShift = shift;
    }

    bool spacingStretched = false;
    try {
        if (binned) {
            mesh.setMeshSpacing(hrStretched);
            spacingStretched = true;
        }

        Inform m("BinnedFieldSolver::executeSolvePass");
        m << level4 << "pass=" << pass.solveLabel()
          << ", binIndex=" << static_cast<int>(unit.ordinal)
          << ", suppressFieldDump=" << (pass.suppressFieldDump ? 1 : 0);
        if (backendRequest.hasShiftedGreenFunction()) {
            m << ", plane=" << pass.planeZ
              << ", shift_z=" << (*backendRequest.greenFunctionShift)[Dim - 1];
        }
        m << endl;

        this->runSolver(backendRequest, pass.suppressFieldDump, diagnostics);

        if (direct && pass.dumpDirichletPlaneAfter) {
            dumpDirichletPlaneDiagnosticsIfRequested(bunch, "legacy", pass.planeZ, diagnostics);
            dumpedDirichletPlaneThisStep = true;
        }

        {
            auto compositionEvent = diagnostics.scopedEvent(
                    opalx::spacecharge::SelfFieldEventKind::FieldComposition,
                    pass.compositionLabel());
            if (direct) {
                if (gatherAttribute_m != GatherAttribute::ElectricFieldE) {
                    throw OpalException(
                            "BinnedFieldSolver::executeSolvePass",
                            "Unsupported gather attribute in whole-bunch solver.");
                }
                fieldComposer_m.gatherElectrostatic(scatterGather_m, pc->E, pc->R, workspace);
            } else {
                typename FieldComposer_t::Policy compositionPolicy;
                compositionPolicy.meanMomentum = unit.meanMomentum;
                compositionPolicy.gamma        = unit.gamma;
                compositionPolicy.magneticSign = pass.magneticSign;
                compositionPolicy.sourceRule   = pass.sourceRule;
                fieldComposer_m.accumulate(workspace, compositionPolicy);
            }
        }

        if (spacingStretched) {
            mesh.setMeshSpacing(hrOrig);
            spacingStretched = false;
        }
    } catch (...) {
        if (spacingStretched) {
            mesh.setMeshSpacing(hrOrig);
        }
        throw;
    }

    if (binned && pass.dumpDirichletPlaneAfter && !dumpedDirichletPlaneThisStep) {
        dumpDirichletPlaneDiagnosticsIfRequested(bunch, "binned", pass.planeZ, diagnostics);
        dumpedDirichletPlaneThisStep = true;
    }
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::prepareRhoForUnit(
        PartBunch_t& bunch, const SolveUnit_t& unit, const SolvePass_t& pass) {
    // Scatter bin charge to rho using the dt*Q deposition workflow.
    Inform m("BinnedFieldSolver::prepareRhoForUnit");
    const std::size_t binIndex    = unit.ordinal;
    const std::size_t nPartGlobal = unit.globalParticleCount;
    const double gammaBin         = unit.gamma;
    m << level4 << "prepareRho: binIndex=" << static_cast<int>(binIndex)
      << ", nPartGlobal=" << static_cast<unsigned long long>(nPartGlobal)
      << ", gammaBin=" << std::setprecision(10) << gammaBin << endl;

    // access particle views and validate scatter support.
    std::shared_ptr<ParticleCtr_t> pc = bunch.getParticleContainer();

    if (scatterAttribute_m != ScatterAttribute::ChargeQ) {
        throw OpalException(
                "BinnedFieldSolver::prepareRhoForUnit",
                "Unsupported scatter attribute in binned solver.");
    }

    const auto& indexedSelection = unit.indexedSelection;
    const auto& policy           = indexedSelection.policy();
    const auto& hash             = indexedSelection.hash();
    const std::size_t pBegin     = static_cast<std::size_t>(policy.begin());
    const std::size_t pEnd       = static_cast<std::size_t>(policy.end());
    const std::size_t hExtent    = static_cast<std::size_t>(hash.extent(0));
    const std::size_t nLocal     = pc->getLocalNum();
    const char* modeName         = "PrimaryAndImage";
    if (pass.depositKind == ScatterGather_t::DepositKind::Primary) {
        modeName = "PrimaryOnly";
    } else if (pass.depositKind == ScatterGather_t::DepositKind::Image) {
        modeName = "ImageOnly";
    }

    if (pEnd > hExtent) {
        throw OpalException(
                "BinnedFieldSolver::prepareRhoForUnit",
                "Bin scatter policy exceeds hash extent: policyEnd=" + std::to_string(pEnd)
                        + ", hashExtent=" + std::to_string(hExtent) + ".");
    }

    // If the selected bin spans every local particle, the hashed subset scatter is
    // equivalent to the all-local scatter. Prefer the all-local path here: it avoids
    // dereferencing the bin hash view in the hot scatter kernel and is the common
    // AWAGun early-emission case after 128 bins merge down to one bin.
    const bool scatterAllLocal = unit.coversAllLocalParticles;
    if (scatterAllLocal != (pBegin == 0 && pEnd == nLocal)) {
        throw OpalException(
                "BinnedFieldSolver::prepareRhoForUnit",
                "The solve-unit all-local flag does not match its indexed selection.");
    }
    m << level5 << "prepareRho: scatter mode=" << modeName
      << ", path=" << (scatterAllLocal ? "all-local" : "subset")
      << ", localP=" << static_cast<unsigned long long>(nLocal) << ", policy=[" << pBegin << ","
      << pEnd << "), hashExtent=" << static_cast<unsigned long long>(hExtent) << endl;

    const auto selection = unit.depositSelection();

    const auto& backendCapabilities = this->getBackendCapabilities();

    typename ScatterGather_t::ChargeNormalization normalization;
    normalization.timeStep         = bunch.getdT();
    normalization.gamma            = gammaBin;
    normalization.selectedCharge   = pc->getChargePerParticle() * static_cast<double>(nPartGlobal);
    normalization.couplingConstant = this->getCouplingConstant();
    normalization.normalizeByCellVolume = backendCapabilities.normalizeChargeByCellVolume;
    normalization.subtractNeutralizingBackground =
            backendCapabilities.subtractNeutralizingBackground;

    scatterGather_m.depositCharge(
            *pc, this->getWorkspace(), pass.depositKind, selection, normalization,
            pass.imagePolicy);
}
