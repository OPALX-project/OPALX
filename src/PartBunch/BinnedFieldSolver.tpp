#include "Structure/DataSink.h"

#include <cstring>
#include <utility>
#include <vector>

template <typename T, unsigned Dim>
BinnedFieldSolver<T, Dim>::BinnedFieldSolver(
        std::string solver, std::shared_ptr<Workspace_t> workspace,
        std::shared_ptr<BCHandler_t> bcHandler, int tablePrintFrequency, bool adaptiveBinning,
        std::string greensFunction)
    : FieldSolver<T, Dim>(solver, std::move(workspace), bcHandler, std::move(greensFunction)) {
    scatterAttribute_m = ScatterAttribute::ChargeQ;
    gatherAttribute_m  = GatherAttribute::ElectricFieldE;
    adaptiveBinning_m  = adaptiveBinning;
    static_cast<void>(tablePrintFrequency);
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::computeSelfFields(
        PartBunch_t& bunch, opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
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

    // Fail fast on a zero per-particle charge. prepareRhoForBin scatters
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

    // decide which solver path to run (binned vs legacy).
    const bool hasBins = bunch.hasBinning();

    m << level4 << "Entry: rank=" << ippl::Comm->rank() << ", localParticles=" << pc->getLocalNum()
      << ", totalParticles=" << pc->getTotalNum() << ", hasBins=" << (hasBins ? 1 : 0)
      << ", stype=" << this->getStype() << endl;

    // Temporarily disable image charges if the step limit has been reached.
    // The controller's enabled flag gates the image scatter pass; by disabling it
    // here both the legacy and binned paths automatically perform primary-only scatter.
    const bool imageWasEnabled     = imageScatterController_m.isEnabled();
    const bool imageActiveThisStep = isImageChargeActiveForStep(bunch.getGlobalTrackStep());

    // Mirror the same step-budget toggling for the shifted Green's path.
    const bool shiftedGreensWasEnabled = shiftedGreensEnabled_m;
    const bool shiftedGreensActiveThisStep =
            isShiftedGreensActiveForStep(bunch.getGlobalTrackStep());

    const auto restoreCorrectionState = [&]() {
        if (imageWasEnabled && !imageActiveThisStep) {
            imageScatterController_m.configure(true, imageScatterController_m.getZPlane());
        }
        if (shiftedGreensWasEnabled && !shiftedGreensActiveThisStep) {
            shiftedGreensEnabled_m = true;
        }
    };

    try {
        if (imageWasEnabled && !imageActiveThisStep) {
            m << level3 << "ZEROFACE_MAXSTEPS reached (step=" << bunch.getGlobalTrackStep()
              << ", maxSteps=" << zerofaceMaxSteps_m << "); disabling image charges for this step."
              << endl;
            imageScatterController_m.configure(false, imageScatterController_m.getZPlane());
        }
        if (shiftedGreensWasEnabled && !shiftedGreensActiveThisStep) {
            m << level3 << "ZEROFACE_MAXSTEPS reached (step=" << bunch.getGlobalTrackStep()
              << ", maxSteps=" << zerofaceMaxSteps_m
              << "); disabling SHIFTED_GREENS_FUNCTION correction for this step." << endl;
            shiftedGreensEnabled_m = false;
        }

        if (hasBins) {
            m << level4 << "Dispatching to computeBinnedSelfFields() (binned path)." << endl;
            computeBinnedSelfFields(bunch, diagnostics);
        } else {
            // Legacy path has no separate correction pass: PicScatterGather deposits
            // primary and image charge together before one standard solve. The shifted Green's
            // correction is only
            // implemented for the binned path, so warn once if the user requested it
            // without binning.
            if (shiftedGreensWasEnabled && shiftedGreensActiveThisStep) {
                m << level3 << "SHIFTED_GREENS_FUNCTION is set but no binning is active; "
                  << "the legacy path does not apply the Dirichlet correction." << endl;
            }
            m << level4 << "Dispatching to computeLegacySelfFields() (legacy path)." << endl;
            computeLegacySelfFields(bunch, diagnostics);
        }
    } catch (...) {
        try {
            restoreCorrectionState();
        } catch (...) {
        }
        throw;
    }

    restoreCorrectionState();
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
void BinnedFieldSolver<T, Dim>::setImageChargeConfiguration(bool enabled, double zPlane) {
    if (enabled && shiftedGreensEnabled_m) {
        throw OpalException(
                "BinnedFieldSolver::setImageChargeConfiguration",
                "Cannot enable image charges while SHIFTED_GREENS_FUNCTION is active: "
                "ZEROFACE_R0Z and SHIFTED_GREENS_FUNCTION are mutually exclusive.");
    }
    imageScatterController_m.configure(enabled, zPlane);
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::setShiftedGreensConfiguration(bool enabled, double zPlane) {
    if (enabled && imageScatterController_m.isEnabled()) {
        throw OpalException(
                "BinnedFieldSolver::setShiftedGreensConfiguration",
                "Cannot enable SHIFTED_GREENS_FUNCTION while image charges are active: "
                "ZEROFACE_R0Z and SHIFTED_GREENS_FUNCTION are mutually exclusive.");
    }
    shiftedGreensEnabled_m = enabled;
    shiftedGreensPlaneZ_m  = zPlane;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::setZerofaceMaxSteps(int maxSteps) {
    if (maxSteps < 0) {
        throw OpalException(
                "BinnedFieldSolver::setZerofaceMaxSteps", "ZEROFACE_MAXSTEPS must be >= 0.");
    }
    zerofaceMaxSteps_m = maxSteps;
}

template <typename T, unsigned Dim>
bool BinnedFieldSolver<T, Dim>::isImageChargeActiveForStep(size_t step) const {
    if (!imageScatterController_m.isEnabled()) {
        return false;
    }
    if (zerofaceMaxSteps_m <= 0) {
        return true;  // unlimited
    }
    return step < static_cast<size_t>(zerofaceMaxSteps_m);
}

template <typename T, unsigned Dim>
bool BinnedFieldSolver<T, Dim>::isShiftedGreensActiveForStep(size_t step) const {
    if (!shiftedGreensEnabled_m) {
        return false;
    }
    if (zerofaceMaxSteps_m <= 0) {
        return true;  // unlimited
    }
    return step < static_cast<size_t>(zerofaceMaxSteps_m);
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::setZeroFacePlaneDumpFrequency(int frequency) {
    if (frequency < 0) {
        throw OpalException(
                "BinnedFieldSolver::setZeroFacePlaneDumpFrequency",
                "ZEROFACEPLANEDUMP frequency must be >= 0.");
    }
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::dumpDirichletPlaneDiagnosticsIfRequested(
        PartBunch_t& bunch, const std::string& solveTag,
        opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
    if (!imageScatterController_m.isEnabled()) {
        return;
    }

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
    const double zPlane = imageScatterController_m.getZPlane();

    DataSink* dataSink = bunch.getDataSink();
    if (!dataSink) {
        return;
    }

    const auto planeDiagnostics =
            dataSink->dumpDirichletPlane(step, bunch.getT(), zPlane, *potentialField, solveTag);
    if (planeDiagnostics.sampleCount == 0) {
        return;
    }

    m << level2 << "Dirichlet-plane potential diagnostics (" << solveTag << ") at step " << step
      << ": z=" << zPlane << " m, mean(phi)=" << planeDiagnostics.mean
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
void BinnedFieldSolver<T, Dim>::computeBinnedSelfFields(
        PartBunch_t& bunch, opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
    // execute full binned self-field algorithm.
    // fetch the adaptive bin structure.
    std::shared_ptr<AdaptBins_t> bins = bunch.getBins();
    if (!bins) {
        // Defensive: runtime selection above should prevent this.
        computeLegacySelfFields(bunch, diagnostics);
        return;
    }

    // build and merge adaptive bins for this step.
    rebinAndPrepare(bunch, bins);

    // obtain the temporary buffers used to accumulate bin contributions.
    Workspace_t& workspace = this->getWorkspace();

    fieldComposer_m.clearAccumulation(workspace);

    // determine the number of bins used for this step.
    const bin_index_type nBins = bins->getCurrentBinCount();
    diagnostics.recordBinCount(static_cast<std::size_t>(nBins));

    // Level-5 debug: per-step overview before entering the bin loop.
    Inform m("BinnedFieldSolver::computeBinnedSelfFields");
    m << level4 << "Binned mode: nBins=" << static_cast<int>(nBins)
      << ", stype=" << this->getStype() << endl;

    // Cache values for the level-3 per-call table.
    std::vector<BinStatsRow> binStats;
    binStats.reserve(static_cast<size_t>(nBins));

    bool dumpedDirichletPlaneThisStep = false;

    // iterate over merged bins and accumulate E contributions.
    for (bin_index_type binIndex = 0; binIndex < nBins; ++binIndex) {
        // process a single merged bin (gamma->rho->solve->accumulate).
        const size_type nPartGlobal = bins->getNPartInBin(binIndex, true);
        if (nPartGlobal == 0) {
            continue;
        }

        auto solveUnitEvent =
                diagnostics.scopedEvent(opalx::spacecharge::SelfFieldEventKind::SolveUnit, "bin");

        m << level4 << "binIndex=" << static_cast<int>(binIndex)
          << " nPartGlobal=" << static_cast<unsigned long long>(nPartGlobal) << endl;

        // compute global average gamma for this bin.
        const BinKinematics kinematics = computeGammaBinGlobal(bunch, bins, binIndex, nPartGlobal);
        const double gammaBin          = kinematics.gammaBin;
        if (gammaBin <= 0.0) {
            throw OpalException(
                    "BinnedFieldSolver::computeBinnedSelfFields",
                    "Computed non-positive gamma for bin.");
        }

        m << level4 << "binIndex=" << static_cast<int>(binIndex)
          << " gammaBin=" << std::setprecision(10) << gammaBin << endl;

        binStats.push_back(
                BinStatsRow{
                        static_cast<long long>(binIndex),
                        static_cast<unsigned long long>(nPartGlobal), gammaBin});

        // Mesh references reused for both primary and correction passes.
        auto& mesh        = this->getRho()->get_mesh();
        const auto hrOrig = mesh.getMeshSpacing();
        auto hrStretched  = hrOrig;
        hrStretched[Dim - 1] *= gammaBin;

        const bool imageActive         = imageScatterController_m.isEnabled();
        const bool shiftedGreensActive = shiftedGreensEnabled_m;
        // Mutual exclusion enforced at config time; defensive assert.
        assert(!(imageActive && shiftedGreensActive));
        const bool correctionActive = imageActive || shiftedGreensActive;

        // --- Primary pass: scatter real charges, solve, accumulate with +B ---
        {
            auto solvePassEvent = diagnostics.scopedEvent(
                    opalx::spacecharge::SelfFieldEventKind::SolvePass, "primary");
            const ImageScatterMode scatterMode = correctionActive
                                                         ? ImageScatterMode::PrimaryOnly
                                                         : ImageScatterMode::PrimaryAndImage;
            prepareRhoForBin(bunch, bins, binIndex, nPartGlobal, gammaBin, scatterMode);

            *(this->getE()) = 0.0;
            mesh.setMeshSpacing(hrStretched);

            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " primary runSolver(true) start"
              << " (hr_z stretched by gamma=" << gammaBin << ")" << endl;
            this->runSolver(true, diagnostics);
            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " primary runSolver(true) done; accumulate->Etmp" << endl;

            {
                auto compositionEvent = diagnostics.scopedEvent(
                        opalx::spacecharge::SelfFieldEventKind::FieldComposition,
                        "primary-accumulation");
                typename FieldComposer_t::Policy compositionPolicy;
                for (unsigned d = 0; d < Dim; ++d) {
                    compositionPolicy.meanMomentum[d] = kinematics.pmean[d];
                }
                compositionPolicy.gamma        = gammaBin;
                compositionPolicy.magneticSign = +1.0;
                compositionPolicy.sourceRule   = opalx::spacecharge::FieldSourceRule::Direct;
                fieldComposer_m.accumulate(workspace, compositionPolicy);
            }

            mesh.setMeshSpacing(hrOrig);
        }

        // --- Dirichlet correction pass: one of two mutually-exclusive paths ---
        //
        // Path A (image charges): legacy path. Scatters mirrored particles onto
        // the same mesh, then solves with the standard Green's function. Only
        // correct while the mesh straddles the cathode plane; degrades silently
        // once the bunch has drifted beyond ZEROFACE_MAXSTEPS.
        //
        // Path B (shifted Green's function): new path. Re-scatters the primary
        // charges, then solves with a translated free-space kernel that encodes
        // the image-charge contribution analytically. The component-wise z-flip
        // + sign-flip on the solver output, applied by FieldComposer,
        // produces the image-charge E field directly. Works at any bunch-to-
        // plane distance and requires the OPEN solver (checked in
        // FieldSolver::runShiftedOpenSolver).
        if (imageActive) {
            auto solvePassEvent = diagnostics.scopedEvent(
                    opalx::spacecharge::SelfFieldEventKind::SolvePass, "image");
            prepareRhoForBin(
                    bunch, bins, binIndex, nPartGlobal, gammaBin, ImageScatterMode::ImageOnly);

            *(this->getE()) = 0.0;
            mesh.setMeshSpacing(hrStretched);

            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " image runSolver(true) start" << endl;
            this->runSolver(true, diagnostics);
            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " image runSolver(true) done; accumulate->Etmp (B negated)" << endl;

            {
                auto compositionEvent = diagnostics.scopedEvent(
                        opalx::spacecharge::SelfFieldEventKind::FieldComposition,
                        "image-accumulation");
                typename FieldComposer_t::Policy compositionPolicy;
                for (unsigned d = 0; d < Dim; ++d) {
                    compositionPolicy.meanMomentum[d] = kinematics.pmean[d];
                }
                compositionPolicy.gamma        = gammaBin;
                compositionPolicy.magneticSign = -1.0;
                compositionPolicy.sourceRule   = opalx::spacecharge::FieldSourceRule::Direct;
                fieldComposer_m.accumulate(workspace, compositionPolicy);
            }
            mesh.setMeshSpacing(hrOrig);

            // Dump phi ~= 0 check on the Dirichlet plane AFTER the correction
            // lands on the mesh. Only safe for the explicit-image path: there
            // the cathode plane always lies inside the computational domain.
            // For the shifted path the domain may be far from z = R0Z, so the
            // interpolated phi on the plane would be meaningless.
            if (!dumpedDirichletPlaneThisStep) {
                dumpDirichletPlaneDiagnosticsIfRequested(bunch, "binned", diagnostics);
                dumpedDirichletPlaneThisStep = true;
            }
        } else if (shiftedGreensActive) {
            auto solvePassEvent = diagnostics.scopedEvent(
                    opalx::spacecharge::SelfFieldEventKind::SolvePass, "shifted-green");
            // Shifted-Green's-function image correction. Multi-rank enabled: the
            // axis-flip source read crosses ranks under PARFFTZ, but that is handled
            // inside FieldComposer (it stages a globally mirrored field in a
            // local view populated via peer-rank MPI exchange).
            //
            // Re-scatter the primary charges (solve() overwrote the RHS in the
            // primary pass). This matches the ImageOnly path's pattern.
            prepareRhoForBin(
                    bunch, bins, binIndex, nPartGlobal, gammaBin, ImageScatterMode::PrimaryOnly);

            *(this->getE()) = 0.0;
            mesh.setMeshSpacing(hrStretched);

            // Shift formula in the bin rest frame. Old OPAL computes the same
            // distance as zshift = -2 * gamma * (z_center - z_plane); IPPL's
            // shiftedGreensFunction uses the opposite sign convention because it
            // evaluates G(r - shift). See TestShiftedGreensFunction.
            const auto origin = mesh.getOrigin();
            const int N_z =
                    static_cast<int>(this->getRho()->getLayout().getDomain()[Dim - 1].length());
            const double zCenter =
                    origin[Dim - 1] + 0.5 * static_cast<double>(N_z) * hrOrig[Dim - 1];
            ippl::Vector<double, Dim> shift(0.0);
            shift[Dim - 1] = 2.0 * gammaBin * (zCenter - shiftedGreensPlaneZ_m);

            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " shifted-GF runSolver start, plane=" << shiftedGreensPlaneZ_m
              << ", shift_z=" << shift[Dim - 1] << endl;
            this->runShiftedOpenSolver(shift, diagnostics);
            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " shifted-GF runSolver done; accumulate->Etmp (B negated, z-flip)" << endl;

            // Keep the real charge sign for the shifted solve. The image-charge
            // sign enters through the z-flip/component-sign rule in
            // FieldComposer. Pre-negating rho here would invert the
            // image field a second time and reinforce the near-cathode
            // transverse field instead of cancelling it.
            {
                auto compositionEvent = diagnostics.scopedEvent(
                        opalx::spacecharge::SelfFieldEventKind::FieldComposition,
                        "shifted-green-accumulation");
                typename FieldComposer_t::Policy compositionPolicy;
                for (unsigned d = 0; d < Dim; ++d) {
                    compositionPolicy.meanMomentum[d] = kinematics.pmean[d];
                }
                compositionPolicy.gamma        = gammaBin;
                compositionPolicy.magneticSign = -1.0;
                compositionPolicy.sourceRule =
                        opalx::spacecharge::FieldSourceRule::ShiftedGreenImageZ;
                fieldComposer_m.accumulate(workspace, compositionPolicy);
            }

            mesh.setMeshSpacing(hrOrig);
        }
    }

    // after all bins, gather the accumulated lab-frame field back to particles.
    {
        auto compositionEvent = diagnostics.scopedEvent(
                opalx::spacecharge::SelfFieldEventKind::FieldComposition, "final-gather");
        if (gatherAttribute_m != GatherAttribute::ElectricFieldE) {
            throw OpalException(
                    "BinnedFieldSolver::computeBinnedSelfFields",
                    "Unsupported gather attribute in binned solver.");
        }

        std::shared_ptr<ParticleCtr_t> pc = bunch.getParticleContainer();
        fieldComposer_m.gatherAccumulated(scatterGather_m, pc->E, pc->B, pc->R, workspace);
    }

    // per-call table: gammaBin / nParticles / binNumber.
    if (diagnostics.shouldPrintBinTable(bunch.getGlobalTrackStep())) {
        printBinStatsTable(bins->getBinningCmdName(), binStats);
    }
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::computeLegacySelfFields(
        PartBunch_t& bunch, opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
    diagnostics.recordBinCount(0);
    auto solveUnitEvent = diagnostics.scopedEvent(
            opalx::spacecharge::SelfFieldEventKind::SolveUnit, "whole-bunch");
    auto solvePassEvent = diagnostics.scopedEvent(
            opalx::spacecharge::SelfFieldEventKind::SolvePass, "primary-and-image");
    // This code is a direct move of the legacy implementation from
    // PartBunch::computeSelfFields (scatter/solve/gather for all particles).

    //  access the particle container for scattering/gathering.
    std::shared_ptr<ParticleCtr_t> pc = bunch.getParticleContainer();
    if (!pc) {
        throw OpalException(
                "BinnedFieldSolver::computeLegacySelfFields",
                "Bunch particle container is not available.");
    }

    // Level-5 debug: legacy mode entry.
    Inform m("BinnedFieldSolver::computeLegacySelfFields");
    m << level4 << "Legacy mode entry: localP=" << pc->getLocalNum()
      << ", totalP=" << pc->getTotalNum() << ", stype=" << this->getStype() << endl;

    if (scatterAttribute_m != ScatterAttribute::ChargeQ) {
        throw OpalException(
                "BinnedFieldSolver::computeLegacySelfFields",
                "Unsupported scatter attribute in legacy solver.");
    }

    Workspace_t& workspace          = this->getWorkspace();
    const auto& backendCapabilities = this->getBackendCapabilities();
    using Selection                 = typename ScatterGather_t::Selection;
    const auto selection            = Selection::direct(
            0, static_cast<typename ScatterGather_t::size_type>(pc->getLocalNum()));
    const auto depositKind = imageScatterController_m.isEnabled()
                                     ? ScatterGather_t::DepositKind::PrimaryAndImage
                                     : ScatterGather_t::DepositKind::Primary;

    typename ScatterGather_t::ChargeNormalization normalization;
    normalization.timeStep              = bunch.getdT();
    normalization.gamma                 = 1.0;
    normalization.selectedCharge        = pc->getTotalCharge();
    normalization.couplingConstant      = this->getCouplingConstant();
    normalization.normalizeByCellVolume = backendCapabilities.normalizeChargeByCellVolume;
    normalization.subtractNeutralizingBackground =
            backendCapabilities.subtractNeutralizingBackground;

    typename ScatterGather_t::ImagePolicy imagePolicy;
    imagePolicy.enabled = imageScatterController_m.isEnabled();
    imagePolicy.planeZ  = imageScatterController_m.getZPlane();

    scatterGather_m.depositCharge(
            *pc, workspace, depositKind, selection, normalization, imagePolicy);

    // Ensure deterministic output even for solver types that do not update `E`.
    *(this->getE()) = 0.0;

    // run the solver once and gather mesh E back to particles.
    m << level4 << "Legacy mode: runSolver() start" << endl;
    this->runSolver(false, diagnostics);
    dumpDirichletPlaneDiagnosticsIfRequested(bunch, "legacy", diagnostics);
    m << level4 << "Legacy mode: gather E->particles" << endl;

    // Gather solver output directly (legacy path does not use Etmp).
    if (gatherAttribute_m == GatherAttribute::ElectricFieldE) {
        auto compositionEvent = diagnostics.scopedEvent(
                opalx::spacecharge::SelfFieldEventKind::FieldComposition, "legacy-gather");
        fieldComposer_m.gatherElectrostatic(scatterGather_m, pc->E, pc->R, workspace);
    } else {
        throw OpalException(
                "BinnedFieldSolver::computeLegacySelfFields",
                "Unsupported gather attribute in legacy solver.");
    }

    // TABLEPRINTFREQ is binned-mode only; legacy mode intentionally prints nothing.
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::rebinAndPrepare(
        PartBunch_t& bunch, std::shared_ptr<AdaptBins_t> bins) {
    Inform m("BinnedFieldSolver::rebinAndPrepare");
    m << level4 << "Rebin start: maxBins=" << static_cast<int>(bins->getMaxBinCount())
      << ", adaptiveBinning=" << (adaptiveBinning_m ? 1 : 0) << endl;
    bins->doFullRebin(bins->getMaxBinCount());
    bunch.dumpBinConfig(true);
    if (adaptiveBinning_m) {
        bins->genAdaptiveHistogram();
        bunch.dumpBinConfig(false);
    }
    bins->sortContainerByBin();
    m << level4 << "Rebin done: currentBins=" << static_cast<int>(bins->getCurrentBinCount())
      << endl;
}

template <typename T, unsigned Dim>
typename BinnedFieldSolver<T, Dim>::BinKinematics BinnedFieldSolver<T, Dim>::computeGammaBinGlobal(
        PartBunch_t& bunch, std::shared_ptr<AdaptBins_t> bins, const bin_index_type binIndex,
        const size_type nPartGlobal) const {
    // compute global mean momentum and gamma for the merged bin.
    Inform m("BinnedFieldSolver::computeGammaBinGlobal");
    m << level4 << "gammaBinGlobal: binIndex=" << static_cast<int>(binIndex)
      << ", nPartGlobal=" << static_cast<unsigned long long>(nPartGlobal) << endl;

    typename particle_position_type::view_type pView = bunch.getParticleContainer()->P.getView();
    typename AdaptBins_t::hash_type indices          = bins->getHashArray();

    // compute local momentum and gamma sums over particles in this bin.
    Vector_t<double, Dim> localPsum(0.0);
    Kokkos::parallel_reduce(
            "BinnedFieldSolver::pmeanPerBin", bins->getBinIterationPolicy(binIndex),
            KOKKOS_LAMBDA(const size_type i, Vector_t<double, Dim>& sum) {
                sum += pView(indices(i));
            },
            localPsum);

    double localGammaSum = 0.0;
    Kokkos::parallel_reduce(
            "BinnedFieldSolver::gammaMeanPerBin", bins->getBinIterationPolicy(binIndex),
            KOKKOS_LAMBDA(const size_type i, double& sum) {
                const Vector_t<double, Dim> p = pView(indices(i));
                sum += Kokkos::sqrt(1.0 + p.dot(p));
            },
            localGammaSum);

    // reduce momentum sums across MPI ranks and normalize by `nPartGlobal`.
    Vector_t<double, Dim> globalPsum(0.0);
    ippl::Comm->allreduce(localPsum, 1, std::plus<Vector_t<double, Dim>>());
    globalPsum = localPsum;

    ippl::Comm->allreduce(localGammaSum, 1, std::plus<double>());

    BinKinematics kinematics;
    if (nPartGlobal == 0) {
        return kinematics;
    }

    kinematics.pmean    = globalPsum / static_cast<double>(nPartGlobal);
    kinematics.gammaBin = localGammaSum / static_cast<double>(nPartGlobal);
    return kinematics;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::prepareRhoForBin(
        PartBunch_t& bunch, std::shared_ptr<AdaptBins_t> bins, const bin_index_type binIndex,
        const size_type nPartGlobal, const double gammaBin, ImageScatterMode mode) {
    // Scatter bin charge to rho using the dt*Q deposition workflow.
    Inform m("BinnedFieldSolver::prepareRhoForBin");
    m << level4 << "prepareRho: binIndex=" << static_cast<int>(binIndex)
      << ", nPartGlobal=" << static_cast<unsigned long long>(nPartGlobal)
      << ", gammaBin=" << std::setprecision(10) << gammaBin << endl;

    // access particle views and validate scatter support.
    std::shared_ptr<ParticleCtr_t> pc = bunch.getParticleContainer();

    if (scatterAttribute_m != ScatterAttribute::ChargeQ) {
        throw OpalException(
                "BinnedFieldSolver::prepareRhoForBin",
                "Unsupported scatter attribute in binned solver.");
    }

    // Scatter bin charge to rho (with bin iteration policy and hash indexing).
    // Master approach: scale dt by Q, scatter dt, then restore dt.
    const auto policy       = bins->getBinIterationPolicy(binIndex);
    const auto hash         = bins->getHashArray();
    const size_type pBegin  = static_cast<size_type>(policy.begin());
    const size_type pEnd    = static_cast<size_type>(policy.end());
    const size_type hExtent = static_cast<size_type>(hash.extent(0));
    const size_type nLocal  = pc->getLocalNum();
    const char* modeName =
            mode == ImageScatterMode::PrimaryOnly
                    ? "PrimaryOnly"
                    : (mode == ImageScatterMode::ImageOnly ? "ImageOnly" : "PrimaryAndImage");

    if (pEnd > hExtent) {
        throw OpalException(
                "BinnedFieldSolver::prepareRhoForBin",
                "Bin scatter policy exceeds hash extent: policyEnd=" + std::to_string(pEnd)
                        + ", hashExtent=" + std::to_string(hExtent) + ".");
    }

    // If the selected bin spans every local particle, the hashed subset scatter is
    // equivalent to the all-local scatter. Prefer the all-local path here: it avoids
    // dereferencing the bin hash view in the hot scatter kernel and is the common
    // AWAGun early-emission case after 128 bins merge down to one bin.
    const bool scatterAllLocal = (pBegin == 0 && pEnd == nLocal);
    m << level5 << "prepareRho: scatter mode=" << modeName
      << ", path=" << (scatterAllLocal ? "all-local" : "subset")
      << ", localP=" << static_cast<unsigned long long>(nLocal) << ", policy=[" << pBegin << ","
      << pEnd << "), hashExtent=" << static_cast<unsigned long long>(hExtent) << endl;

    using Selection = typename ScatterGather_t::Selection;
    const Selection selection =
            scatterAllLocal ? Selection::direct(0, nLocal) : Selection::indexed(policy, hash);

    typename ScatterGather_t::DepositKind depositKind =
            ScatterGather_t::DepositKind::PrimaryAndImage;
    if (mode == ImageScatterMode::PrimaryOnly) {
        depositKind = ScatterGather_t::DepositKind::Primary;
    } else if (mode == ImageScatterMode::ImageOnly) {
        depositKind = ScatterGather_t::DepositKind::Image;
    }

    const auto& backendCapabilities = this->getBackendCapabilities();

    typename ScatterGather_t::ChargeNormalization normalization;
    normalization.timeStep         = bunch.getdT();
    normalization.gamma            = gammaBin;
    normalization.selectedCharge   = pc->getChargePerParticle() * static_cast<double>(nPartGlobal);
    normalization.couplingConstant = this->getCouplingConstant();
    normalization.normalizeByCellVolume = backendCapabilities.normalizeChargeByCellVolume;
    normalization.subtractNeutralizingBackground =
            backendCapabilities.subtractNeutralizingBackground;

    typename ScatterGather_t::ImagePolicy imagePolicy;
    imagePolicy.enabled = imageScatterController_m.isEnabled();
    imagePolicy.planeZ  = imageScatterController_m.getZPlane();

    scatterGather_m.depositCharge(
            *pc, this->getWorkspace(), depositKind, selection, normalization, imagePolicy);
}
