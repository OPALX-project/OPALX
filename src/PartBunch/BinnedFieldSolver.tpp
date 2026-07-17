#include "Structure/DataSink.h"

#include <cstring>
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
            // Legacy path has no separate correction pass: PicScatterGather deposits
            // primary and image charge together before one standard solve. The shifted Green's
            // correction is only
            // implemented for the binned path, so warn once if the user requested it
            // without binning.
            if (shiftedGreensWasEnabled && shiftedGreensActiveThisStep) {
                m << level3 << "SHIFTED_GREENS_FUNCTION is set but no binning is active; "
                  << "the legacy path does not apply the Dirichlet correction." << endl;
            }
        }

        m << level4 << "Dispatching prepared iteration plan." << endl;
        executeIterationPlan(bunch, iterationPlan, prepared, particleGeneration, diagnostics);
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
void BinnedFieldSolver<T, Dim>::executeIterationPlan(
        PartBunch_t& bunch, IterationPlan_t& iterationPlan, const PreparedIteration_t& prepared,
        std::uint64_t particleGeneration, opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
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
            computeLegacySelfFields(bunch, unit, diagnostics);
            continue;
        }
        if (!binned
            || unit.fieldMode != opalx::spacecharge::SolveUnitFieldMode::LorentzTransformed) {
            throw OpalException(
                    "BinnedFieldSolver::executeIterationPlan",
                    "A Lorentz-transformed solve unit requires a binning plan.");
        }
        if (unit.ordinal >= nBins) {
            throw OpalException(
                    "BinnedFieldSolver::executeIterationPlan",
                    "A binning plan emitted a solve-unit ordinal outside its prepared range.");
        }

        const std::size_t binIndex    = unit.ordinal;
        const std::size_t nPartGlobal = unit.globalParticleCount;

        m << level4 << "binIndex=" << static_cast<int>(binIndex)
          << " nPartGlobal=" << static_cast<unsigned long long>(nPartGlobal) << endl;

        const double gammaBin = unit.gamma;
        if (gammaBin <= 0.0) {
            throw OpalException(
                    "BinnedFieldSolver::executeIterationPlan",
                    "Computed non-positive gamma for bin.");
        }

        m << level4 << "binIndex=" << static_cast<int>(binIndex)
          << " gammaBin=" << std::setprecision(10) << gammaBin << endl;

        binStats_m.push_back(
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
            prepareRhoForUnit(bunch, unit, scatterMode);

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
                    compositionPolicy.meanMomentum[d] = unit.meanMomentum[d];
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
            prepareRhoForUnit(bunch, unit, ImageScatterMode::ImageOnly);

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
                    compositionPolicy.meanMomentum[d] = unit.meanMomentum[d];
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
            prepareRhoForUnit(bunch, unit, ImageScatterMode::PrimaryOnly);

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
                    compositionPolicy.meanMomentum[d] = unit.meanMomentum[d];
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
void BinnedFieldSolver<T, Dim>::computeLegacySelfFields(
        PartBunch_t& bunch, const SolveUnit_t& unit,
        opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
    if (unit.kind != opalx::spacecharge::IterationKind::WholeBunch
        || unit.fieldMode != opalx::spacecharge::SolveUnitFieldMode::Direct) {
        throw OpalException(
                "BinnedFieldSolver::computeLegacySelfFields",
                "A whole-bunch plan emitted an incompatible solve unit.");
    }

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
    const auto selection            = unit.depositSelection();
    const auto depositKind          = imageScatterController_m.isEnabled()
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
void BinnedFieldSolver<T, Dim>::prepareRhoForUnit(
        PartBunch_t& bunch, const SolveUnit_t& unit, ImageScatterMode mode) {
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
    const char* modeName =
            mode == ImageScatterMode::PrimaryOnly
                    ? "PrimaryOnly"
                    : (mode == ImageScatterMode::ImageOnly ? "ImageOnly" : "PrimaryAndImage");

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
