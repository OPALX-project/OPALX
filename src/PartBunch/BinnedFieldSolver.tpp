#include "Structure/DataSink.h"

#include <cstring>
#include <utility>
#include <vector>

template <typename T, unsigned Dim>
BinnedFieldSolver<T, Dim>::BinnedFieldSolver(
        std::string solver,
        Field_t<Dim>* rho,
        VField_t<T, Dim>* E,
        Field_t<Dim>* phi,
        std::shared_ptr<BCHandler_t> bcHandler,
        int tablePrintFrequency)
    : FieldSolver<T, Dim>(solver, rho, E, phi, bcHandler) {
    scatterAttribute_m    = ScatterAttribute::ChargeQ;
    gatherAttribute_m     = GatherAttribute::ElectricFieldE;
    tablePrintFrequency_m = tablePrintFrequency;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::computeSelfFields(std::shared_ptr<PartBunch_t> bunch) {
    // validate inputs and decide between binned vs legacy solver.
    if (!bunch) {
        throw OpalException("BinnedFieldSolver::computeSelfFields", "Passed nullptr bunch.");
    }

    // access the particle container for trivial-case checks.
    std::shared_ptr<ParticleCtr_t> pc = bunch->getParticleContainer();
    if (!pc) {
        throw OpalException(
                "BinnedFieldSolver::computeSelfFields",
                "Bunch particle container is not available.");
    }

    Inform m("BinnedFieldSolver::computeSelfFields");
    // TYPE=NONE is a true no-op: skip all binning/scatter/solve/gather work.
    if (this->getStype() == "NONE") {
        // Already called in ParallelTracker::resetFields()
        // pc->E = 0.0;
        // pc->B = 0.0;
        m << level5 << "Skipping scatter/gather and self-field computation for NONE solver."
          << endl;
        return;
    }

    // trivial case where self-field has no effect.
    if (ippl::Comm->size() == 1 && pc->getLocalNum() <= 1) {
        pc->E = 0.0;
        return;
    }

    // decide which solver path to run (binned vs legacy).
    const bool hasBins = bunch->hasBinning();

    m << level4 << "Entry: rank=" << ippl::Comm->rank() << ", localParticles=" << pc->getLocalNum()
      << ", totalParticles=" << pc->getTotalNum() << ", hasBins=" << (hasBins ? 1 : 0)
      << ", stype=" << this->getStype() << endl;

    // Temporarily disable image charges if the step limit has been reached.
    // The controller's enabled flag gates the image scatter pass; by disabling it
    // here both the legacy and binned paths automatically perform primary-only scatter.
    const bool imageWasEnabled = imageScatterController_m.isEnabled();
    const bool imageActiveThisStep = isImageChargeActiveForStep(bunch->getGlobalTrackStep());
    if (imageWasEnabled && !imageActiveThisStep) {
        m << level3 << "ZEROFACE_MAXSTEPS reached (step=" << bunch->getGlobalTrackStep()
          << ", maxSteps=" << zerofaceMaxSteps_m
          << "); disabling image charges for this step." << endl;
        imageScatterController_m.configure(false, imageScatterController_m.getZPlane());
    }

    // Mirror the same step-budget toggling for the shifted Green's path.
    const bool shiftedGreensWasEnabled   = shiftedGreensEnabled_m;
    const bool shiftedGreensActiveThisStep =
            isShiftedGreensActiveForStep(bunch->getGlobalTrackStep());
    if (shiftedGreensWasEnabled && !shiftedGreensActiveThisStep) {
        m << level3 << "ZEROFACE_MAXSTEPS reached (step=" << bunch->getGlobalTrackStep()
          << ", maxSteps=" << zerofaceMaxSteps_m
          << "); disabling SHIFTED_GREENS_FUNCTION correction for this step." << endl;
        shiftedGreensEnabled_m = false;
    }

    if (hasBins) {
        m << level4 << "Dispatching to computeBinnedSelfFields() (binned path)." << endl;
        computeBinnedSelfFields(bunch);
    } else {
        // Legacy path has no separate correction pass: it scatters primary+image
        // in one shot via ImageChargeScatterController::scatterPrimaryAndImage
        // and does one standard solve. The shifted Green's correction is only
        // implemented for the binned path, so warn once if the user requested it
        // without binning.
        if (shiftedGreensWasEnabled && shiftedGreensActiveThisStep) {
            m << level3 << "SHIFTED_GREENS_FUNCTION is set but no binning is active; "
              << "the legacy path does not apply the Dirichlet correction." << endl;
        }
        m << level4 << "Dispatching to computeLegacySelfFields() (legacy path)." << endl;
        computeLegacySelfFields(bunch);
    }

    // Restore image-charge controller state if it was temporarily disabled.
    if (imageWasEnabled && !imageActiveThisStep) {
        imageScatterController_m.configure(true, imageScatterController_m.getZPlane());
    }
    if (shiftedGreensWasEnabled && !shiftedGreensActiveThisStep) {
        shiftedGreensEnabled_m = true;
    }
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
                "BinnedFieldSolver::setZerofaceMaxSteps",
                "ZEROFACE_MAXSTEPS must be >= 0.");
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
    zeroFacePlaneDumpFrequency_m = frequency;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::dumpDirichletPlaneDiagnosticsIfRequested(
        std::shared_ptr<PartBunch_t> bunch, const std::string& solveTag) {
    if (!imageScatterController_m.isEnabled() || zeroFacePlaneDumpFrequency_m <= 0) {
        return;
    }

    const long long step = bunch->getGlobalTrackStep();
    if (step < 0 || (step % zeroFacePlaneDumpFrequency_m) != 0) {
        return;
    }

    Inform m("BinnedFieldSolver::dumpDirichletPlaneDiagnosticsIfRequested");

    if (ippl::Comm->size() != 1) {
        if (!warnedPlaneDumpParallelUnsupported_m) {
            warnedPlaneDumpParallelUnsupported_m = true;
            m << level3
              << "Dirichlet-plane diagnostics currently support only single-rank runs. "
              << "Skipping dump and statistics output." << endl;
        }
        return;
    }

    Field_t<Dim>* potentialField =
            (this->getStype() == "CG") ? this->getPhi() : this->getRho();
    if (!potentialField) {
        return;
    }
    const double zPlane      = imageScatterController_m.getZPlane();

    std::shared_ptr<DataSink> dataSink = bunch->getDataSink();
    if (!dataSink) {
        return;
    }

    const auto diagnostics =
            dataSink->dumpDirichletPlane(step, bunch->getT(), zPlane, *potentialField, solveTag);
    if (diagnostics.sampleCount == 0) {
        return;
    }

    m << level2 << "Dirichlet-plane potential diagnostics (" << solveTag << ") at step " << step
      << ": z=" << zPlane << " m, mean(phi)=" << diagnostics.mean
      << " V, var(phi)=" << diagnostics.variance << " V^2" << endl;
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
void BinnedFieldSolver<T, Dim>::computeBinnedSelfFields(std::shared_ptr<PartBunch_t> bunch) {
    // execute full binned self-field algorithm.
    // fetch the adaptive bin structure.
    std::shared_ptr<AdaptBins_t> bins = bunch->getBins();
    if (!bins) {
        // Defensive: runtime selection above should prevent this.
        computeLegacySelfFields(bunch);
        return;
    }

    // build and merge adaptive bins for this step.
    // build and merge adaptive bins for this step.
    rebinAndPrepare(bunch, bins);

    // obtain the temporary E buffer used to accumulate bin contributions.
    std::shared_ptr<VField_t<T, Dim>> EtmpSP = bunch->getTempEField();
    if (!EtmpSP) {
        throw OpalException(
                "BinnedFieldSolver::computeBinnedSelfFields",
                "Temporary E field (Etmp) is not initialized.");
    }
    std::shared_ptr<VField_t<T, Dim>> BtmpSP = bunch->getTempBField();
    if (!BtmpSP) {
        throw OpalException(
                "BinnedFieldSolver::computeBinnedSelfFields",
                "Temporary B field (Btmp) is not initialized.");
    }

    VField_t<T, Dim>& Etmp = *EtmpSP;
    VField_t<T, Dim>& Btmp = *BtmpSP;
    // clear the accumulation buffer.
    Etmp = 0.0;
    Btmp = 0.0;

    // determine the number of bins used for this step.
    const bin_index_type nBins = bins->getCurrentBinCount();

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
        auto& mesh = this->getRho()->get_mesh();
        const auto hrOrig = mesh.getMeshSpacing();
        auto hrStretched = hrOrig;
        hrStretched[Dim - 1] *= gammaBin;

        const bool imageActive         = imageScatterController_m.isEnabled();
        const bool shiftedGreensActive = shiftedGreensEnabled_m;
        // Mutual exclusion enforced at config time; defensive assert.
        assert(!(imageActive && shiftedGreensActive));
        const bool correctionActive = imageActive || shiftedGreensActive;

        // --- Primary pass: scatter real charges, solve, accumulate with +B ---
        {
            const ImageScatterMode scatterMode = correctionActive
                    ? ImageScatterMode::PrimaryOnly
                    : ImageScatterMode::PrimaryAndImage;
            prepareRhoForBin(bunch, bins, binIndex, nPartGlobal, gammaBin, scatterMode);

            *(this->getE()) = 0.0;
            mesh.setMeshSpacing(hrStretched);

            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " primary runSolver(true) start"
              << " (hr_z stretched by gamma=" << gammaBin << ")" << endl;
            this->runSolver(true);
            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " primary runSolver(true) done; accumulate->Etmp" << endl;

            // DEBUG: report max|E| after the primary solve AND Etmp before/after.
            auto e_max_primary = [&](const char* tag, auto& field) {
                double localMax = 0.0;
                auto vE = field.getView();
                ippl::parallel_reduce(
                    tag, field.getFieldRangePolicy(),
                    KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx,
                                  double& lmax) {
                        auto e = apply(vE, idx);
                        const double m2 = e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
                        if (m2 > lmax) lmax = m2;
                    },
                    Kokkos::Max<double>(localMax));
                double globalMax = 0.0;
                ippl::Comm->allreduce(localMax, globalMax, 1, std::greater<double>());
                return std::sqrt(globalMax);
            };
            m << level5 << "[DEBUG] binIndex=" << static_cast<int>(binIndex)
              << " primary: max|E| = " << e_max_primary("primary |E|", *this->getE())
              << " max|Etmp before| = " << e_max_primary("Etmp pre-primary", *EtmpSP) << endl;

            accumulateFieldToTemp(gammaBin, kinematics.pmean, EtmpSP, BtmpSP, +1.0);

            m << level5 << "[DEBUG] binIndex=" << static_cast<int>(binIndex)
              << " primary: max|Etmp after| = "
              << e_max_primary("Etmp post-primary", *EtmpSP) << endl;

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
        // + sign-flip on the solver output, baked into accumulateFieldToTemp,
        // produces the image-charge E field directly. Works at any bunch-to-
        // plane distance and requires the OPEN solver (checked in
        // FieldSolver::runShiftedOpenSolver).
        if (imageActive) {
            prepareRhoForBin(bunch, bins, binIndex, nPartGlobal, gammaBin,
                             ImageScatterMode::ImageOnly);

            *(this->getE()) = 0.0;
            mesh.setMeshSpacing(hrStretched);

            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " image runSolver(true) start" << endl;
            this->runSolver(true);
            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " image runSolver(true) done; accumulate->Etmp (B negated)" << endl;

            accumulateFieldToTemp(gammaBin, kinematics.pmean, EtmpSP, BtmpSP, -1.0);
            mesh.setMeshSpacing(hrOrig);

            // Dump phi ~= 0 check on the Dirichlet plane AFTER the correction
            // lands on the mesh. Only safe for the explicit-image path: there
            // the cathode plane always lies inside the computational domain.
            // For the shifted path the domain may be far from z = R0Z, so the
            // interpolated phi on the plane would be meaningless.
            if (!dumpedDirichletPlaneThisStep) {
                dumpDirichletPlaneDiagnosticsIfRequested(bunch, "binned");
                dumpedDirichletPlaneThisStep = true;
            }
        }
        else if (shiftedGreensActive) {
            // Shifted-Green's-function image correction. Multi-rank enabled: the
            // axis-flip source read crosses ranks under PARFFTZ, but that is handled
            // inside accumulateFieldToTemp (it calls buildFlippedZSlab to stage a
            // local view populated via peer-rank MPI exchange).
            //
            // Re-scatter the primary charges (solve() overwrote the RHS in the
            // primary pass). This matches the ImageOnly path's pattern.
            prepareRhoForBin(bunch, bins, binIndex, nPartGlobal, gammaBin,
                             ImageScatterMode::PrimaryOnly);

            *(this->getE()) = 0.0;
            mesh.setMeshSpacing(hrStretched);

            // Shift formula in stretched (rest-frame) coordinates:
            //   shift_z = L + 2*origin_z - 2*R0Z = 2 * (z_center_rest - R0Z).
            // Origin is in lab-frame z; hrStretched[Dim-1] is the rest-frame
            // z-spacing. See the TestShiftedGreensFunction derivation.
            const auto origin = mesh.getOrigin();
            const int    N_z  = static_cast<int>(
                    this->getRho()->getLayout().getDomain()[Dim - 1].length());
            const double z_center_rest =
                    origin[Dim - 1] + 0.5 * static_cast<double>(N_z) * hrStretched[Dim - 1];
            ippl::Vector<double, Dim> shift(0.0);
            shift[Dim - 1] = 2.0 * (z_center_rest - shiftedGreensPlaneZ_m);

            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " shifted-GF runSolver start, plane=" << shiftedGreensPlaneZ_m
              << ", shift_z=" << shift[Dim - 1] << endl;
            this->runShiftedOpenSolver(shift);
            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " shifted-GF runSolver done; accumulate->Etmp (B negated, z-flip)"
              << endl;

            // Axis-flip + component-wise sign is handled inside
            // accumulateFieldToTemp when flipAxis >= 0 (see method doc).
            constexpr int zFlipAxis = static_cast<int>(Dim) - 1;
            accumulateFieldToTemp(gammaBin, kinematics.pmean, EtmpSP, BtmpSP, -1.0,
                                  zFlipAxis);

            mesh.setMeshSpacing(hrOrig);
        }
    }

    // after all bins, gather the accumulated lab-frame field back to particles.
    gatherFromTempToParticles(bunch, EtmpSP, BtmpSP);

    // per-call table: gammaBin / nParticles / binNumber.
    if (tablePrintFrequency_m > 0) {
        const long long step = bunch->getGlobalTrackStep();
        if (step >= 0 && (step % tablePrintFrequency_m) == 0) {
            printBinStatsTable(bins->getBinningCmdName(), binStats);
        }
    }
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::computeLegacySelfFields(std::shared_ptr<PartBunch_t> bunch) {
    // This code is a direct move of the legacy implementation from
    // PartBunch::computeSelfFields (scatter/solve/gather for all particles).

    //  access the particle container for scattering/gathering.
    std::shared_ptr<ParticleCtr_t> pc = bunch->getParticleContainer();
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

    typename PartBunch_t::Base::particle_position_type* R = &pc->R;

    Field_t<Dim>* rho = this->getRho();
    *rho              = 0.0;

    // Scatter charge to mesh rho using dt-weighted deposition (master approach):
    // scale dt by Q, scatter dt, then restore dt.
    imageScatterController_m.scatterPrimaryAndImage(pc, *R, *rho);

    // Rho normalization for fractional time steps.
    (*rho) = (*rho) / bunch->getdT();

    //  apply mesh normalization, background subtraction, and rho scaling.
    const std::string stype = this->getStype();
    if (stype != "FEM" && stype != "FEM_PRECON") {
        const double cellVolume =
                std::reduce(bunch->hr_m.begin(), bunch->hr_m.end(), 1.0, std::multiplies<double>());
        (*rho) = (*rho) / cellVolume;
    }

    // Alpine uses net-0 charge for non-OPEN solvers (periodic BCs).
    if (stype != "OPEN") {
        double size = 1.0;
        for (size_t d = 0; d < Dim; ++d) {
            size *= bunch->rmax_m[d] - bunch->rmin_m[d];
        }

        const double totalQ = bunch->getParticleContainer()->getTotalCharge();
        (*rho)              = (*rho) - (totalQ / size);
    }

    (*rho) = (*rho) * this->getCouplingConstant();

    // Ensure deterministic output even for solver types that do not update `E`.
    *(this->getE()) = 0.0;

    // run the solver once and gather mesh E back to particles.
    m << level4 << "Legacy mode: runSolver() start" << endl;
    this->runSolver();
    dumpDirichletPlaneDiagnosticsIfRequested(bunch, "legacy");
    m << level4 << "Legacy mode: gather E->particles" << endl;

    // Gather solver output directly (legacy path does not use Etmp).
    if (gatherAttribute_m == GatherAttribute::ElectricFieldE) {
        gather(pc->E, *this->getE(), *R);
    } else {
        throw OpalException(
                "BinnedFieldSolver::computeLegacySelfFields",
                "Unsupported gather attribute in legacy solver.");
    }

    // TABLEPRINTFREQ is binned-mode only; legacy mode intentionally prints nothing.
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::rebinAndPrepare(
        std::shared_ptr<PartBunch_t> bunch, std::shared_ptr<AdaptBins_t> bins) {
    // adaptive histogram configuration.
    // execute full rebin and generate the adaptive histogram (bin merging).
    Inform m("BinnedFieldSolver::rebinAndPrepare");
    m << level4 << "Rebin start: maxBins=" << static_cast<int>(bins->getMaxBinCount()) << endl;
    bins->doFullRebin(bins->getMaxBinCount());
    bunch->dumpBinConfig(true);
    bins->sortContainerByBin();
    bins->genAdaptiveHistogram();
    bunch->dumpBinConfig(false);
    m << level4 << "Rebin done: currentBins=" << static_cast<int>(bins->getCurrentBinCount())
      << endl;
}

template <typename T, unsigned Dim>
typename BinnedFieldSolver<T, Dim>::BinKinematics BinnedFieldSolver<T, Dim>::computeGammaBinGlobal(
        std::shared_ptr<PartBunch_t> bunch, std::shared_ptr<AdaptBins_t> bins,
        const bin_index_type binIndex, const size_type nPartGlobal) const {
    // compute global mean momentum and gamma for the merged bin.
    Inform m("BinnedFieldSolver::computeGammaBinGlobal");
    m << level4 << "gammaBinGlobal: binIndex=" << static_cast<int>(binIndex)
      << ", nPartGlobal=" << static_cast<unsigned long long>(nPartGlobal) << endl;

    typename particle_position_type::view_type pView = bunch->getParticleContainer()->P.getView();
    typename AdaptBins_t::hash_type indices          = bins->getHashArray();

    // compute local momentum sums over particles in this bin.
    Vector_t<double, Dim> localPsum(0.0);
    Kokkos::parallel_reduce(
            "BinnedFieldSolver::pmeanPerBin", bins->getBinIterationPolicy(binIndex),
            KOKKOS_LAMBDA(const size_type i, Vector_t<double, Dim>& sum) {
                sum += pView(indices(i));
            },
            localPsum);

    // reduce momentum sums across MPI ranks and normalize by `nPartGlobal`.
    Vector_t<double, Dim> globalPsum(0.0);
    ippl::Comm->allreduce(localPsum, 1, std::plus<Vector_t<double, Dim>>());
    globalPsum = localPsum;

    BinKinematics kinematics;
    if (nPartGlobal == 0) {
        return kinematics;
    }

    kinematics.pmean    = globalPsum / static_cast<double>(nPartGlobal);
    kinematics.gammaBin = Kokkos::sqrt(1.0 + kinematics.pmean.dot(kinematics.pmean));
    return kinematics;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::prepareRhoForBin(
        std::shared_ptr<PartBunch_t> bunch, std::shared_ptr<AdaptBins_t> bins,
        const bin_index_type binIndex, const size_type nPartGlobal, const double gammaBin,
        ImageScatterMode mode) {
    // Scatter bin charge to rho using dt-weighted deposition.
    // If the ParticleContainer supports scaleDtByCharge(), use the master approach:
    // scale dt by charge, scatter dt, then unscale.
    // Otherwise, fall back to temporarily scaling/scattering Q by dt.
    Inform m("BinnedFieldSolver::prepareRhoForBin");
    m << level4 << "prepareRho: binIndex=" << static_cast<int>(binIndex)
      << ", nPartGlobal=" << static_cast<unsigned long long>(nPartGlobal)
      << ", gammaBin=" << std::setprecision(10) << gammaBin << endl;

    Field_t<Dim>* rho = this->getRho();
    *rho              = 0.0;

    // access particle views and validate scatter support.
    std::shared_ptr<ParticleCtr_t> pc                     = bunch->getParticleContainer();
    typename PartBunch_t::Base::particle_position_type* R = &pc->R;

    if (scatterAttribute_m != ScatterAttribute::ChargeQ) {
        throw OpalException(
                "BinnedFieldSolver::prepareRhoForBin",
                "Unsupported scatter attribute in binned solver.");
    }

    // Scatter bin charge to rho (with bin iteration policy and hash indexing).
    // Master approach: scale dt by Q, scatter dt, then restore dt.
    if (mode == ImageScatterMode::PrimaryOnly) {
        imageScatterController_m.scatterPrimaryOnly(
                pc, *R, *rho, bins->getBinIterationPolicy(binIndex), bins->getHashArray());
    } else if (mode == ImageScatterMode::ImageOnly) {
        imageScatterController_m.scatterImageOnly(
                pc, *R, *rho, bins->getBinIterationPolicy(binIndex), bins->getHashArray());
    } else {
        imageScatterController_m.scatterPrimaryAndImage(
                pc, *R, *rho, bins->getBinIterationPolicy(binIndex), bins->getHashArray());
    }

    // normalize rho for fractional time steps and mesh conventions.
    (*rho) = (*rho) / bunch->getdT();

    const std::string stype = this->getStype();
    if (stype != "FEM" && stype != "FEM_PRECON") {
        const double cellVolume =
                std::reduce(bunch->hr_m.begin(), bunch->hr_m.end(), 1.0, std::multiplies<double>());
        (*rho) = (*rho) / cellVolume;
    }

    // subtract non-OPEN background and apply Lorentz rest-frame scaling.
    // Background subtraction for non-OPEN solvers. Here we subtract only the bin's charge.
    if (stype != "OPEN") {
        double size = 1.0;
        for (size_t d = 0; d < Dim; ++d) {
            size *= bunch->rmax_m[d] - bunch->rmin_m[d];
        }

        const double totalQBin =
            bunch->getParticleContainer()->getChargePerParticle() * static_cast<double>(nPartGlobal);
        (*rho)                 = (*rho) - (totalQBin / size);
    }

    // Lorentz transform of charge density to the bin rest frame (thesis Eq. step 7).
    (*rho) = (*rho) / gammaBin;
    (*rho) = (*rho) * this->getCouplingConstant();
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::accumulateFieldToTemp(
        const double gammaBin, const Vector_t<double, Dim>& pmean,
        std::shared_ptr<VField_t<T, Dim>> EtmpSP, std::shared_ptr<VField_t<T, Dim>> BtmpSP,
        double bFieldSign, int flipAxis) {
    // transform rest-frame fields to lab-frame fields and accumulate.
    Inform m("BinnedFieldSolver::accumulateFieldToTemp");
    m << level4 << "accumulate: gammaBin=" << std::setprecision(10) << gammaBin
      << ", flipAxis=" << flipAxis << endl;

    const double invGamma         = (gammaBin > 0.0) ? (1.0 / gammaBin) : 0.0;
    const Vector_t<double, Dim> v = Physics::c * pmean * invGamma;
    const double vNorm            = Kokkos::sqrt(v.dot(v));
    const Vector_t<double, Dim> w = (vNorm > 0.0) ? (v / vNorm) : Vector_t<double, Dim>(0.0);
    const double gammaMinusOne    = gammaBin - 1.0;
    const double gammaOverCSq     = gammaBin / (Physics::c * Physics::c);

    const VField_t<T, Dim>& Eprime = *(this->getE());
    VField_t<T, Dim>& Etmp         = *EtmpSP;
    VField_t<T, Dim>& Btmp         = *BtmpSP;
    auto ePrimeView                = Eprime.getView();
    auto eTmpView                  = Etmp.getView();
    auto bTmpView                  = Btmp.getView();

    const int nghost = Eprime.getNghost();
    // Axis-flip geometry: for physical index k_phys in [0, N), the mirrored
    // physical index is N-1-k_phys. In view (ghost-padded) coordinates that is
    //   k_view_flipped = (N - 1 - k_phys) + nghost = 2*nghost + N - 1 - k_view.
    // view.extent(d) equals N + 2*nghost, so (view.extent(d) - 1 - idx) gives
    // the flipped view index directly.
    const int capturedFlipAxis = flipAxis;

    // parallel element-wise transformation and accumulation into temporaries.
    if (capturedFlipAxis < 0) {
        // Fast path: no flip, original behavior.
        ippl::parallel_for(
                "BinnedFieldSolver::accumulateFieldToTemp", Eprime.getFieldRangePolicy(),
                KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx) {
                    Vector_t<T, Dim> ePrime = apply(ePrimeView, idx);
                    const T ePrimeDotW      = ePrime.dot(w);
                    Vector_t<T, Dim> eLab   = gammaBin * ePrime - gammaMinusOne * ePrimeDotW * w;
                    Vector_t<T, Dim> bLab   = bFieldSign * gammaOverCSq * cross(v, ePrime);
                    Vector_t<T, Dim> eTotal = apply(eTmpView, idx);
                    Vector_t<T, Dim> bTotal = apply(bTmpView, idx);
                    eTotal += eLab;
                    bTotal += bLab;
                    apply(eTmpView, idx) = eTotal;
                    apply(bTmpView, idx) = bTotal;
                });
    } else {
        // Axis-flipped path: read E' at the GLOBALLY-flipped source index and apply
        // the component-wise image-charge sign rule (flip all components except the
        // one parallel to the flip axis) BEFORE the Lorentz transform.
        //
        // Under PARFFTZ=true the flipped z-index generally lives on a peer rank, so
        // we first populate flippedZSlab_m via a peer-rank MPI exchange. After the
        // call, flippedZSlab_m(i, j, k) == ePrimeView(i, j, flipped_k_viewindex)
        // for each local (i, j, k) — no cross-rank read required in the lambda.
        const int flipAxisCap = capturedFlipAxis;
        if (flipAxisCap != static_cast<int>(Dim) - 1) {
            throw OpalException(
                    "BinnedFieldSolver::accumulateFieldToTemp",
                    "flipAxis != Dim-1 not supported (only z-axis flip implemented).");
        }
        (void)nghost;

        this->buildFlippedZSlab(Eprime);
        auto flippedView = flippedZSlab_m;

        ippl::parallel_for(
                "BinnedFieldSolver::accumulateFieldToTemp[flipped]",
                Eprime.getFieldRangePolicy(),
                KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx) {
                    // Read pre-flipped E' at the local (i, j, k).
                    Vector_t<T, Dim> ePrime = flippedView(idx[0], idx[1], idx[2]);

                    // Component-wise sign rule (derivation: phi_image(r) = -phi_shifted(R(r))
                    // so E_image_i = -E_shifted_i(R(r)) for i != flipAxis and
                    // E_image_i = +E_shifted_i(R(r)) for i == flipAxis).
                    for (unsigned d = 0; d < Dim; ++d) {
                        if (static_cast<int>(d) != flipAxisCap) {
                            ePrime[d] = -ePrime[d];
                        }
                    }

                    const T ePrimeDotW      = ePrime.dot(w);
                    Vector_t<T, Dim> eLab   = gammaBin * ePrime - gammaMinusOne * ePrimeDotW * w;
                    Vector_t<T, Dim> bLab   = bFieldSign * gammaOverCSq * cross(v, ePrime);
                    Vector_t<T, Dim> eTotal = apply(eTmpView, idx);
                    Vector_t<T, Dim> bTotal = apply(bTmpView, idx);
                    eTotal += eLab;
                    bTotal += bLab;
                    apply(eTmpView, idx) = eTotal;
                    apply(bTmpView, idx) = bTotal;
                });
    }
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::buildFlippedZSlab(const VField_t<T, Dim>& src) {
    // Populate flippedZSlab_m with src z-globally-flipped, staged through host mirrors
    // for MPI. After this call:
    //   flippedZSlab_m(i, j, k_view) == src_view(i, j, flipped_k_view)
    // for each local (i, j) and each local k in the physical range [nghost, nghost+Nz_phys),
    // where flipped_k_view corresponds to the GLOBAL flip k_glob -> N_z_global - 1 - k_glob.
    // Ghost z-indices in flippedZSlab_m are zero-initialized but not used by the caller
    // (the accumulate lambda iterates src.getFieldRangePolicy() which excludes ghosts).

    const auto& layout = src.getLayout();
    const auto& ldom   = layout.getLocalNDIndex();
    const auto& gdom   = layout.getDomain();
    const int   nghost = src.getNghost();

    const int Nz_glob = static_cast<int>(gdom[Dim - 1].length());
    const int k0_r    = static_cast<int>(ldom[Dim - 1].first());
    const int Nz_phys = static_cast<int>(ldom[Dim - 1].length());
    const int k1_r    = k0_r + Nz_phys;

    auto srcView        = src.getView();
    const size_t Xf     = srcView.extent(0);
    const size_t Yf     = srcView.extent(1);
    const size_t Zf     = srcView.extent(2);
    const size_t xyCells = Xf * Yf;

    // Resize scratch if needed, zero-init (only ghost z-slices will remain zero).
    if (flippedZSlab_m.extent(0) != Xf || flippedZSlab_m.extent(1) != Yf
        || flippedZSlab_m.extent(2) != Zf) {
        Kokkos::realloc(flippedZSlab_m, Xf, Yf, Zf);
    }
    Kokkos::deep_copy(flippedZSlab_m, Vector_t<T, Dim>(0.0));

    // Host mirror of src (copy from device if needed).
    auto srcHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, srcView);
    // Host mirror of destination — we build on host then deep_copy back to device.
    auto dstHost = Kokkos::create_mirror_view(Kokkos::HostSpace{}, flippedZSlab_m);
    Kokkos::deep_copy(dstHost, Vector_t<T, Dim>(0.0));

    const int P = ippl::Comm->size();
    const int r = ippl::Comm->rank();

    const MPI_Comm comm = ippl::Comm->getCommunicator();
    const int tag       = 0x5F1F;  // "shifted-GF flip"

    const auto& hDomains = layout.getHostLocalDomains();

    // Pack helpers --------------------------------------------------------------------
    // Buffer format per-slab: (Xf * Yf * slab_len) cells of Vector_t<T, Dim>.
    const size_t bytesPerKSlab = sizeof(Vector_t<T, Dim>) * xyCells;

    std::vector<MPI_Request> requests;
    std::vector<std::vector<char>> sendBufs;
    std::vector<std::vector<char>> recvBufs;
    // For each recv buffer, remember (recvGlobStart, recvLen) so we can unpack into
    // the correct dst k range after MPI_Waitall.
    std::vector<std::pair<int, int>> recvMeta;

    requests.reserve(static_cast<size_t>(2 * P));
    sendBufs.reserve(static_cast<size_t>(P));
    recvBufs.reserve(static_cast<size_t>(P));
    recvMeta.reserve(static_cast<size_t>(P));

    for (int p = 0; p < P; ++p) {
        const int k0_p = static_cast<int>(hDomains(p)[Dim - 1].first());
        const int k1_p = k0_p + static_cast<int>(hDomains(p)[Dim - 1].length());

        // --- RECV: what k_glob range do I need to receive from p?
        // I need: k_glob_src = Nz_glob - 1 - (k0_r + phys_k_target), phys_k_target in [0, Nz_phys)
        //       = range [Nz_glob - k1_r, Nz_glob - k0_r)
        // Intersect with p's owned range [k0_p, k1_p).
        const int recvA = std::max(k0_p, Nz_glob - k1_r);
        const int recvB = std::min(k1_p, Nz_glob - k0_r);

        if (recvA < recvB) {
            const int recvLen   = recvB - recvA;
            const size_t bytes  = bytesPerKSlab * static_cast<size_t>(recvLen);
            if (p == r) {
                // Self-pair: local copy, no MPI. Source phys_k = k - k0_r,
                // destination phys_k = Nz_glob - 1 - k0_r - k.
                for (int s = 0; s < recvLen; ++s) {
                    const int kGlob      = recvA + s;
                    const int phys_k_src = kGlob - k0_r;
                    const int phys_k_dst = Nz_glob - 1 - k0_r - kGlob;
                    const size_t srcK    = static_cast<size_t>(phys_k_src + nghost);
                    const size_t dstK    = static_cast<size_t>(phys_k_dst + nghost);
                    for (size_t i = 0; i < Xf; ++i) {
                        for (size_t j = 0; j < Yf; ++j) {
                            dstHost(i, j, dstK) = srcHost(i, j, srcK);
                        }
                    }
                }
            } else {
                recvBufs.emplace_back(bytes);
                MPI_Request req;
                MPI_Irecv(recvBufs.back().data(), static_cast<int>(bytes), MPI_BYTE,
                          p, tag, comm, &req);
                requests.push_back(req);
                recvMeta.emplace_back(recvA, recvLen);
            }
        }

        // --- SEND: what of my owned k_glob do I send to p? (skip self)
        // p needs k_glob_src in [Nz_glob - k1_p, Nz_glob - k0_p). Intersect with [k0_r, k1_r).
        if (p != r) {
            const int sendA = std::max(k0_r, Nz_glob - k1_p);
            const int sendB = std::min(k1_r, Nz_glob - k0_p);
            if (sendA < sendB) {
                const int sendLen  = sendB - sendA;
                const size_t bytes = bytesPerKSlab * static_cast<size_t>(sendLen);
                sendBufs.emplace_back(bytes);
                auto* typedPtr = reinterpret_cast<Vector_t<T, Dim>*>(sendBufs.back().data());
                for (int s = 0; s < sendLen; ++s) {
                    const int kGlob      = sendA + s;
                    const int phys_k_src = kGlob - k0_r;
                    const size_t srcK    = static_cast<size_t>(phys_k_src + nghost);
                    for (size_t i = 0; i < Xf; ++i) {
                        for (size_t j = 0; j < Yf; ++j) {
                            *typedPtr++ = srcHost(i, j, srcK);
                        }
                    }
                }
                MPI_Request req;
                MPI_Isend(sendBufs.back().data(), static_cast<int>(bytes), MPI_BYTE,
                          p, tag, comm, &req);
                requests.push_back(req);
            }
        }
    }

    if (!requests.empty()) {
        std::vector<MPI_Status> statuses(requests.size());
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), statuses.data());
    }

    // Unpack received slabs into dstHost.
    for (size_t b = 0; b < recvMeta.size(); ++b) {
        const int recvA   = recvMeta[b].first;
        const int recvLen = recvMeta[b].second;
        const auto* typedPtr =
                reinterpret_cast<const Vector_t<T, Dim>*>(recvBufs[b].data());
        for (int s = 0; s < recvLen; ++s) {
            const int kGlob      = recvA + s;
            const int phys_k_dst = Nz_glob - 1 - k0_r - kGlob;
            const size_t dstK    = static_cast<size_t>(phys_k_dst + nghost);
            for (size_t i = 0; i < Xf; ++i) {
                for (size_t j = 0; j < Yf; ++j) {
                    dstHost(i, j, dstK) = *typedPtr++;
                }
            }
        }
    }

    // Mirror back to device (no-op when device == host).
    Kokkos::deep_copy(flippedZSlab_m, dstHost);
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::gatherFromTempToParticles(
        std::shared_ptr<PartBunch_t> bunch, std::shared_ptr<VField_t<T, Dim>> EtmpSP,
        std::shared_ptr<VField_t<T, Dim>> BtmpSP) {
    // gather accumulated lab-frame E and B from mesh back to particles.
    Inform m("BinnedFieldSolver::gatherFromTempToParticles");
    m << level4 << "gather Etmp/Btmp->particles" << endl;

    VField_t<T, Dim>& Etmp            = *EtmpSP;
    VField_t<T, Dim>& Btmp            = *BtmpSP;
    std::shared_ptr<ParticleCtr_t> pc = bunch->getParticleContainer();

    // gather only the supported field attribute back to particles.
    if (gatherAttribute_m == GatherAttribute::ElectricFieldE) {
        gather(pc->E, Etmp, pc->R);
        gather(pc->B, Btmp, pc->R);
    } else {
        throw OpalException(
                "BinnedFieldSolver::gatherFromTempToParticles",
                "Unsupported gather attribute in binned solver.");
    }
}
