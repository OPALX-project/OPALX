#include "Structure/DataSink.h"

#include "PartBunch/FieldOps.hpp"

#include <cstring>
#include <utility>
#include <vector>

template <typename T, unsigned Dim>
BinnedFieldSolver<T, Dim>::BinnedFieldSolver(
        std::string solver, Field_t<Dim>* rho, VField_t<T, Dim>* E, Field_t<Dim>* phi,
        std::shared_ptr<BCHandler_t> bcHandler, int tablePrintFrequency, bool adaptiveBinning)
    : FieldSolver<T, Dim>(solver, rho, E, phi, bcHandler) {
    scatterAttribute_m    = ScatterAttribute::ChargeQ;
    gatherAttribute_m     = GatherAttribute::ElectricFieldE;
    tablePrintFrequency_m = tablePrintFrequency;
    adaptiveBinning_m     = adaptiveBinning;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::computeSelfFields(PartBunch_t& bunch) {
    // Validate inputs and decide between binned vs legacy solver.
    std::shared_ptr<ParticleCtr_t> pc = bunch.getParticleContainer();
    if (!pc) {
        throw OpalException(
                "BinnedFieldSolver::computeSelfFields",
                "Bunch particle container is not available.");
    }

    Inform m("BinnedFieldSolver::computeSelfFields");
    // TYPE=NONE is a true no-op: skip all binning/scatter/solve/gather work.
    if (this->getSolverCapabilities().isNoOp) {
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

    const size_t step = static_cast<size_t>(bunch.getGlobalTrackStep());
    const BoundaryCorrection<T, Dim>* activeCorrection = activeBoundaryCorrectionForStep(step);
    const BoundaryCorrection<T, Dim>* configuredCorrection = configuredBoundaryCorrection();
    if (configuredCorrection && !activeCorrection) {
        m << level3 << "ZEROFACE_MAXSTEPS reached (step=" << bunch.getGlobalTrackStep()
          << ", maxSteps=" << zerofaceMaxSteps_m << "); " << configuredCorrection->name()
          << " correction inactive for this step." << endl;
    }

    if (hasBins) {
        m << level4 << "Dispatching to computeBinnedSelfFields() (binned path)." << endl;
        computeBinnedSelfFields(bunch, activeCorrection);
    } else {
        // Legacy path has no separate correction pass: it scatters primary+image
        // in one shot via ImageChargeScatterController::scatterPrimaryAndImage
        // and does one standard solve. The shifted Green's correction is only
        // implemented for the binned path, so warn once if the user requested it
        // without binning.
        if (activeCorrection && activeCorrection->kind() == BoundaryCorrectionKind::ShiftedGreens) {
            m << level3 << "SHIFTED_GREENS_FUNCTION is set but no binning is active; "
              << "the legacy path does not apply the Dirichlet correction." << endl;
        }
        m << level4 << "Dispatching to computeLegacySelfFields() (legacy path)." << endl;
        computeLegacySelfFields(bunch, activeCorrection);
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
    if (enabled && shiftedGreensCorrection_m.isEnabled()) {
        throw OpalException(
                "BinnedFieldSolver::setImageChargeConfiguration",
                "Cannot enable image charges while SHIFTED_GREENS_FUNCTION is active: "
                "ZEROFACE_R0Z and SHIFTED_GREENS_FUNCTION are mutually exclusive.");
    }
    imageCorrection_m.configure(enabled, zPlane);
    imageScatterController_m.configure(enabled, zPlane);
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::setShiftedGreensConfiguration(bool enabled, double zPlane) {
    if (enabled && imageCorrection_m.isEnabled()) {
        throw OpalException(
                "BinnedFieldSolver::setShiftedGreensConfiguration",
                "Cannot enable SHIFTED_GREENS_FUNCTION while image charges are active: "
                "ZEROFACE_R0Z and SHIFTED_GREENS_FUNCTION are mutually exclusive.");
    }
    shiftedGreensCorrection_m.configure(enabled, zPlane, this->getSolverCapabilities());
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
    return imageCorrection_m.isActiveForStep(step, zerofaceMaxSteps_m);
}

template <typename T, unsigned Dim>
bool BinnedFieldSolver<T, Dim>::isShiftedGreensActiveForStep(size_t step) const {
    return shiftedGreensCorrection_m.isActiveForStep(step, zerofaceMaxSteps_m);
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
        PartBunch_t& bunch, const std::string& solveTag) {
    if (!imageCorrection_m.isEnabled() || zeroFacePlaneDumpFrequency_m <= 0) {
        return;
    }

    const long long step = bunch.getGlobalTrackStep();
    if (step < 0 || (step % zeroFacePlaneDumpFrequency_m) != 0) {
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

    Field_t<Dim>* potentialField =
            this->getSolverCapabilities().usesSeparatePotentialField ? this->getPhi()
                                                                     : this->getRho();
    if (!potentialField) {
        return;
    }
    const double zPlane = imageCorrection_m.planeZ();

    DataSink* dataSink = bunch.getDataSink();
    if (!dataSink) {
        return;
    }

    const auto diagnostics =
            dataSink->dumpDirichletPlane(step, bunch.getT(), zPlane, *potentialField, solveTag);
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
const BoundaryCorrection<T, Dim>* BinnedFieldSolver<T, Dim>::configuredBoundaryCorrection() const {
    if (imageCorrection_m.isEnabled()) {
        return &imageCorrection_m;
    }
    if (shiftedGreensCorrection_m.isEnabled()) {
        return &shiftedGreensCorrection_m;
    }
    return nullptr;
}

template <typename T, unsigned Dim>
const BoundaryCorrection<T, Dim>* BinnedFieldSolver<T, Dim>::activeBoundaryCorrectionForStep(
        size_t step) const {
    const BoundaryCorrection<T, Dim>* correction = configuredBoundaryCorrection();
    if (!correction || !correction->isActiveForStep(step, zerofaceMaxSteps_m)) {
        return nullptr;
    }
    return correction;
}

template <typename T, unsigned Dim>
ImageScatterMode BinnedFieldSolver<T, Dim>::primaryScatterModeForStep(
        const BoundaryCorrection<T, Dim>* activeCorrection) const {
    if (activeCorrection) {
        return ImageScatterMode::PrimaryOnly;
    }

    // If explicit image charges are configured but outside the step budget, avoid
    // scatterPrimaryAndImage because the controller itself remains configured.
    if (imageCorrection_m.isEnabled()) {
        return ImageScatterMode::PrimaryOnly;
    }

    return ImageScatterMode::PrimaryAndImage;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::applyBoundaryCorrectionPass(
        PartBunch_t& bunch, std::shared_ptr<AdaptBins_t> bins, const BinContext<T, Dim>& context,
        const BoundaryCorrectionPass<T, Dim>& pass, bool& dumpedDirichletPlaneThisStep) {
    Inform m("BinnedFieldSolver::applyBoundaryCorrectionPass");

    auto& mesh = this->getRho()->get_mesh();
    prepareRhoForBin(bunch, bins, context, pass.scatterMode);

    opalx::fieldops::setVectorField(*(this->getE()), Vector_t<T, Dim>(0.0));
    mesh.setMeshSpacing(context.solveSpacing);

    m << level4 << "binIndex=" << context.binIndex << " correction runSolver start";
    if (pass.solveRequest.hasShiftedGreens()) {
        m << ", greensShift=" << *pass.solveRequest.greensShift;
    }
    m << endl;

    this->runSolver(pass.solveRequest, pass.forceSkipFieldDump);

    m << level4 << "binIndex=" << context.binIndex
      << " correction runSolver done; accumulate->Etmp (B sign=" << pass.bFieldSign
      << ", flipAxis=" << pass.flipAxis << ")" << endl;

    fieldAccumulator_m.accumulate(
            *(this->getE()), context.gamma, context.pmean, pass.bFieldSign, pass.flipAxis);

    mesh.setMeshSpacing(context.meshSpacing);

    if (pass.dumpDirichletPlane && !dumpedDirichletPlaneThisStep) {
        dumpDirichletPlaneDiagnosticsIfRequested(bunch, "binned");
        dumpedDirichletPlaneThisStep = true;
    }
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::computeBinnedSelfFields(
        PartBunch_t& bunch, const BoundaryCorrection<T, Dim>* activeCorrection) {
    // execute full binned self-field algorithm.
    // fetch the adaptive bin structure.
    std::shared_ptr<AdaptBins_t> bins = bunch.getBins();
    if (!bins) {
        // Defensive: runtime selection above should prevent this.
        computeLegacySelfFields(bunch, activeCorrection);
        return;
    }

    // build and merge adaptive bins for this step.
    // build and merge adaptive bins for this step.
    rebinAndPrepare(bunch, bins);

    // obtain the temporary E buffer used to accumulate bin contributions.
    std::shared_ptr<VField_t<T, Dim>> EtmpSP = bunch.getTempEField();
    if (!EtmpSP) {
        throw OpalException(
                "BinnedFieldSolver::computeBinnedSelfFields",
                "Temporary E field (Etmp) is not initialized.");
    }
    std::shared_ptr<VField_t<T, Dim>> BtmpSP = bunch.getTempBField();
    if (!BtmpSP) {
        throw OpalException(
                "BinnedFieldSolver::computeBinnedSelfFields",
                "Temporary B field (Btmp) is not initialized.");
    }

    fieldAccumulator_m.bind(EtmpSP, BtmpSP);
    fieldAccumulator_m.clear();

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
        if (kinematics.gammaBin <= 0.0) {
            throw OpalException(
                    "BinnedFieldSolver::computeBinnedSelfFields",
                    "Computed non-positive gamma for bin.");
        }

        m << level4 << "binIndex=" << static_cast<int>(binIndex)
          << " gammaBin=" << std::setprecision(10) << kinematics.gammaBin << endl;

        binStats.push_back(
                BinStatsRow{
                        static_cast<long long>(binIndex),
                        static_cast<unsigned long long>(nPartGlobal), kinematics.gammaBin});

        // Mesh references reused for both primary and correction passes.
        auto& mesh        = this->getRho()->get_mesh();
        const auto hrOrig = mesh.getMeshSpacing();
        auto hrStretched  = hrOrig;
        hrStretched[Dim - 1] *= kinematics.gammaBin;

        BinContext<T, Dim> context;
        context.binIndex                     = static_cast<long long>(binIndex);
        context.particleCountGlobal          = static_cast<std::uint64_t>(nPartGlobal);
        context.gamma                        = kinematics.gammaBin;
        context.pmean                        = kinematics.pmean;
        context.meshSpacing                  = hrOrig;
        context.solveSpacing                 = hrStretched;
        context.imageCorrectionActive =
                activeCorrection
                && activeCorrection->kind() == BoundaryCorrectionKind::ImageCharge;
        context.shiftedGreensCorrectionActive =
                activeCorrection
                && activeCorrection->kind() == BoundaryCorrectionKind::ShiftedGreens;

        // --- Primary pass: scatter real charges, solve, accumulate with +B ---
        {
            const ImageScatterMode scatterMode = primaryScatterModeForStep(activeCorrection);
            prepareRhoForBin(bunch, bins, context, scatterMode);

            opalx::fieldops::setVectorField(*(this->getE()), Vector_t<T, Dim>(0.0));
            mesh.setMeshSpacing(context.solveSpacing);

            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " primary runSolver(true) start"
              << " (hr_z stretched by gamma=" << context.gamma << ")" << endl;
            this->runSolver(true);
            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " primary runSolver(true) done; accumulate->Etmp" << endl;

            fieldAccumulator_m.accumulate(*(this->getE()), context.gamma, context.pmean, +1.0);

            mesh.setMeshSpacing(context.meshSpacing);
        }

        if (activeCorrection) {
            const auto correctionPass = activeCorrection->makePass(
                    context, mesh.getOrigin(),
                    this->getRho()->getLayout().getDomain()[Dim - 1].length());
            if (correctionPass) {
                applyBoundaryCorrectionPass(
                        bunch, bins, context, *correctionPass, dumpedDirichletPlaneThisStep);
            }
        }
    }

    // after all bins, gather the accumulated lab-frame field back to particles.
    gatherFromTempToParticles(bunch);

    // per-call table: gammaBin / nParticles / binNumber.
    if (tablePrintFrequency_m > 0) {
        const long long step = bunch.getGlobalTrackStep();
        if (step >= 0 && (step % tablePrintFrequency_m) == 0) {
            printBinStatsTable(bins->getBinningCmdName(), binStats);
        }
    }
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::computeLegacySelfFields(
        PartBunch_t& bunch, const BoundaryCorrection<T, Dim>* activeCorrection) {
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

    typename PartBunch_t::Base::particle_position_type* R = &pc->R;

    Field_t<Dim>& rho = *(this->getRho());
    opalx::fieldops::setScalarField(rho, 0.0);

    // Scatter charge to mesh rho using dt-weighted deposition (master approach):
    // scale dt by Q, scatter dt, then restore dt.
    if (activeCorrection && activeCorrection->kind() == BoundaryCorrectionKind::ImageCharge) {
        imageScatterController_m.scatterPrimaryAndImage(pc, *R, rho);
    } else {
        imageScatterController_m.scatterPrimaryOnly(pc, *R, rho);
    }

    //  apply mesh normalization, background subtraction, and rho scaling.
    const std::string stype = this->getStype();
    double normalizer       = bunch.getdT();
    if (stype != "FEM" && stype != "FEM_PRECON") {
        const double cellVolume =
                std::reduce(bunch.hr_m.begin(), bunch.hr_m.end(), 1.0, std::multiplies<double>());
        normalizer *= cellVolume;
    }

    // Alpine uses net-0 charge for non-OPEN solvers (periodic BCs).
    double shift = 0.0;
    if (stype != "OPEN") {
        double size = 1.0;
        for (size_t d = 0; d < Dim; ++d) {
            size *= bunch.rmax_m[d] - bunch.rmin_m[d];
        }

        const double totalQ = bunch.getParticleContainer()->getTotalCharge();
        shift               = -(totalQ / size) * this->getCouplingConstant();
    }

    opalx::fieldops::scaleAndShiftScalarField(
            rho, this->getCouplingConstant() / normalizer, shift);

    // Ensure deterministic output even for solver types that do not update `E`.
    opalx::fieldops::setVectorField(*(this->getE()), Vector_t<T, Dim>(0.0));

    // run the solver once and gather mesh E back to particles.
    m << level4 << "Legacy mode: runSolver() start" << endl;
    this->runSolver();
    if (activeCorrection && activeCorrection->kind() == BoundaryCorrectionKind::ImageCharge) {
        dumpDirichletPlaneDiagnosticsIfRequested(bunch, "legacy");
    }
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
        PartBunch_t& bunch, std::shared_ptr<AdaptBins_t> bins,
        const BinContext<T, Dim>& context, ImageScatterMode mode) {
    // Scatter bin charge to rho using dt-weighted deposition.
    // If the ParticleContainer supports scaleDtByCharge(), use the master approach:
    // scale dt by charge, scatter dt, then unscale.
    // Otherwise, fall back to temporarily scaling/scattering Q by dt.
    Inform m("BinnedFieldSolver::prepareRhoForBin");
    m << level4 << "prepareRho: binIndex=" << context.binIndex
      << ", nPartGlobal=" << context.particleCountGlobal
      << ", gammaBin=" << std::setprecision(10) << context.gamma << endl;

    Field_t<Dim>& rho = *(this->getRho());
    opalx::fieldops::setScalarField(rho, 0.0);

    // access particle views and validate scatter support.
    std::shared_ptr<ParticleCtr_t> pc                     = bunch.getParticleContainer();
    typename PartBunch_t::Base::particle_position_type* R = &pc->R;

    if (scatterAttribute_m != ScatterAttribute::ChargeQ) {
        throw OpalException(
                "BinnedFieldSolver::prepareRhoForBin",
                "Unsupported scatter attribute in binned solver.");
    }

    // Scatter bin charge to rho (with bin iteration policy and hash indexing).
    // Master approach: scale dt by Q, scatter dt, then restore dt.
    const auto binIndex    = static_cast<bin_index_type>(context.binIndex);
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

    if (mode == ImageScatterMode::PrimaryOnly) {
        if (scatterAllLocal) {
            imageScatterController_m.scatterPrimaryOnly(pc, *R, rho);
        } else {
            imageScatterController_m.scatterPrimaryOnly(pc, *R, rho, policy, hash);
        }
    } else if (mode == ImageScatterMode::ImageOnly) {
        if (scatterAllLocal) {
            imageScatterController_m.scatterImageOnly(pc, *R, rho);
        } else {
            imageScatterController_m.scatterImageOnly(pc, *R, rho, policy, hash);
        }
    } else {
        if (scatterAllLocal) {
            imageScatterController_m.scatterPrimaryAndImage(pc, *R, rho);
        } else {
            imageScatterController_m.scatterPrimaryAndImage(pc, *R, rho, policy, hash);
        }
    }

    // normalize rho for fractional time steps and mesh conventions.
    const std::string stype = this->getStype();
    double normalizer       = bunch.getdT();
    if (stype != "FEM" && stype != "FEM_PRECON") {
        const double cellVolume =
                std::reduce(bunch.hr_m.begin(), bunch.hr_m.end(), 1.0, std::multiplies<double>());
        normalizer *= cellVolume;
    }

    // subtract non-OPEN background and apply Lorentz rest-frame scaling.
    // Background subtraction for non-OPEN solvers. Here we subtract only the bin's charge.
    double shift = 0.0;
    if (stype != "OPEN") {
        double size = 1.0;
        for (size_t d = 0; d < Dim; ++d) {
            size *= bunch.rmax_m[d] - bunch.rmin_m[d];
        }

        const double totalQBin = bunch.getParticleContainer()->getChargePerParticle()
                                 * static_cast<double>(context.particleCountGlobal);
        shift = -(totalQBin / size) * this->getCouplingConstant() / context.gamma;
    }

    // Lorentz transform of charge density to the bin rest frame (thesis Eq. step 7).
    normalizer *= context.gamma;
    opalx::fieldops::scaleAndShiftScalarField(
            rho, this->getCouplingConstant() / normalizer, shift);
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::gatherFromTempToParticles(PartBunch_t& bunch) {
    // gather accumulated lab-frame E and B from mesh back to particles.
    Inform m("BinnedFieldSolver::gatherFromTempToParticles");
    m << level4 << "gather Etmp/Btmp->particles" << endl;

    std::shared_ptr<ParticleCtr_t> pc = bunch.getParticleContainer();

    // gather only the supported field attribute back to particles.
    if (gatherAttribute_m == GatherAttribute::ElectricFieldE) {
        fieldAccumulator_m.gatherToParticles(*pc);
    } else {
        throw OpalException(
                "BinnedFieldSolver::gatherFromTempToParticles",
                "Unsupported gather attribute in binned solver.");
    }
}
