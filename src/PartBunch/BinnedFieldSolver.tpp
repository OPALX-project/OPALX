template <typename T, unsigned Dim>
BinnedFieldSolver<T, Dim>::BinnedFieldSolver(
    std::string solver, Field_t<Dim>* rho, VField_t<T, Dim>* E, Field_t<Dim>* phi,
    std::shared_ptr<BCHandler_t> bcHandler, int tablePrintFrequency)
    : FieldSolver<T, Dim>(solver, rho, E, phi, bcHandler) {
    scatterAttribute_m = ScatterAttribute::ChargeQ;
    gatherAttribute_m  = GatherAttribute::ElectricFieldE;
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
            "BinnedFieldSolver::computeSelfFields", "Bunch particle container is not available.");
    }

    Inform m("BinnedFieldSolver::computeSelfFields");
    // TYPE=NONE is a true no-op: skip all binning/scatter/solve/gather work.
    if (this->getStype() == "NONE") {
        // Already called in ParallelTracker::resetFields()
        // pc->E = 0.0;
        // pc->B = 0.0;
        m << level5 << "Skipping scatter/gather and self-field computation for NONE solver." << endl;
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

    if (hasBins) {
        m << level4 << "Dispatching to computeBinnedSelfFields() (binned path)." << endl;
        computeBinnedSelfFields(bunch);
    } else {
        m << level4 << "Dispatching to computeLegacySelfFields() (legacy path)." << endl;
        computeLegacySelfFields(bunch);
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
void BinnedFieldSolver<T, Dim>::printBinStatsTable(
    const std::string& binningCmdName,
    const std::vector<BinStatsRow>& rows) {
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
                static_cast<long long>(binIndex), static_cast<unsigned long long>(nPartGlobal),
                gammaBin});

        // build rho for this bin and apply lab->solver corrections.
        prepareRhoForBin(bunch, bins, binIndex, nPartGlobal, gammaBin);

        // Ensure deterministic accumulation even for solver types that do not update `E`.
        *(this->getE()) = 0.0;

        m << level4 << "binIndex=" << static_cast<int>(binIndex) << " runSolver(true) start"
          << endl;
        this->runSolver(true);
        m << level4 << "binIndex=" << static_cast<int>(binIndex)
          << " runSolver(true) done; accumulate->Etmp" << endl;

        accumulateFieldToTemp(gammaBin, kinematics.pmean, EtmpSP, BtmpSP);
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
    pc->scaleDtByCharge();
    ippl::ParticleAttrib<T>* dtAttrib = &pc->dt;
    scatter(*dtAttrib, *rho, *R);
    pc->unscaleDtByCharge();

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

    kinematics.pmean = globalPsum / static_cast<double>(nPartGlobal);
    kinematics.gammaBin =
        Kokkos::sqrt(1.0 + kinematics.pmean.dot(kinematics.pmean));
    return kinematics;
}

template <typename T, unsigned Dim>
void BinnedFieldSolver<T, Dim>::prepareRhoForBin(
    std::shared_ptr<PartBunch_t> bunch, std::shared_ptr<AdaptBins_t> bins,
    const bin_index_type binIndex, const size_type nPartGlobal, const double gammaBin) {
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
    pc->scaleDtByCharge();
    ippl::ParticleAttrib<T>* dtAttrib = &pc->dt;
    scatter(*dtAttrib,
            *rho,
            *R,
            bins->getBinIterationPolicy(binIndex),
            bins->getHashArray());
    pc->unscaleDtByCharge();

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
    std::shared_ptr<VField_t<T, Dim>> EtmpSP, std::shared_ptr<VField_t<T, Dim>> BtmpSP) {
    // transform rest-frame fields to lab-frame fields and accumulate.
    Inform m("BinnedFieldSolver::accumulateFieldToTemp");
    m << level4 << "accumulate: gammaBin=" << std::setprecision(10) << gammaBin << endl;

    const double invGamma = (gammaBin > 0.0) ? (1.0 / gammaBin) : 0.0;
    const Vector_t<double, Dim> v = Physics::c * pmean * invGamma;
    const double vNorm             = Kokkos::sqrt(v.dot(v));
    const Vector_t<double, Dim> w =
        (vNorm > 0.0) ? (v / vNorm) : Vector_t<double, Dim>(0.0);
    const double gammaMinusOne = gammaBin - 1.0;
    const double gammaOverCSq  = gammaBin / (Physics::c * Physics::c);

    const VField_t<T, Dim>& Eprime = *(this->getE());
    VField_t<T, Dim>& Etmp         = *EtmpSP;
    VField_t<T, Dim>& Btmp         = *BtmpSP;
    auto ePrimeView                = Eprime.getView();
    auto eTmpView                  = Etmp.getView();
    auto bTmpView                  = Btmp.getView();

    // parallel element-wise transformation and accumulation into temporaries.
    ippl::parallel_for(
        "BinnedFieldSolver::accumulateFieldToTemp", Eprime.getFieldRangePolicy(),
        KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx) {
            Vector_t<T, Dim> ePrime = apply(ePrimeView, idx);
            const T ePrimeDotW = ePrime.dot(w);
            Vector_t<T, Dim> eLab =
                gammaBin * ePrime + gammaMinusOne * ePrimeDotW * w;
            Vector_t<T, Dim> bLab = gammaOverCSq * cross(v, ePrime);
            Vector_t<T, Dim> eTotal = apply(eTmpView, idx);
            Vector_t<T, Dim> bTotal = apply(bTmpView, idx);
            eTotal += eLab;
            bTotal += bLab;
            apply(eTmpView, idx) = eTotal;
            apply(bTmpView, idx) = bTotal;
        });
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
