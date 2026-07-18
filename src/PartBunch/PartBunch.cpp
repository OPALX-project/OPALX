/**
 * @file PartBunch.cpp
 * @brief Template method definitions for PartBunch.
 */

#include "PartBunch/PartBunch.h"
#include <algorithm>
#include <cmath>
#include "Algorithms/Matrix.h"
#include "PartBunch/BinnedFieldSolver.h"
#include "Particle/ParticleAttrib.h"
#include "Physics/ParticleProperties.h"
#include "Structure/Beam.h"
#include "Utilities/Util.h"

#undef doDEBUG

/**
 * @copybrief PartBunch::PartBunch
 */
template <typename T, unsigned Dim>
PartBunch<T, Dim>::PartBunch(
        std::vector<double> qi, std::vector<double> mi, const std::vector<Beam*>& beams,
        std::vector<size_t> totalParticlesPerBeam, double lbt, std::string integration_method,
        FieldSolverCmd* OPALFieldSolver, DataSink* dataSink)
    : ippl::PicManager<
              T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>, LoadBalancer<T, Dim>>(),
      dt_m(0),
      it_m(0),
      integration_method_m(integration_method),
      solver_m(""),
      lbt_m(lbt),
      OPALFieldSolver_m(OPALFieldSolver),
      dataSink_m(dataSink),
      globalTrackStep_m(0),
      rmsDensity_m(0.0) {
    qi_m         = qi;
    mi_m         = mi;
    bunchState_m = std::make_shared<BunchStateHandler>();

    Inform m("PartBunch::PartBunch");
    m << level4 << "PartBunch Constructor" << endl;

    const size_t num_containers = beams.size();
    if (num_containers == 0) {
        throw OpalException("PartBunch::PartBunch", "num_containers must be > 0.");
    }
    if (OPALFieldSolver_m == nullptr) {
        throw OpalException("PartBunch::PartBunch", "OPALFieldSolver must not be null.");
    }
    if (dataSink_m == nullptr) {
        throw OpalException("PartBunch::PartBunch", "dataSink must not be null.");
    }
    if (qi.size() != num_containers) {
        throw OpalException("PartBunch::PartBunch", "qi size must match num_containers.");
    }
    if (mi.size() != num_containers) {
        throw OpalException("PartBunch::PartBunch", "mi size must match num_containers.");
    }
    if (totalParticlesPerBeam.size() != num_containers) {
        throw OpalException(
                "PartBunch::PartBunch", "totalParticlesPerBeam size must match num_containers.");
    }
    for (size_t i = 0; i < num_containers; ++i) {
        if (beams[i] == nullptr) {
            throw OpalException("PartBunch::PartBunch", "beams must not contain null pointers.");
        }
    }

    //  get the needed information from OPAL FieldSolver command

    nr_m = Vector_t<int, Dim>(
            OPALFieldSolver_m->getNX(), OPALFieldSolver_m->getNY(), OPALFieldSolver_m->getNZ());

    const Vector_t<bool, 3> domainDecomposition = OPALFieldSolver_m->getDomainDecomposition();

    for (unsigned i = 0; i < Dim; i++) {
        this->domain_m[i] = ippl::Index(nr_m[i]);
        this->decomp_m[i] = domainDecomposition[i];
    }

    this->setBCHandler(std::make_shared<BCHandler_t>(OPALFieldSolver_m->constructBCHandler()));

    // TODO: support mixed periodic/open per axis; currently all periodic or all open.
    bool isAllPeriodic = this->getBCHandler()->isAll(BCHandler_t::PERIODIC);
    m << level5 << "* FieldContainer set to isAllPeriodic = " << isAllPeriodic << endl;

    //      set stuff for pre_run i.e. warmup
    //      this will be reset when the correct computational
    //      domain is set

    Vector_t<double, Dim> length(6.0);
    this->hr_m     = length / this->nr_m;
    this->origin_m = -3.0;
    this->dt_m     = 0.5 / this->nr_m[2];

    rmin_m = origin_m;
    rmax_m = origin_m + length;

    this->setFieldContainer(
            std::make_shared<FieldContainer_t>(
                    hr_m, rmin_m, rmax_m, decomp_m, domain_m, origin_m, isAllPeriodic));

    this->setParticleContainer(
            std::make_shared<ParticleContainer_t>(
                    this->fcontainer_m->getMesh(), this->fcontainer_m->getFL(),
                    beams[0]->hasPolarization()));
    this->pcontainer_m->setBunchStateHandler(bunchState_m);
    /// \todo if we want, we could also have a separate BunchStateHandler for each container later?
    /// But I think it could also make sense to only have one global handler.
    for (size_t i = 1; i < num_containers; ++i) {
        auto pc = std::make_shared<ParticleContainer_t>(
                this->fcontainer_m->getMesh(), this->fcontainer_m->getFL(),
                beams[i]->hasPolarization());
        pc->setBunchStateHandler(bunchState_m);
        this->addParticleContainer(pc);
    }
    const auto& containers = this->getParticleContainers();
    particleNames_m.resize(containers.size());
    for (size_t i = 0; i < containers.size(); ++i) {
        containers[i]->setQ(qi[i]);
        containers[i]->setM(mi[i]);
        containers[i]->setReference(&beams[i]->getReference());
        particleNames_m[i] = beams[i]->getParticleName();
        containers[i]->Sp  = static_cast<short>(
                ParticleProperties::getParticleType(beams[i]->getParticleName()));
    }

    // Pre-allocate per-rank capacity without bumping localNum_m. Subsequent emission /
    // distribution loaders call createParticles() to fill the buffer non-destructively.
    const double nRanks = static_cast<double>(ippl::Comm->size());
    for (size_t i = 0; i < num_containers; ++i) {
        const size_t maxLocalNum =
                static_cast<size_t>(totalParticlesPerBeam[i] / nRanks + 2 * nRanks + 1);
        containers[i]->allocateParticles(maxLocalNum);
        *gmsg << level3 << "* Container " << i << ": capacity for " << maxLocalNum
              << " particles allocated." << endl;
    }

    setSolver();

    pre_run();
    this->setT(0.0);

    globalPartPerNode_m = std::make_unique<size_t[]>(ippl::Comm->size());

    resetPcActive();

    m << level5 << "* PartBunch constructor done." << endl;
}

/**
 * @copybrief PartBunch::resetPcActive
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::resetPcActive() {
    const auto& containers = this->getParticleContainers();
    const size_t n         = containers.size();
    pcActive_m.resize(n);
    pcAtSStop_m.resize(n);
    for (size_t i = 0; i < n; ++i) {
        pcAtSStop_m[i] = false;
        const auto& pc = containers[i];
        if (!pc || pc->getTotalNum() == 0) {
            pcActive_m[i] = false;
        } else {
            pcActive_m[i] = true;
        }
    }
}

/**
 * @copybrief PartBunch::setPcAtSStop
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::setPcAtSStop(size_t i) {
    if (i >= pcActive_m.size()) {
        return;
    }
    pcActive_m[i]  = false;
    pcAtSStop_m[i] = true;
}

/**
 * @copybrief PartBunch::refreshPcActiveAfterEmit
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::refreshPcActiveAfterEmit() {
    const auto& containers = this->getParticleContainers();
    const size_t n         = containers.size();
    if (pcActive_m.size() != n) {
        return;
    }
    for (size_t i = 0; i < n; ++i) {
        if (pcAtSStop_m[i]) {
            continue;
        }
        const auto& pc = containers[i];
        if (pc && pc->getTotalNum() > 0) {
            pcActive_m[i] = true;
        }
    }
}

/**
 * @copybrief PartBunch::gatherLoadBalanceStatistics
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::gatherLoadBalanceStatistics() {
    std::fill_n(globalPartPerNode_m.get(), ippl::Comm->size(), 0);  // Fill the array with zeros
    globalPartPerNode_m[ippl::Comm->rank()] = this->getParticleContainer()->getLocalNum();
    ippl::Comm->allreduce(globalPartPerNode_m.get(), ippl::Comm->size(), std::plus<size_t>());
}

/**
 * @copybrief PartBunch::setSolver
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::setSolver() {
    Inform m("PartBunch::setSolver");
    m << level2 << "Initializing solver: " << OPALFieldSolver_m->getType() << endl;
    if (this->solver_m != "")
        m << level1 << "Warning solver already initiated but overwrite ..." << endl;

    this->solver_m = OPALFieldSolver_m->getType();

    this->fcontainer_m->initializeFields(this->solver_m);

    auto binnedSolver = std::make_shared<BinnedFieldSolver<T, Dim>>(
            this->solver_m, this->fcontainer_m->sharedWorkspace(), this->getBCHandler(),
            OPALFieldSolver_m->getGreensFunction());
    this->setFieldSolver(binnedSolver);
    m << level4 << "Binned field solver set (binned or legacy at runtime)." << endl;

    this->fsolver_m->initSolver();
    m << level4 << "Field solver initialized." << endl;

    this->setLoadBalancer(std::make_shared<LoadBalancer_t>(this->lbt_m));
    m << level3 << "Solver and Load Balancer set." << endl;
}

template <typename T, unsigned Dim>
typename PartBunch<T, Dim>::BinnedFieldSolver_t* PartBunch<T, Dim>::getFieldSolver() {
    return static_cast<BinnedFieldSolver_t*>(this->fsolver_m.get());
}

template <typename T, unsigned Dim>
const typename PartBunch<T, Dim>::BinnedFieldSolver_t* PartBunch<T, Dim>::getFieldSolver() const {
    return static_cast<const BinnedFieldSolver_t*>(this->fsolver_m.get());
}

template <typename T, unsigned Dim>
std::string PartBunch<T, Dim>::getFieldSolverType() {
    return this->getFieldSolver()->getStype();
}

/**
 * @copybrief PartBunch::calcBeamParameters
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::calcBeamParameters() {
    Inform m("PartBunch::calcBeamParameters");
    std::shared_ptr<ParticleContainer_t> pc = this->pcontainer_m;

    using view_type = ippl::ParticleAttrib<Vector_t<double, 3>>::view_type;
    view_type Rview = pc->R.getView();
    view_type Pview = pc->P.getView();
    this->getParticleContainer()->updateMoments();
    m << level5 << "Moments updated." << endl;

    ////////////////////////////////////
    //// Calculate Moments of R and P //
    ////////////////////////////////////

    using MomentsVec = ippl::Vector<double, 2 * Dim>;
    using MomentsMat = ippl::Vector<MomentsVec, 2 * Dim>;

    MomentsVec loc_centroid(0.0);
    MomentsMat loc_moment(MomentsVec(0.0));

    MomentsVec centroid(0.0);
    MomentsMat moment(MomentsVec(0.0));

    for (unsigned i = 0; i < 2 * Dim; ++i) {
        const size_t nLocal = pc->getLocalNum();
        Kokkos::parallel_reduce(
                "calc moments of particle distr.", nLocal,
                KOKKOS_LAMBDA(
                        const size_t k, double& cent, double& mom0, double& mom1, double& mom2,
                        double& mom3, double& mom4, double& mom5) {
                    double part[2 * Dim];
                    part[0] = Rview(k)[0];
                    part[1] = Pview(k)[0];
                    part[2] = Rview(k)[1];
                    part[3] = Pview(k)[1];
                    part[4] = Rview(k)[2];
                    part[5] = Pview(k)[2];

                    cent += part[i];
                    mom0 += part[i] * part[0];
                    mom1 += part[i] * part[1];
                    mom2 += part[i] * part[2];
                    mom3 += part[i] * part[3];
                    mom4 += part[i] * part[4];
                    mom5 += part[i] * part[5];
                },
                Kokkos::Sum<T>(loc_centroid[i]), Kokkos::Sum<T>(loc_moment[i][0]),
                Kokkos::Sum<T>(loc_moment[i][1]), Kokkos::Sum<T>(loc_moment[i][2]),
                Kokkos::Sum<T>(loc_moment[i][3]), Kokkos::Sum<T>(loc_moment[i][4]),
                Kokkos::Sum<T>(loc_moment[i][5]));
        Kokkos::fence();
    }
    m << level5 << "Local moments calculated." << endl;

    moment   = loc_moment;
    centroid = loc_centroid;
    ippl::Comm->allreduce(moment, 1, std::plus<MomentsMat>());
    ippl::Comm->allreduce(centroid, 1, std::plus<MomentsVec>());
    ippl::Comm->barrier();
    m << level5 << "Global moments calculated." << endl;

    ippl::Vector<double, Dim> rmax_loc(0.0);
    ippl::Vector<double, Dim> rmin_loc(0.0);
    ippl::Vector<double, Dim> rmax(0.0);
    ippl::Vector<double, Dim> rmin(0.0);

    // TODO: fuse min/max reductions with ippl::Vector reductions.
    for (unsigned d = 0; d < Dim; ++d) {
        Kokkos::parallel_reduce(
                "rel max", pc->getLocalNum(),
                KOKKOS_LAMBDA(const int i, double& mm) {
                    double tmp_vel = Rview(i)[d];
                    mm             = tmp_vel > mm ? tmp_vel : mm;
                },
                Kokkos::Max<T>(rmax_loc[d]));

        Kokkos::parallel_reduce(
                "rel min", pc->getLocalNum(),
                KOKKOS_LAMBDA(const int i, double& mm) {
                    double tmp_vel = Rview(i)[d];
                    mm             = tmp_vel < mm ? tmp_vel : mm;
                },
                Kokkos::Min<T>(rmin_loc[d]));
    }
    m << level5 << "Local min/max calculated." << endl;
    Kokkos::fence();
    rmax = rmax_loc;
    rmin = rmin_loc;
    ippl::Comm->allreduce(rmax, 1, std::greater<ippl::Vector<double, Dim>>());
    ippl::Comm->allreduce(rmin, 1, std::less<ippl::Vector<double, Dim>>());
    ippl::Comm->barrier();
    m << level5 << "Global min/max calculated." << endl;

    rmax_m = rmax;
    rmin_m = rmin;
}

/**
 * @copybrief PartBunch::pre_run
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::pre_run() {
    Inform m("PartBunch::pre_run");
    m << level2 << "Starting pre_run..." << endl;
    auto rhoView = this->fcontainer_m->getRho().getView();
    Kokkos::deep_copy(rhoView, 0.0);
    m << level4 << "Rho initialized to zero." << endl;

    /*
     * Skip full field dumps during warmup: runSolver(true) is implemented on the
     * concrete solver type, not on the IPPL base class.
     */
    this->getFieldSolver()->runSolver(true);
    m << level4 << "Field solver ran during pre_run." << endl;
    m << level4 << "Warm-up solve excluded from runtime diagnostics. pre_run done." << endl;
}

/**
 * @copybrief PartBunch::print
 */
template <typename T, unsigned Dim>
Inform& PartBunch<T, Dim>::print(Inform& os) {
    // if (pc->getLocalNum() != 0) {  // to suppress Nans
    Inform::FmtFlags_t ff = os.flags();

    const auto& containers = this->getParticleContainers();
    for (size_t ci = 0; ci < containers.size(); ++ci) {
        const auto& pc = containers[ci];
        if (!pc) {
            os << level1 << "Skipping null particle container: " << ci << endl;
            continue;
        }

        const double ek  = pc->getMeanKineticEnergy();
        const double dek = pc->getStdKineticEnergy();

        // ParticleContainer tracks charge/mass storage mode for QM attributes.
        std::string qmStorageModeStr = "SINGLE";
        const auto qmMode            = pc->getQMStorageMode();
        if (qmMode == ParticleContainer_t::QMStorageMode::Attributes) {
            qmStorageModeStr = "ATTRIBUTES";
        }

        os << level1 << std::scientific << "\n"
           << "* ************** B U N C H "
              "********************************************************* \n"
           << "* CONTAINER       = " << ci << "\n"
           << "* PARTICLES       = " << pc->getTotalNum() << "\n"
           << "* CHARGE          = " << pc->getTotalCharge() << " (Cb) \n"
           << "* QM STORAGE MODE = " << qmStorageModeStr << "\n"
           << "* <EKIN>          = " << Util::getEnergyString(ek) << "\n"
           << "* <dEKIN>         = " << Util::getEnergyString(dek) << "\n"
           << "* INTEGRATOR      = " << integration_method_m << "\n"
           << "* LB Threshold    = " << lbt_m << "\n"
           << "* MIN R (origin)  = " << Util::getLengthString(pc->getMinR(), 5) << "\n"
           << "* MAX R (max ext) = " << Util::getLengthString(pc->getMaxR(), 5) << "\n"
           << "* RMS R           = " << Util::getLengthString(pc->getRmsR(), 5) << "\n"
           << "* RMS P           = " << pc->getRmsP() << " [beta gamma]\n"
           << "* Mean R          = " << pc->getMeanR() << " [m]\n"
           << "* Mean P          = " << pc->getMeanP() << " [beta gamma]\n"
           << "* MESH SPACING    = "
           << Util::getLengthString(this->fcontainer_m->getMesh().getMeshSpacing(), 5) << "\n"
           << "* COMPDOM INCR    = " << this->OPALFieldSolver_m->getBoxIncr() << " (%) \n"
           << "* FIELD LAYOUT    = " << this->fcontainer_m->getFL() << "\n"
           << "* Centroid : \n* ";
        for (unsigned int i = 0; i < 2 * Dim; i++) {
            os << level1 << pc->getCentroid()[i] << " ";
        }
        os << level1 << endl << "* Cov Matrix : \n* ";
        for (unsigned int i = 0; i < 2 * Dim; i++) {
            for (unsigned int j = 0; j < 2 * Dim; j++) {
                os << level1 << pc->getCovMatrix()(i, j) << " ";
            }
            os << level1 << "\n* ";
        }
        os << level1
           << "* "
              "********************************************************************************"
              "** \n"
           << endl;
    }

    os.flags(ff);
    return os;
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::updateAllParticleMoments() {
    const auto& containers = this->getParticleContainers();
    for (size_t i = 0; i < containers.size(); ++i) {
        if (!containers[i]) {
            continue;
        }
        this->getParticleContainer(i)->updateMoments();
    }
}

/**
 * @copybrief PartBunch::computeSelfFields
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::computeSelfFields(
        opalx::spacecharge::IterationPlan<T, Dim>& iterationPlan, std::uint64_t particleGeneration,
        const opalx::spacecharge::PreparedCorrection<T, Dim>& correction,
        opalx::spacecharge::BinConfigurationObserver* binConfigurationObserver,
        opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
    BinnedFieldSolver_t* bsolver = this->getFieldSolver();

    bsolver->computeSelfFields(
            *this, iterationPlan, particleGeneration, correction, binConfigurationObserver,
            diagnostics);
}

template <typename T, unsigned Dim>
int PartBunch<T, Dim>::getCurrentNBins() const {
    return this->getFieldSolver()->legacyReportedBinCount();
}

/**
 * @copybrief PartBunch::performBunchSanityChecks
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::performBunchSanityChecks() const {
    Inform ms("PartBunch::performBunchSanityChecks");
    ms << level4 << "========== Performing sanity checks on PartBunch... ==========" << endl;
    // TODO: extend checks; prefer throwing OpalException with clear messages.

    // Check if bc handler was initialized properly
    if (!this->getBCHandler()) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks", "BC Handler not initialized properly.");
    }
    ms << level4 << "BC Handler initialized properly." << endl;

    if (!this->getBunchStateHandler()) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks", "BunchStateHandler not initialized.");
    }
    ms << level4 << "BunchStateHandler initialized." << endl;

    if (!hasFieldSolver()) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks", "Field Solver was not initialized.");
    }
    ms << level4 << "Field Solver object was initialized." << endl;

    // Verify we can access the concrete FieldSolver and its internals.
    auto fs = std::dynamic_pointer_cast<BinnedFieldSolver_t>(this->fsolver_m);
    if (!fs) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks", "FieldSolver is not set in PartBunch.");
    }

    // cannot use getFieldContainer, since this getter cannot be const!
    const std::shared_ptr<FieldContainer<T, Dim>> fctr = this->fcontainer_m;
    if (!fctr) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks",
                "FieldContainer isn't initialized correctly.");
    }

    if (&fs->getWorkspace() != &fctr->workspace()) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks",
                "FieldSolver and FieldContainer do not share the same PIC workspace.");
    }
    ms << level4 << "FieldSolver and FieldContainer share one PIC workspace." << endl;

    // Check internal field pointers are set
    if (fs->getRho() == nullptr || fs->getE() == nullptr || fs->getPhi() == nullptr) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks",
                "FieldSolver internal fields (rho/E/phi) not assigned.");
    }
    ms << level4 << "FieldSolver internal field pointers are set." << endl;

    // Ensure FieldSolver fields point to our FieldContainer's fields
    if (fs->getRho() != &fctr->getRho() || fs->getE() != &fctr->getE()
        || fs->getPhi() != &fctr->getPhi()) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks",
                "FieldSolver fields do not match FieldContainer.");
    }
    ms << level4 << "FieldSolver fields match FieldContainer." << endl;

    /*
    // Check if all three fields (rho, E, phi) have the same mesh and layout
    auto rhoMesh = fs->getRho()->get_mesh();
    auto EMesh   = fs->getE()->get_mesh();
    auto phiMesh = fs->getPhi()->get_mesh();
    if (rhoMesh->getOrigin() != EMesh->getOrigin() ||
        rhoMesh->getOrigin() != phiMesh->getOrigin() ||
        rhoMesh->getMeshSpacing() != EMesh->getMeshSpacing() ||
        rhoMesh->getMeshSpacing() != phiMesh->getMeshSpacing()) {
        throw OpalException("PartBunch::performBunchSanityChecks",
                            "FieldSolver fields do not share the same mesh.");
    }
    ms << "FieldSolver fields share the same mesh." << endl;*/

    // Check solver type string and that a backend was emplaced
    const std::string stype = fs->getStype();
    if (stype.empty()) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks", "FieldSolver type string is empty.");
    }
    if (stype != "FFT" && stype != "OPEN" && stype != "CG" && stype != "NONE") {
        throw OpalException(
                "PartBunch::performBunchSanityChecks", "Unsupported FieldSolver type: " + stype);
    }
    ms << level4 << "FieldSolver type: " << stype << endl;

    // Basic check that the E-field layout has non-zero extent
    auto Eview = fctr->getE().getView();
    if (stype != "NONE" && (Eview.extent(0) == 0 || Eview.extent(1) == 0 || Eview.extent(2) == 0)) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks",
                "E-field layout not initialized (zero extent). ");
    }
    ms << level4 << "E-field layout initialized." << endl;

    // Temporary E/B accumulation fields (binned solver path)
    auto Etmp = fctr->getTempEField();
    auto Btmp = fctr->getTempBField();
    if (!Etmp || !Btmp) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks",
                "Temporary E field (Etmp) and/or B field (Btmp) not initialized.");
    }
    auto EtmpView = Etmp->getView();
    auto BtmpView = Btmp->getView();
    if (EtmpView.extent(0) == 0 || EtmpView.extent(1) == 0 || EtmpView.extent(2) == 0) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks",
                "Etmp field layout not initialized (zero extent). ");
    }
    if (BtmpView.extent(0) == 0 || BtmpView.extent(1) == 0 || BtmpView.extent(2) == 0) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks",
                "Btmp field layout not initialized (zero extent). ");
    }
    if (&Etmp->get_mesh() != &fctr->getMesh() || &Btmp->get_mesh() != &fctr->getMesh()) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks",
                "Etmp/Btmp fields do not use the FieldContainer mesh.");
    }
    auto mirrorScratch = fctr->getFlippedZSlabField();
    if (!mirrorScratch || &mirrorScratch->get_mesh() != &fctr->getMesh()
        || &mirrorScratch->getLayout() != &fctr->getFL()) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks",
                "Persistent mirror scratch does not use the PIC workspace mesh and layout.");
    }
    ms << level4 << "Persistent scratch fields use the PIC workspace mesh and layout." << endl;

    if (!this->pcontainer_m) {
        throw OpalException(
                "PartBunch::performBunchSanityChecks",
                "Primary ParticleContainer not initialized.");
    }
    ms << level4 << "Primary ParticleContainer present." << endl;

    ms << level2 << "========= Done performing PartBunch sanity checks... =========" << endl;
}

/** Explicit instantiation for 3D double (OPAL-T). */
template class PartBunch<double, 3>;
