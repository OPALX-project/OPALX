#include "PartBunch/PartBunch.h"
#include "PartBunch/BinnedFieldSolver.h"
#include "Algorithms/Matrix.h"
#include "Particle/ParticleAttrib.h"
#include "Utilities/Util.h"
#include "Structure/DataSink.h"

#undef doDEBUG

template <typename T, unsigned Dim>
PartBunch<T, Dim>::PartBunch(std::vector<double> qi,
                             std::vector<double> mi,
                             size_t num_containers,
                             /*int nt,*/
                             double lbt,
                             std::string integration_method,
                             std::shared_ptr<FieldSolverCmd> OPALFieldSolver,
                             std::shared_ptr<DataSink> dataSink)
    : ippl::PicManager<T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>, LoadBalancer<T, Dim>>(),
      RefPartR_m(0.0),
      RefPartP_m(0.0),
      dt_m(0),
      it_m(0),
      lbt_m(lbt),
      isFirstRepartition_m(true),
      integration_method_m(integration_method),
      solver_m(""),
      //nt_m(nt),
      qi_m(std::move(qi)),
      mi_m(std::move(mi)),
      OPALFieldSolver_m(OPALFieldSolver),
      dataSink_m(std::move(dataSink)),
      globalTrackStep_m(0),
      rmsDensity_m(0.0) {

    Inform m("PartBunch::PartBunch");
    m << level4 << "PartBunch Constructor" << endl;

    if (num_containers == 0) {
        throw OpalException("PartBunch::PartBunch",
                            "num_containers must be > 0.");
    }
    if (qi_m.size() != num_containers) {
        throw OpalException("PartBunch::PartBunch",
                            "qi size must match num_containers.");
    }
    if (mi_m.size() != num_containers) {
        throw OpalException("PartBunch::PartBunch",
                            "mi size must match num_containers.");
    }

    //  get the needed information from OPAL FieldSolver command

    nr_m = Vector_t<int, Dim>(
            OPALFieldSolver_m->getNX(), 
            OPALFieldSolver_m->getNY(), 
            OPALFieldSolver_m->getNZ()
    );

    const Vector_t<bool, 3> domainDecomposition = OPALFieldSolver_m->getDomDec();

    for (unsigned i = 0; i < Dim; i++) {
        this->domain_m[i] = ippl::Index(nr_m[i]);
        this->decomp_m[i] = domainDecomposition[i];
    }

    this->setBCHandler(std::make_shared<BCHandler_t>(
        OPALFieldSolver_m->constructBCHandler()
    ));

    /// \todo so far, we only use true for all periodic and false for all open.
    bool isAllPeriodic = this->getBCHandler()->isAll(BCHandler_t::PERIODIC);
    m << level5 << "* FieldContainer set to isAllPeriodic = " << isAllPeriodic << endl;

    //      set stuff for pre_run i.e. warmup
    //      this will be reset when the correct computational
    //      domain is set

    Vector_t<double, Dim> length (6.0);
    this->hr_m = length / this->nr_m;
    this->origin_m = -3.0;
    this->dt_m = 0.5 / this->nr_m[2];

    rmin_m = origin_m;
    rmax_m = origin_m + length;

    this->setFieldContainer(std::make_shared<FieldContainer_t>(
        hr_m, rmin_m, rmax_m, decomp_m, domain_m, origin_m, isAllPeriodic
    ));

    this->setParticleContainer(std::make_shared<ParticleContainer_t>(
        this->fcontainer_m->getMesh(), this->fcontainer_m->getFL()));
    for (size_t i = 1; i < num_containers; ++i) {
        this->addParticleContainer(std::make_shared<ParticleContainer_t>(
            this->fcontainer_m->getMesh(), this->fcontainer_m->getFL()));
    }

    setSolver();

    // Build temporary accumulation fields after solver/field initialization so they
    // match the current mesh/layout and backing storage configuration.
    this->setTempEField(std::make_shared<VField_t<T, Dim>>());
    this->getTempEField()->initialize(this->fcontainer_m->getMesh(), this->fcontainer_m->getFL());
    this->setTempBField(std::make_shared<VField_t<T, Dim>>());
    this->getTempBField()->initialize(this->fcontainer_m->getMesh(), this->fcontainer_m->getFL());
    // -----------------------------------------------

    pre_run();
    
    globalPartPerNode_m = std::make_unique<size_t[]>(ippl::Comm->size());

    m << level5 << "* PartBunch constructor done." << endl;
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::do_binaryRepart() {
    using FieldContainer_t = FieldContainer<T, Dim>;
    std::shared_ptr<FieldContainer_t> fc = this->fcontainer_m;

    size_type totalP = this->getTotalNum();

    if (this->loadbalancer_m->balance(totalP)) {
        auto* mesh = &fc->getRho().get_mesh();
        auto* FL = &fc->getFL();
        this->loadbalancer_m->repartition(FL, mesh, isFirstRepartition_m);
    }
}

template <typename T, unsigned Dim>
void  PartBunch<T, Dim>::gatherLoadBalanceStatistics() {
        std::fill_n(globalPartPerNode_m.get(), ippl::Comm->size(), 0);  // Fill the array with zeros
        globalPartPerNode_m[ippl::Comm->rank()] = getLocalNum();
        ippl::Comm->allreduce(globalPartPerNode_m.get(), 
                              ippl::Comm->size(), 
                              std::plus<size_t>());
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::setSolver() {
    Inform m("PartBunch::setSolver");
    m << level2 << "Initializing solver: " << OPALFieldSolver_m->getType() << endl;
    if (this->solver_m != "")
        m << level1 << "Warning solver already initiated but overwrite ..." << endl;

    this->solver_m = OPALFieldSolver_m->getType();

    this->fcontainer_m->initializeFields(this->solver_m);

    // Needs to happen before setting the field solver, since the field solver needs the bins.
    setBins();

    this->setFieldSolver(
        std::make_shared<BinnedFieldSolver<T, Dim>>(
            this->solver_m, &this->fcontainer_m->getRho(), &this->fcontainer_m->getE(),
            &this->fcontainer_m->getPhi(), this->getBCHandler(),
            hasBinning() ? OPALFieldSolver_m->getBinningCmd()->getTablePrintFrequency() : 0
        )
    );
    m << level4 << "Binned field solver set (binned or legacy at runtime)." << endl;

    this->fsolver_m->initSolver();
    m << level4 << "Field solver initialized." << endl;

    /// ADA we need to be able to set a load balancer when not having a field solver
    this->setLoadBalancer(std::make_shared<LoadBalancer_t>(
        this->lbt_m, 
        this->fcontainer_m, 
        this->pcontainer_m, 
        this->fsolver_m
    ));
    m << level3 << "Solver and Load Balancer set." << endl;
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::setBins() {
    Inform m("PartBunch::setBins");

    BinningCmd* binningCmd = OPALFieldSolver_m->getBinningCmd();

    if (!OPALFieldSolver_m->hasBinningCmd()) {
        m << level2 << "Solver " << OPALFieldSolver_m->getOpalName() << " has no binning command attached, not using binning." << endl;
        return;
    }

    m << level4 << "Using binning command: " << binningCmd->getOpalName() << endl;

    switch (binningCmd->getParameterType()) {
        case BinningParameter::VELOCITYZ:
            this->setBins(
                std::make_shared<ParticleBinning::AdaptBins<ParticleContainer_t, CoordinateSelector_t>>(
                    this->getParticleContainer(), 
                    CoordinateSelector_t(2), 
                    binningCmd->getMaxBins(),
                    binningCmd->getBinningAlpha(), 
                    binningCmd->getBinningBeta(),
                    binningCmd->getDesiredWidth(), 
                    binningCmd->getOpalName()
                )
            );
            break;
        case BinningParameter::GAMMAZ:
            this->setBins(
                std::make_shared<ParticleBinning::AdaptBins<ParticleContainer_t, GammaSelector_t>>(
                    this->getParticleContainer(),
                    GammaSelector_t(2), 
                    binningCmd->getMaxBins(),
                    binningCmd->getBinningAlpha(), 
                    binningCmd->getBinningBeta(),
                    binningCmd->getDesiredWidth(), 
                    binningCmd->getOpalName()
                )
            );
            break;
        default:
            throw OpalException(
                "PartBunch::setBins", "Binning parameter " + binningCmd->getParameter()
                                          + " not supported yet! Only VELOCITYZ and GAMMAZ.");
    }
    m << level3 << "Bins set." << endl;
    this->getBins()->debug();
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::calcBeamParameters() {
    Inform m("PartBunch::calcBeamParameters");
    std::shared_ptr<ParticleContainer_t> pc = this->pcontainer_m;
    
    using view_type = ippl::ParticleAttrib<Vector_t<double,3>>::view_type;
    view_type Rview = pc->R.getView();
    view_type Pview = pc->P.getView();
    this->updateMoments();
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
        const size_t nLocal = this->getLocalNum();
        Kokkos::parallel_reduce("calc moments of particle distr.", nLocal,
            KOKKOS_LAMBDA(const size_t k, double& cent, double& mom0, double& mom1, 
                          double& mom2, double& mom3, double& mom4, 
                          double& mom5) {
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

    /// \todo do this in one step much nicer with ippl::Vector...
    for (unsigned d = 0; d < Dim; ++d) {
        Kokkos::parallel_reduce("rel max", this->getLocalNum(),
            KOKKOS_LAMBDA(const int i, double& mm) {
                double tmp_vel = Rview(i)[d];
                mm             = tmp_vel > mm ? tmp_vel : mm;
            }, Kokkos::Max<T>(rmax_loc[d]));
        
        Kokkos::parallel_reduce("rel min", this->getLocalNum(),
            KOKKOS_LAMBDA(const int i, double& mm) {
                double tmp_vel = Rview(i)[d];
                mm             = tmp_vel < mm ? tmp_vel : mm;
            }, Kokkos::Min<T>(rmin_loc[d]));
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

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::pre_run() {
    Inform m("PartBunch::pre_run");
    m << level2 << "Starting pre_run..." << endl;
    this->fcontainer_m->getRho() = 0.0;
    m << level4 << "Rho initialized to zero." << endl;

    /*
    Force skip field dump during pre_run/warmup!
    In order to call runSolver without field dump, we need to dynamic cast 
    fsolver_m to FieldSolver_t, since this addition is not possible in the base
    class (without changing ippl).
    */
    this->getFieldSolver()->runSolver(true);
    m << level4 << "Field solver ran during pre_run." << endl;
    this->getFieldSolver()->resetCallCounter();
    m << level4 << "Call counter reset. pre_run done." << endl;
}

template <typename T, unsigned Dim>
double PartBunch<T, Dim>::get_meanKineticEnergy() {
    // Single source of truth: computed in DistributionMoments during updateMoments().
    // Unit: MeV (see DistributionMoments implementation).
    return this->pcontainer_m->getMeanKineticEnergy();
}

template <typename T, unsigned Dim>
Inform& PartBunch<T, Dim>::print(Inform& os) {
    // if (this->getLocalNum() != 0) {  // to suppress Nans
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
        const auto qmMode             = pc->getQMStorageMode();
        if (qmMode == ParticleContainer_t::QMStorageMode::Attributes) {
            qmStorageModeStr = "ATTRIBUTES";
        }

        os << level1 << std::scientific << "\n"
           << "* ************** B U N C H "
            "********************************************************* \n"
           << "* CONTAINER       = " << ci << "\n"
           << "* PARTICLES       = " << pc->getTotalNum() << "\n"
           << "* CHARGE          = " << this->getCharge(ci) << " (Cb) \n"
           << "* QM STORAGE MODE = " << qmStorageModeStr << "\n"
           << "* <EKIN>          = " << Util::getEnergyString(ek) << "\n"
           << "* <dEKIN>         = " << Util::getEnergyString(dek) << "\n"
           << "* INTEGRATOR      = " << integration_method_m << "\n"
           << "* MIN R (origin)  = " << Util::getLengthString(pc->getMinR(), 5) << "\n"
           << "* MAX R (max ext) = " << Util::getLengthString(pc->getMaxR(), 5) << "\n"
           << "* RMS R           = " << Util::getLengthString(pc->getRmsR(), 5) << "\n"
           << "* RMS P           = " << pc->getRmsP() << " [beta gamma]\n"
           << "* Mean R          = " << pc->getMeanR() << " [m]\n"
           << "* Mean P          = " << pc->getMeanP() << " [beta gamma]\n"
           << "* MESH SPACING    = " << Util::getLengthString(this->fcontainer_m->getMesh().getMeshSpacing(), 5) << "\n"
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
        os << level1 << "* "
           "********************************************************************************"
           "** \n"
           << endl;
    }

    os.flags(ff);
    return os;
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::bunchUpdate() {
    Inform m ("PartBunch::bunchUpdate");
    m << level4 << "Updating bunch and doing repartitioning if needed." << endl;
    /* \brief
       1. calculates and set hr
       2. do repartitioning
    */

    auto *mesh = &this->fcontainer_m->getMesh();
    auto *FL   = &this->fcontainer_m->getFL();

    std::shared_ptr<ParticleContainer_t> pc = this->getParticleContainer();

    pc->computeMinMaxR();

    ippl::Vector<double, 3> o = pc->getMinR();
    ippl::Vector<double, 3> e = pc->getMaxR();
    ippl::Vector<double, 3> l = e - o;

    /*
    If a coordinate of l is too close to zero, set it to 1e-12.
    This avoids having a mesh spacing of zero, which would crash ippl and allows
    empty simulations - especially important for emission sources.
    */
    for (int i = 0; i < 3; i++) {
        if (l[i] < 1e-6) { 
            l[i] = 1e-6; 
            m << level3 << "Mesh spacing in dimension " << i << " too small. Set to 1e-6." << endl;
            //return;
        }
    }

    /*
    Now matches OPAL: domain + incr% on each side.
    Note that there is still a mismatch: OPAL only resizes in z direction and
    keeps x/y the same. But this doesn't make too much sense in my opinion...
    */
    // Update origin and extent for the FieldContainer (not for the particles!)
    o = o - l*this->OPALFieldSolver_m->getBoxIncr()/100.;
    e = e + l*this->OPALFieldSolver_m->getBoxIncr()/100.;
    l = e - o;

    hr_m = l / this->nr_m;

    mesh->setMeshSpacing(hr_m);
    mesh->setOrigin(o);

    /*
    I think these in the field container should reflect mesh boundaries, not 
    particle boundaries, since the field solver needs to know the mesh and solve
    */
    this->getFieldContainer()->setRMin(o);
    this->getFieldContainer()->setRMax(e);
    this->getFieldContainer()->setHr(hr_m);

    m << level3 << "Field Container updated with new mesh boundaries and spacing:" << endl;
    m << level3 << "\t\t> Mesh origin:   " << mesh->getOrigin() << endl;
    m << level3 << "\t\t> Mesh spacing:  " << hr_m << endl;
    m << level3 << "\t\t> Box increment: " << this->OPALFieldSolver_m->getBoxIncr() << "%" << endl;

    pc->getLayout().updateLayout(*FL, *mesh);
    pc->update();
    m << level5 << "Particle container updated with new layout." << endl;

    this->isFirstRepartition_m = true;
    //this->loadbalancer_m->initializeORB(FL, mesh);
    //this->loadbalancer_m->repartition(FL, mesh, this->isFirstRepartition_m);
    m << level5 << "Load balancer repartitioning done." << endl;

    this->updateMoments();
    m << level5 << "Moments updated." << endl;
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::computeSelfFields() {
    using BinnedSolver_t = BinnedFieldSolver<T, Dim>;

    BinnedSolver_t* bsolver = dynamic_cast<BinnedSolver_t*>(this->fsolver_m.get());
    if (!bsolver) {
        throw OpalException("PartBunch::computeSelfFields",
                            "Field solver is not a BinnedFieldSolver instance.");
    }

    // Non-owning shared_ptr alias: solver does not manage bunch lifetime.
    std::shared_ptr<PartBunch<T, Dim>> bunchPtr(this, [](PartBunch<T, Dim>*) { });
    bsolver->computeSelfFields(bunchPtr);
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::dumpBinConfig(bool preMerge) {
    if (!hasBinning() || !dataSink_m) {
        throw OpalException("PartBunch::dumpBinConfig", 
                            "No binning or data sink set, but dumpBinConfig() was called.");
    }

    Inform m("PartBunch::dumpBinConfig");

    BinningCmd* binningCmd = OPALFieldSolver_m->getBinningCmd();
    if (!binningCmd) {
        return;
    }

    const long long step = getGlobalTrackStep();
    const int dumpFreq   = binningCmd->getDumpBinsFrequency();
    if (dumpFreq <= 0 || (step % dumpFreq) != 0) {
        return;
    }

    std::shared_ptr<AdaptBins_t> bins = getBins();
    if (!bins) {
        return;
    }

    std::vector<typename AdaptBins_t::size_type> countsHost;
    std::vector<typename AdaptBins_t::value_type> widthsHost;
    const auto xMin = bins->getBinConfigHost(countsHost, widthsHost);

    std::vector<std::size_t> counts(countsHost.begin(), countsHost.end());
    std::vector<double> widths(widthsHost.begin(), widthsHost.end());

    m << level5 << "Dumping bin configuration (preMerge=" << (preMerge ? 1 : 0)
      << ") at globalTrackStep=" << step
      << " with nBins=" << counts.size()
      << " to file \"" << binningCmd->getDumpBinsFileName() << "\"." << endl;

    dataSink_m->dumpBinConfig(
        step,
        getT(),
        preMerge,
        counts,
        widths,
        static_cast<double>(xMin),
        binningCmd->getDumpBinsFileName());
}

template <typename T, unsigned Dim>
void PartBunch<T,Dim>::performBunchSanityChecks() const {
    Inform ms("PartBunch::performBunchSanityChecks");
    ms << level4 << "========== Performing sanity checks on PartBunch... ==========" << endl;
    /// \todo always try to add more checks here! Best practice: throw explanatory exceptions and give output when passed.

    // Check if bc handler was initialized properly
    if (!this->getBCHandler()) {
        throw OpalException("PartBunch::performBunchSanityChecks", 
                            "BC Handler not initialized properly.");
    }
    ms << level4 << "BC Handler initialized properly." << endl;

    if (!hasFieldSolver()) {
        throw OpalException("PartBunch::performBunchSanityChecks", 
                            "Field Solver was not initialized.");
    }
    ms << level4 << "Field Solver object was initialized." << endl;

    // Verify we can access the concrete FieldSolver and its internals
    auto fs = std::dynamic_pointer_cast<FieldSolver_t>(this->fsolver_m);

    // cannot use getFieldContainer, since this getter cannot be const!
    const std::shared_ptr<FieldContainer<T, Dim>> fctr = this->fcontainer_m;
    if (!fctr) {
        throw OpalException("PartBunch::performBunchSanityChecks",
                            "FieldContainer isn't initialized correctly.");
    }

    // Check internal field pointers are set
    if (fs->getRho() == nullptr || 
        fs->getE()   == nullptr || 
        fs->getPhi() == nullptr) {
        throw OpalException("PartBunch::performBunchSanityChecks",
                            "FieldSolver internal fields (rho/E/phi) not assigned.");
    }
    ms << level4 << "FieldSolver internal field pointers are set." << endl;

    // Ensure FieldSolver fields point to our FieldContainer's fields
    if (fs->getRho() != &fctr->getRho() ||
        fs->getE()   != &fctr->getE()   ||
        fs->getPhi() != &fctr->getPhi()) {
        throw OpalException("PartBunch::performBunchSanityChecks",
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
        throw OpalException("PartBunch::performBunchSanityChecks",
                            "FieldSolver type string is empty.");
    }
    if (stype != "FFT" && 
        stype != "OPEN" && 
        stype != "CG" && 
        stype != "NONE") {
        throw OpalException("PartBunch::performBunchSanityChecks",
                            "Unsupported FieldSolver type: " + stype);
    }
    ms << level4 << "FieldSolver type: " << stype << endl;

    // Basic check that the E-field layout has non-zero extent
    auto Eview = fctr->getE().getView();
    if (Eview.extent(0) == 0 || Eview.extent(1) == 0 || Eview.extent(2) == 0) {
        throw OpalException("PartBunch::performBunchSanityChecks",
                            "E-field layout not initialized (zero extent). ");
    }
    ms << level4 << "E-field layout initialized." << endl;

    ms << level2 << "========= Done performing PartBunch sanity checks... =========" << endl;
}


// Explicit instantiations
template class PartBunch<double, 3>;
