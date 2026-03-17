#include "PartBunch/PartBunch.h"
#include "Algorithms/Matrix.h"
#include "Utilities/Util.h"
#include "Structure/DataSink.h"

#undef doDEBUG

template <typename T, unsigned Dim>
PartBunch<T, Dim>::PartBunch(double qi,
                             double mi,
                             size_t totalP,
                             /*int nt,*/
                             double lbt,
                             std::string integration_method,
                             std::shared_ptr<FieldSolverCmd>& OPALFieldSolver,
                             std::shared_ptr<DataSink> dataSink)
    : ippl::PicManager<T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>, LoadBalancer<T, Dim>>(),
      time_m(0.0),
      totalP_m(totalP),
      //nt_m(nt),
      lbt_m(lbt),
      dt_m(0),
      it_m(0),
      integration_method_m(integration_method),
      solver_m(""),
      isFirstRepartition_m(true),
      qi_m(qi),
      mi_m(mi),
      rmsDensity_m(0.0),
      RefPartR_m(0.0),
      RefPartP_m(0.0),
      localTrackStep_m(0),
      globalTrackStep_m(0),
      OPALFieldSolver_m(OPALFieldSolver),
      dataSink_m(std::move(dataSink)) {

    Inform m("PartBunch::PartBunch");
    m << level4 << "PartBunch Constructor" << endl;

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
        this->fcontainer_m->getMesh(), this->fcontainer_m->getFL()
    ));

    this->setTempEField(std::make_shared<VField_t<T, Dim>>(
        this->fcontainer_m->getE()
    )); // user copy constructor
    this->getTempEField()->initialize(this->fcontainer_m->getMesh(), 
                                      this->fcontainer_m->getFL());
    // -----------------------------------------------

    setSolver();

    pre_run();
    
    globalPartPerNode_m = std::make_unique<size_t[]>(ippl::Comm->size());

    m << level5 << "* PartBunch constructor done." << endl;
}

template <typename T, unsigned Dim>
T PartBunch<T, Dim>::getCouplingConstant() const {
    /*
    This function needs to be here, since FieldSoler_t is only fully defined
    at instanciation of PartBunch, so not yet in the header file.
    */
    
    if (!hasFieldSolver()) {
        throw OpalException("PartBunch::getCouplingConstant",
                            "Cannot return coupling if fsolver_m is not a "
                            "FieldSolver instance");
    }
    return this->getFieldSolver()->getCouplingConstant();
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::gatherCIC() {
    using Base = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;
    typename Base::particle_position_type* Ep = &this->pcontainer_m->E;
    typename Base::particle_position_type* R = &this->pcontainer_m->R;
    VField_t<T, Dim>* Ef = &this->fcontainer_m->getE();
    gather(*Ep, *Ef, *R);
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
    
    this->setFieldSolver(std::make_shared<FieldSolver_t>(
        this->solver_m, 
        &this->fcontainer_m->getRho(), 
        &this->fcontainer_m->getE(),
        &this->fcontainer_m->getPhi(),
        this->getBCHandler()
    ));
    m << level4 << "Field solver set." << endl;

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

    setBins();
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

    std::string parameterName = binningCmd->getParameter();
    if (parameterName != "VELOCITYZ") {
        throw OpalException("PartBunch::setBins",
                            "Binning parameter " + parameterName + " not supported yet! Only VELOCITYZ.");
    }

    this->setBins(std::make_shared<AdaptBins_t>(
        this->getParticleContainer(),
        BinningSelector_t(2),
        binningCmd->getMaxBins(),
        binningCmd->getBinningAlpha(),
        binningCmd->getBinningBeta(),
        binningCmd->getDesiredWidth() // Cost function parameters
    ));
    m << level3 << "Bins set." << endl;
    this->getBins()->debug();
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::spaceChargeEFieldCheck(Vector_t<double, 3> /*efScale*/) {
    Inform msg("PartBunch::spaceChargeEFieldCheck");

    auto pE_view   = this->pcontainer_m->E.getView();
    auto fphi_view = this->fcontainer_m->getPhi().getView();


    
    double avgphi        = 0.0;
    double avgE          = 0.0;
    double minEComponent = std::numeric_limits<T>::max();
    double maxEComponent = std::numeric_limits<T>::min();
    double minE          = std::numeric_limits<T>::max();
    double maxE          = std::numeric_limits<T>::min();
    double cc            = getCouplingConstant();
    
    int myRank = ippl::Comm->rank();

    Kokkos::parallel_reduce(
        "check e-field", this->getLocalNum(),
        KOKKOS_LAMBDA(const int i, double& loc_avgE, double& loc_minEComponent,
                      double& loc_maxEComponent, double& loc_minE, 
                      double& loc_maxE) {
            double EX    = pE_view[i][0]*cc;
            double EY    = pE_view[i][1]*cc;
            double EZ    = pE_view[i][2]*cc;
        
            double ENorm = Kokkos::sqrt(EX*EX + EY*EY + EZ*EZ);
            
            loc_avgE += ENorm;

            loc_minEComponent = EX < loc_minEComponent ? EX : loc_minEComponent;
            loc_minEComponent = EY < loc_minEComponent ? EY : loc_minEComponent;
            loc_minEComponent = EZ < loc_minEComponent ? EZ : loc_minEComponent;
            
            loc_maxEComponent = EX > loc_maxEComponent ? EX : loc_maxEComponent;
            loc_maxEComponent = EY > loc_maxEComponent ? EY : loc_maxEComponent;
            loc_maxEComponent = EZ > loc_maxEComponent ? EZ : loc_maxEComponent;   

            loc_minE = ENorm < loc_minE ? ENorm : loc_minE;
            loc_maxE = ENorm > loc_maxE ? ENorm : loc_maxE;
        },
        Kokkos::Sum<T>(avgE), Kokkos::Min<T>(minEComponent),
        Kokkos::Max<T>(maxEComponent), Kokkos::Min<T>(minE),
        Kokkos::Max<T>(maxE)
    );

    if (this->getLocalNum() == 0) {
        minEComponent = maxEComponent = minE = maxE = avgE = 0.0;
    }

    MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &avgE, &avgE, 1, MPI_DOUBLE, MPI_SUM, 0, ippl::Comm->getCommunicator());
    MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &minEComponent, &minEComponent, 1, MPI_DOUBLE, MPI_MIN, 0, ippl::Comm->getCommunicator());
    MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &maxEComponent, &maxEComponent, 1, MPI_DOUBLE, MPI_MAX, 0, ippl::Comm->getCommunicator());
    MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &minE, &minE, 1, MPI_DOUBLE, MPI_MIN, 0, ippl::Comm->getCommunicator());
    MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &maxE, &maxE, 1, MPI_DOUBLE, MPI_MAX, 0, ippl::Comm->getCommunicator());
 
    size_t Np = this->getTotalNum();
    avgE /= (Np == 0) ? 1 : Np; // avoid division by zero for empty simulations (see also DistributionMoments::computeMeans implementation) 

    msg << level4 << "avgENorm = " << avgE << endl;
    
    using mdrange_type = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

    Kokkos::parallel_reduce(
        "check phi", mdrange_type({0,0,0}, {fphi_view.extent(0),fphi_view.extent(1),fphi_view.extent(2)}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, double& loc_avgphi) {
            double phi = fphi_view(i,j,k);
            loc_avgphi += phi;
        },
        Kokkos::Sum<T>(avgphi)
    );

    MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &avgphi, &avgphi, 1, MPI_DOUBLE, MPI_SUM, 0, ippl::Comm->getCommunicator());
    avgphi /= this->getTotalNum(); 
    msg << level4 << "avgphi = " << avgphi << endl;

}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::calcBeamParameters() {
    Inform m("PartBunch::calcBeamParameters");
    std::shared_ptr<ParticleContainer_t> pc = this->pcontainer_m;
    
    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();
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
        Kokkos::parallel_reduce(
            "calc moments of particle distr.", ippl::getRangePolicy(Rview),
            KOKKOS_LAMBDA(const int k, double& cent, double& mom0, double& mom1, 
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
double PartBunch<T, Dim>::getdE() const {
    // Single source of truth: computed in DistributionMoments during updateMoments().
    // Unit: MeV (see DistributionMoments implementation).
    return this->pcontainer_m->getStdKineticEnergy();
}

template <typename T, unsigned Dim>
Inform& PartBunch<T, Dim>::print(Inform& os) {
    // if (this->getLocalNum() != 0) {  // to suppress Nans
    Inform::FmtFlags_t ff = os.flags();

    const double ek  = this->get_meanKineticEnergy();
    const double dek = this->getdE();
    
    os << level1 << std::scientific << "\n"
       << "* ************** B U N C H "
        "********************************************************* \n"
       << "* PARTICLES       = " << this->getTotalNum() << "\n"
       << "* CHARGE          = " << this->qi_m*this->getTotalNum() << " (Cb) \n"
       << "* <EKIN>          = " << ek << " (MeV) \n"
       << "* <dEKIN>         = " << dek << " (MeV) \n"
       << "* INTEGRATOR      = " << integration_method_m << "\n"
       << "* MIN R (origin)  = " << Util::getLengthString( this->pcontainer_m->getMinR(), 5) << "\n"
       << "* MAX R (max ext) = " << Util::getLengthString( this->pcontainer_m->getMaxR(), 5) << "\n"
       << "* RMS R           = " << Util::getLengthString( this->pcontainer_m->getRmsR(), 5) << "\n"
       << "* RMS P           = " << this->pcontainer_m->getRmsP() << " [beta gamma]\n"
       << "* Mean R          = " << this->pcontainer_m->getMeanR() << " [m]\n"
       << "* Mean P          = " << this->pcontainer_m->getMeanP() << " [beta gamma]\n"
       << "* MESH SPACING    = " << Util::getLengthString( this->fcontainer_m->getMesh().getMeshSpacing(), 5) << "\n"
       << "* COMPDOM INCR    = " << this->OPALFieldSolver_m->getBoxIncr() << " (%) \n"
       << "* FIELD LAYOUT    = " << this->fcontainer_m->getFL() << "\n"
       << "* Centroid : \n* ";
    for (unsigned int i=0; i<2*Dim; i++) {
        os << level1 << this->pcontainer_m->getCentroid()[i] << " ";
    }
    os << level1 << endl << "* Cov Matrix : \n* ";
    for (unsigned int i=0; i<2*Dim; i++) {
        for (unsigned int j=0; j<2*Dim; j++) {
            os << level1 << this->pcontainer_m->getCovMatrix()(i,j) << " ";
        }
        os << level1 << "\n* ";
    }
    os << level1 << "* "
        "********************************************************************************"
        "** \n"
       << endl;
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
    Inform m("PartBunch::computeSelfFields");

    if (this->hasBinning()) {
        static IpplTimings::TimerRef completeBinningT = IpplTimings::getTimer("bTotalBinningT");

        // Start binning and sorting by bins
        std::shared_ptr<AdaptBins_t> bins = this->getBins();

        IpplTimings::startTimer(completeBinningT);
        bins->doFullRebin(bins->getMaxBinCount());
        dumpBinConfig(true);   // pre-merge configuration
        bins->sortContainerByBin();    // Sort BEFORE merging to reduce atomics overhead.
        bins->genAdaptiveHistogram();  // merge bins adaptively
        dumpBinConfig(false);  // post-merge configuration

        IpplTimings::stopTimer(completeBinningT);
        m << level4 << "Binning routine done." << endl;
    } else {
        m << level4 << "No AdaptBins object present, not using binning." << endl;
    }

    /*
    I would guess that ths bunchUpdate is only necessary after a push (where we
    need it anyways, since positions have changed). However, when removing it,
    the total energy of the FODO example quickly diverges to "-inf". I don't
    know why this is, but particle positions shouldn't have changed. I would
    therefore assume that we could separate bunchUpdate from pc->update() in 
    order to save some computation. 
    */
    this->bunchUpdate();
    m << level5 << "Bunch updated." << endl;

    /// \todo Add binned field solver here (needs iteration over bins, scatterPerBin calls and Etmp build up)! See https://gitlab.psi.ch/OPAL/opal-x/src/-/blame/binnedFieldSolver/src/PartBunch/PartBunch.cpp?ref_type=heads#L376

    ippl::ParticleAttrib<T>* Q               = &this->pcontainer_m->Q;
    typename Base::particle_position_type* R = &this->pcontainer_m->R;

    this->fcontainer_m->getRho()             = 0.0;
    Field_t<Dim>* rho                        = &this->fcontainer_m->getRho();

    /// \todo replace with scatterCIC? --> later with scatterPerBin!
    // Charge "unit" here is "charge per macroparticle" [C]!
    *Q = (*Q) * this->pcontainer_m->dt; // Scale by time step
    scatter(*Q, *rho, *R); 
    *Q = (*Q) / this->pcontainer_m->dt; // Rescale back to charge per macroparticle
    m << level4 << "Scatter done." << endl;

    /*
    Now rho is in units of [C * s] -- need to divide by dt to get back to [C].
    Note By using the global timestep getdT(), we account for the possibility of 
    "fractional timesteps", meaning the particle was virtually "created in the
    middle" of a full timestep. As of now, this might not be necessary.
    */
    (*rho) = (*rho) / getdT(); 
    m << level4 << "Rho scale by dt done." << endl;

#ifdef doDEBUG
    const double qtot                        = this->qi_m * this->getTotalNum();
    size_type TotalParticles                 = 0;
    size_type localParticles                 = this->pcontainer_m->getLocalNum();
   
    double relError                          = std::fabs((qtot - (*rho).sum()) / qtot);
    
    ippl::Comm->reduce(localParticles, TotalParticles, 1, std::plus<size_type>());
    
    if ((ippl::Comm->rank() == 0) && (relError > 1.0E-13)) {
            Inform m2("PartBunch::computeSelfFields2", INFORM_ALL_NODES);
            m2 << "Time step: " << it_m
               << " total particles in the sim. " << totalP_m 
               << " missing : " << totalP_m-TotalParticles 
               << " rel. error in charge conservation: " << relError << endl;
    }
#endif

    // At this point, the units of rho need to be corrected: rho = rho / cellVolume
    if (this->fsolver_m->getStype() != "FEM" && this->fsolver_m->getStype() != "FEM_PRECON") {
        // CG solver already accounts for cell volume internally in FD,
        // other solvers need explicit normalization here
        double cellVolume = std::reduce(hr_m.begin(), hr_m.end(), 1., std::multiplies<double>());
        (*rho)            = (*rho) / cellVolume;
        m << level4 << "Rho normalized by cell volume: " << cellVolume << "." << endl;
    }

    // Alpine uses net 0 charge density for periodic BCs, so we need to subtract background charge here (?TODO: check)
    // Note: otherwise solvers like the CG solver will "explode" and give unnormalized potentials?
    double totalQ = getCharge();
    if (this->fsolver_m->getStype() != "OPEN") {
        double size = 1;
        for (size_t d = 0; d < 3; d++) {
            size *= rmax_m[d] - rmin_m[d];
        }

        (*rho) = (*rho) - (totalQ / size);
        m << level4 << "Net-0 charge generation with factor " << (totalQ / size) << " done." << endl;
    }

    /*
    This concludes the scatter step ( \todo perhaps we want to put this in its own function, like scatterCIC)

    Next is the solver step. The solver computes the E-field from rho directly. However,
    at the moment, the potential has units of [C/m^3].

    From OPAL:
    The scalar potential is given back with rho_m in units [C/m] = [F*V/m] and must be divided by
    4*pi*\epsilon_0 [F/m] resulting in [V].

    @note rho is overloaded and becomse phi when runSolver is called!
    */
    (*rho) = (*rho) * this->getCouplingConstant(); // now rho_m has units of [V]
    m << level5 << "Rho coupling applied." << endl;

    this->fsolver_m->runSolver();
    m << level4 << "Field solver ran." << endl;

    /*
    Now, with E=-grad(phi), E has units of [V/m] (note, phi is a scalar potential).

    Note that we do this step separately, even though the solvers can provide
    it firectlt, since some field solver then don't produce phi, which might
    be necessary elsewhere or for debugging.
    
    This can be optimized later by changing the output type of the solvers and 
    removing this line! As it is implemented now, it will always provide both 
    (phi and the E field that is actually used for kicking the particles).
    */
    
    gather(this->pcontainer_m->E, this->fcontainer_m->getE(), this->pcontainer_m->R);
    m << level4 << "Gather done." << endl;

    /// \todo put back in
    /*Vector_t<double, 3> efScale = Vector_t<double,3>(
        gammaz*cc/hr_scaled[0], 
        gammaz*cc/hr_scaled[1], 
        cc / gammaz / hr_scaled[2]
    );
    m << "E-field scale = " << efScale << endl;    
    spaceChargeEFieldCheck(efScale);*/
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
void PartBunch<T,Dim>::scatterCICPerBin(PartBunch<T,Dim>::binIndex_t binIndex) {
    /**
     * Scatters only particles in bin binIndex. Scatters all particles if binIndex=-1
     */

    throw OpalException("PartBunch::scatterCICPerBin", 
        "This function is not implemented yet! Please use scatterCIC for now.");

    Inform m("PartBunch::scatterCICPerBin");
    m << "Scattering binIndex = " << binIndex << " to grid." << endl;

    this->fcontainer_m->getRho() = 0.0;

    ippl::ParticleAttrib<T>* q               = &this->pcontainer_m->Q;
    typename Base::particle_position_type* R = &this->pcontainer_m->R;
    Field_t<Dim>* rho                        = &this->fcontainer_m->getRho();
    
    double Q;
    Vector_t<double, 3> rmin                 = rmin_m;
    Vector_t<double, 3> rmax                 = rmax_m;
    Vector_t<double, 3> hr                   = hr_m;

    if (binIndex == -1) {
        // Use original scatterCIC logic for all particles
        Q = this->qi_m * this->getTotalNum();
        scatter(*q, *rho, *R);
    } else {
        // Use per-bin scattering logic
        Q = this->qi_m * this->bins_m->getNPartInBin(binIndex, true);
        scatter(*q, *rho, *R, this->bins_m->getBinIterationPolicy(binIndex), this->bins_m->getHashArray());
    }

    m << "gammz= " << this->pcontainer_m->getMeanP()[2] << endl;
    
#ifdef doDEBUG
    double relError = std::fabs((Q - (*rho).sum()) / Q);
    size_type TotalParticles = 0;
    size_type localParticles = this->pcontainer_m->getLocalNum();
    size_type totalP_check = (binIndex == -1) ? totalP_m : this->pcontainer_m->getTotalNum();

    m << "computeSelfFields sum rho = " << (*rho).sum() << ", relError = " << relError << endl;
    
    ippl::Comm->reduce(localParticles, TotalParticles, 1, std::plus<size_type>());

    if (ippl::Comm->rank() == 0) {
        if (TotalParticles != totalP_check || relError > 1e-10) {
            m << "Time step: " << it_m << endl;
            m << "Total particles in the sim. " << totalP_check << " "
              << "after update: " << TotalParticles << endl;
            m << "Rel. error in charge conservation: " << relError << endl;
            ippl::Comm->abort();
        }
    }
#endif
    
    double cellVolume = std::reduce(hr.begin(), hr.end(), 1., std::multiplies<double>());
        
    m << "cellVolume= " << cellVolume << endl;

    (*rho) = (*rho) / cellVolume;

    // rho = rho_e - rho_i (only if periodic BCs)
    if (this->fsolver_m->getStype() != "OPEN") {
        double size = 1;
        for (unsigned d = 0; d < 3; d++) {
            size *= rmax[d] - rmin[d];
        }
        *rho = *rho - (Q / size);
    }
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
