#include "PartBunch/PartBunch.h"
#include <boost/numeric/ublas/io.hpp>
#include "Utilities/Util.h"

#undef doDEBUG

template <typename T, unsigned Dim>
PartBunch<T, Dim>::PartBunch(double qi, double mi, size_t totalP, int nt, double lbt,
                             std::string integration_method, std::shared_ptr<Distribution> &OPALdistribution,
                             std::shared_ptr<FieldSolverCmd> &OPALFieldSolver)
    : ippl::PicManager<T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>, LoadBalancer<T, Dim>>(),
      time_m(0.0),
      totalP_m(totalP),
      nt_m(nt),
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
      OPALdist_m(OPALdistribution),
      OPALFieldSolver_m(OPALFieldSolver) {

    static IpplTimings::TimerRef gatherInfoPartBunch = IpplTimings::getTimer("gatherInfoPartBunch");
    IpplTimings::startTimer(gatherInfoPartBunch);

    *gmsg << "PartBunch Constructor" << endl;

    //  get the needed information from OPAL FieldSolver command

    nr_m = Vector_t<int, Dim>(
            OPALFieldSolver_m->getNX(), OPALFieldSolver_m->getNY(), OPALFieldSolver_m->getNZ());

    const Vector_t<bool, 3> domainDecomposition = OPALFieldSolver_m->getDomDec();

    for (unsigned i = 0; i < Dim; i++) {
        this->domain_m[i] = ippl::Index(nr_m[i]);
        this->decomp_m[i] = domainDecomposition[i];
    }

    bool isAllPeriodic = true;  // \fixme need to get BCs from OPAL Fieldsolver

    //      set stuff for pre_run i.e. warmup
    //      this will be reset when the correct computational
    //      domain is set

    Vector_t<double, Dim> length (6.0);
    this->hr_m = length / this->nr_m;
    this->origin_m = -3.0;
    this->dt_m = 0.5 / this->nr_m[2];

    rmin_m = origin_m;
    rmax_m = origin_m + length;

    this->setFieldContainer( std::make_shared<FieldContainer_t>(hr_m, rmin_m, rmax_m, decomp_m, domain_m, origin_m, isAllPeriodic) );

    this->setParticleContainer(std::make_shared<ParticleContainer_t>(
        this->fcontainer_m->getMesh(), this->fcontainer_m->getFL()));

    IpplTimings::stopTimer(gatherInfoPartBunch);

    // ---------------- binning setup ----------------
    using bin_index_type = typename ParticleContainer_t::bin_index_type;
    bin_index_type maxBins = Options::maxBins;
    this->setBins(std::make_shared<AdaptBins_t>(
        this->getParticleContainer(), 
        BinningSelector_t(2), // TODO: hardcode z axis with coordinate selector at axis index 2
        static_cast<bin_index_type>(maxBins),
        Options::binningAlpha, Options::binningBeta, Options::desiredWidth // Cost function parameters
    ));
    this->getBins()->debug();

    this->setTempEField(std::make_shared<VField_t<T, Dim>>(this->fcontainer_m->getE())); // user copy constructor
    this->getTempEField()->initialize(this->fcontainer_m->getMesh(), this->fcontainer_m->getFL());
    // -----------------------------------------------

    static IpplTimings::TimerRef setSolverT = IpplTimings::getTimer("setSolver");
    IpplTimings::startTimer(setSolverT);
    setSolver(OPALFieldSolver_m->getType());
    IpplTimings::stopTimer(setSolverT);

    static IpplTimings::TimerRef prerun = IpplTimings::getTimer("prerun");
    IpplTimings::startTimer(prerun);
    pre_run();
    IpplTimings::stopTimer(prerun);

    globalPartPerNode_m = std::make_unique<size_t[]>(ippl::Comm->size());
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
    int it = this->it_m;

    if (this->loadbalancer_m->balance(totalP, it + 1)) {
        auto* mesh = &fc->getRho().get_mesh();
        auto* FL = &fc->getFL();
        this->loadbalancer_m->repartition(FL, mesh, isFirstRepartition_m);
    }
}

template <typename T, unsigned Dim>
void  PartBunch<T, Dim>::gatherLoadBalanceStatistics() {
        std::fill_n(globalPartPerNode_m.get(), ippl::Comm->size(), 0);  // Fill the array with zeros
        globalPartPerNode_m[ippl::Comm->rank()] = getLocalNum();
        ippl::Comm->allreduce(globalPartPerNode_m.get(), ippl::Comm->size(), std::plus<size_t>());
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::setSolver(std::string solver) {
    if (this->solver_m != "")
        *gmsg << "* Warning solver already initiated but overwrite ..." << endl;

    this->solver_m = solver;

    this->fcontainer_m->initializeFields(this->solver_m);
    
    this->setFieldSolver(std::make_shared<FieldSolver_t>(
                                                         this->solver_m, &this->fcontainer_m->getRho(), &this->fcontainer_m->getE(),
                                                         &this->fcontainer_m->getPhi()));

    this->fsolver_m->initSolver();
        
    /// ADA we need to be able to set a load balancer when not having a field solver
    this->setLoadBalancer(std::make_shared<LoadBalancer_t>(this->lbt_m, this->fcontainer_m, this->pcontainer_m, this->fsolver_m));
    
    *gmsg << "* Solver and Load Balancer set" << endl;
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::spaceChargeEFieldCheck(Vector_t<double, 3> efScale) {
    Inform msg("EParticleStats");

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
                                        double& loc_maxEComponent, double& loc_minE, double& loc_maxE) {
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
                            Kokkos::Max<T>(maxE));

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

 msg << "avgENorm = " << avgE << endl;
 
 using mdrange_type             = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

 Kokkos::parallel_reduce(
                         "check phi", mdrange_type({0,0,0}, {fphi_view.extent(0),fphi_view.extent(1),fphi_view.extent(2)}),
                         KOKKOS_LAMBDA(const int i, const int j, const int k, double& loc_avgphi) {
                             double phi = fphi_view(i,j,k);
                             loc_avgphi += phi;
                         },
                         Kokkos::Sum<T>(avgphi));

 MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &avgphi, &avgphi, 1, MPI_DOUBLE, MPI_SUM, 0, ippl::Comm->getCommunicator());
 avgphi /= this->getTotalNum(); 
 msg << "avgphi = " << avgphi << endl;

}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::calcBeamParameters() {
    std::shared_ptr<ParticleContainer_t> pc = this->pcontainer_m;
    
    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();

    ////////////////////////////////////
    //// Calculate Moments of R and P //
    ////////////////////////////////////


    double loc_centroid[2 * Dim]        = {};
    double loc_moment[2 * Dim][2 * Dim] = {};
        
    double centroid[2 * Dim]         = {};
    double moment[2 * Dim][2 * Dim]  = {};

    for (unsigned i = 0; i < 2 * Dim; i++) {
        loc_centroid[i] = 0.0;
        for (unsigned j = 0; j <= i; j++) {
            loc_moment[i][j] = 0.0;
            loc_moment[j][i] = 0.0;
        }
    }

    for (unsigned i = 0; i < 2 * Dim; ++i) {
        Kokkos::parallel_reduce(
                                "calc moments of particle distr.", ippl::getRangePolicy(Rview),
                                KOKKOS_LAMBDA(
                                              const int k, double& cent, double& mom0, double& mom1, double& mom2,
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

    MPI_Allreduce(loc_moment, moment, 2 * Dim * 2 * Dim, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
    MPI_Allreduce(loc_centroid, centroid, 2 * Dim, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    double rmax_loc[Dim];
    double rmin_loc[Dim];
    double rmax[Dim];
    double rmin[Dim];

    for (unsigned d = 0; d < Dim; ++d) {
        Kokkos::parallel_reduce(
                                "rel max", this->getLocalNum(),
                                KOKKOS_LAMBDA(const int i, double& mm) {
                                    double tmp_vel = Rview(i)[d];
                                    mm             = tmp_vel > mm ? tmp_vel : mm;
                                },
                                Kokkos::Max<T>(rmax_loc[d]));
        
        Kokkos::parallel_reduce(
                                "rel min", this->getLocalNum(),
                                KOKKOS_LAMBDA(const int i, double& mm) {
                                    double tmp_vel = Rview(i)[d];
                                    mm             = tmp_vel < mm ? tmp_vel : mm;
                                },
                                Kokkos::Min<T>(rmin_loc[d]));
    }
    Kokkos::fence();
    MPI_Allreduce(rmax_loc, rmax, Dim, MPI_DOUBLE, MPI_MAX, ippl::Comm->getCommunicator());
    MPI_Allreduce(rmin_loc, rmin, Dim, MPI_DOUBLE, MPI_MIN, ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    // \todo can we do this nicer? 
    for (unsigned int i=0; i<Dim; i++) {
        rmax_m(i) = rmax[i];
        rmin_m(i) = rmin[i];
    }
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::pre_run() {
    this->fcontainer_m->getRho() = 0.0;
    this->fsolver_m->runSolver();
}

template <typename T, unsigned Dim>
Inform& PartBunch<T, Dim>::print(Inform& os) {
    // if (this->getLocalNum() != 0) {  // to suppress Nans
    Inform::FmtFlags_t ff = os.flags();

    os << std::scientific;
    os << level1 << "\n";
    os << "* ************** B U N C H "
        "********************************************************* \n";
    os << "* PARTICLES       = " << this->getTotalNum() << "\n";
    os << "* CHARGE          = " << this->qi_m*this->getTotalNum() << " (Cb) \n";
    os << "* INTEGRATOR      = " << integration_method_m << "\n";
    os << "* MIN R (origin)  = " << Util::getLengthString( this->pcontainer_m->getMinR(), 5) << "\n";
    os << "* MAX R (max ext) = " << Util::getLengthString( this->pcontainer_m->getMaxR(), 5) << "\n";
    os << "* RMS R           = " << Util::getLengthString( this->pcontainer_m->getRmsR(), 5) << "\n";
    os << "* RMS P           = " << this->pcontainer_m->getRmsP() << " [beta gamma]\n";
    os << "* Mean R: " << this->pcontainer_m->getMeanR() << " [m]\n";
    os << "* Mean P: " << this->pcontainer_m->getMeanP() << " [beta gamma]\n";
    os << "* MESH SPACING    = " << Util::getLengthString( this->fcontainer_m->getMesh().getMeshSpacing(), 5) << "\n";
    os << "* COMPDOM INCR    = " << this->OPALFieldSolver_m->getBoxIncr() << " (%) \n";
    os << "* FIELD LAYOUT    = " << this->fcontainer_m->getFL() << "\n";
    os << "* Centroid : \n* ";
    for (unsigned int i=0; i<2*Dim; i++) {
        os << this->pcontainer_m->getCentroid()[i] << " ";
    }
    os << endl << "* Cov Matrix : \n* ";
    for (unsigned int i=0; i<2*Dim; i++) {
        for (unsigned int j=0; j<2*Dim; j++) {
            os << this->pcontainer_m->getCovMatrix()(i,j) << " ";
        }
        os << "\n* ";
    }
    os << "* "
        "********************************************************************************"
        "** "
       << endl;
    os.flags(ff);
    return os;
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::bunchUpdate(ippl::Vector<double, 3> hr) {
    /* \brief
       1. calculates and set hr
       2. do repartitioning
    */
    Inform m ("bunchUpdate ");
    
    auto *mesh = &this->fcontainer_m->getMesh();
    auto *FL   = &this->fcontainer_m->getFL();

    std::shared_ptr<ParticleContainer_t> pc = this->getParticleContainer();

    pc->computeMinMaxR();

    /// \brief assume o < 0.0?

    ippl::Vector<double, 3> o = pc->getMinR() - std::numeric_limits<T>::lowest();
    ippl::Vector<double, 3> e = pc->getMaxR() + std::numeric_limits<T>::lowest();
    ippl::Vector<double, 3> l = e - o;

    hr_m = (1.0+this->OPALFieldSolver_m->getBoxIncr()/100.)*(l / this->nr_m);
    mesh->setMeshSpacing(hr);
    mesh->setOrigin(o-0.5*hr_m*this->OPALFieldSolver_m->getBoxIncr()/100.);
    
    pc->getLayout().updateLayout(*FL, *mesh);
    pc->update();
    
    this->getFieldContainer()->setRMin(o);
    this->getFieldContainer()->setRMax(e);
    this->getFieldContainer()->setHr(hr);
    
    this->isFirstRepartition_m = true;
    this->loadbalancer_m->initializeORB(FL, mesh);
    this->loadbalancer_m->repartition(FL, mesh, this->isFirstRepartition_m);
    this->updateMoments();
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::bunchUpdate() {

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

    hr_m = (1.0+this->OPALFieldSolver_m->getBoxIncr()/100.)*(l / this->nr_m);
    mesh->setMeshSpacing(hr_m);
    mesh->setOrigin(o-0.5*hr_m*this->OPALFieldSolver_m->getBoxIncr()/100.);

    this->getFieldContainer()->setRMin(o);
    this->getFieldContainer()->setRMax(e);
    this->getFieldContainer()->setHr(hr_m);

    pc->getLayout().updateLayout(*FL, *mesh);
    pc->update();

    this->isFirstRepartition_m = true;
    this->loadbalancer_m->initializeORB(FL, mesh);
    //this->loadbalancer_m->repartition(FL, mesh, this->isFirstRepartition_m);

    this->updateMoments();
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::computeSelfFields() {
    static IpplTimings::TimerRef completeBinningT = IpplTimings::getTimer("bTotalBinningT");

    // Start binning and sorting by bins
    std::shared_ptr<AdaptBins_t> bins = this->getBins();
    VField_t<double, 3>& Etmp = *(this->getTempEField());

    IpplTimings::startTimer(completeBinningT);
    bins->doFullRebin(bins->getMaxBinCount()); // rebin with 128 bins // bins->getMaxBinCount()
    bins->print(); // For debugging...
    bins->sortContainerByBin(); // Sort BEFORE, since it generates less atomics overhead with more bins!

    
    bins->genAdaptiveHistogram(); // merge bins with width/N_part ratio of 1.0
    IpplTimings::stopTimer(completeBinningT);

    bins->print(); // For debugging...


    static IpplTimings::TimerRef SolveTimer = IpplTimings::getTimer("SolveTimer");
    IpplTimings::startTimer(SolveTimer);

    /*
      \todo check if Lorentz transform is needed

    double gammaz = this->pcontainer_m->getMeanGammaZ();
    gammaz *= gammaz;
    gammaz = std::sqrt(gammaz + 1.0);

    Vector_t<double, 3> hr_scaled = hr_m;
    //    hr_scaled[2] *= gammaz;

    hr_m = hr_scaled;    

    */

    /*
      particles have moved need to adjust grid
    */
    
    this->bunchUpdate();
    
    /*

     scatterCIC start

    */

    /// \todo Add binned field solver here (needs iteration over bins, scatterPerBin calls and Etmp build up)! See https://gitlab.psi.ch/OPAL/opal-x/src/-/blame/binnedFieldSolver/src/PartBunch/PartBunch.cpp?ref_type=heads#L376


    ippl::ParticleAttrib<T>* Q               = &this->pcontainer_m->Q;
    typename Base::particle_position_type* R = &this->pcontainer_m->R;

    this->fcontainer_m->getRho()             = 0.0;
    Field_t<Dim>* rho                        = &this->fcontainer_m->getRho();

    scatter(*Q, *rho, *R); /// \todo replace with scatterCIC? --> later with scatterPerBin!

#ifdef doDEBUG
    const double qtot                        = this->qi_m * this->getTotalNum();
    size_type TotalParticles                 = 0;
    size_type localParticles                 = this->pcontainer_m->getLocalNum();
   
    double relError                          = std::fabs((qtot - (*rho).sum()) / qtot);
    
    ippl::Comm->reduce(localParticles, TotalParticles, 1, std::plus<size_type>());
    
    if ((ippl::Comm->rank() == 0) && (relError > 1.0E-13)) {
            Inform m("computeSelfFields w CICScatter ", INFORM_ALL_NODES);
            m << "Time step: " << it_m
              << " total particles in the sim. " << totalP_m 
              << " missing : " << totalP_m-TotalParticles 
              << " rel. error in charge conservation: " << relError << endl;
    }
#endif
    
   /*

     scatterCIC end

    */

    this->fsolver_m->runSolver();    
    
    gather(this->pcontainer_m->E, this->fcontainer_m->getE(), this->pcontainer_m->R);

#ifdef doDEBUG
    Inform m("computeSelfFields w CICScatter ", INFORM_ALL_NODES);
    double cellVolume = std::reduce(hr_m.begin(), hr_m.end(), 1., std::multiplies<double>());
    m << "cellVolume= " << cellVolume << endl;
    m << "Sum over E-field after gather = " << this->fcontainer_m->getE().sum() << endl;
#endif

    /*
    Vector_t<double, 3> efScale = Vector_t<double,3>(gammaz*cc/hr_scaled[0], gammaz*cc/hr_scaled[1], cc / gammaz / hr_scaled[2]);
    m << "efScale = " << efScale << endl;    
    spaceChargeEFieldCheck(efScale);
    */

    //IpplTimings::stopTimer(SolveTimer);
}

template <typename T, unsigned Dim>
void PartBunch<T,Dim>::scatterCICPerBin(PartBunch<T,Dim>::binIndex_t binIndex) {
    /**
     * Scatters only particles in bin binIndex. Scatters all particles if binIndex=-1
     */
    Inform m("scatterCICPerBin");
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


// Explicit instantiations
template class PartBunch<double, 3>;
