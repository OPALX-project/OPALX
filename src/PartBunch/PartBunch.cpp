#include "PartBunch/PartBunch.hpp"
#include <boost/numeric/ublas/io.hpp>
#include "Utilities/Util.h"

extern Inform* gmsg;

template<>
void PartBunch<double,3>::setSolver(std::string solver) {

        if (this->solver_m != "")
            *gmsg << "* Warning solver already initiated but overwrite ..." << endl;

        this->solver_m = solver;

        this->fcontainer_m->initializeFields(this->solver_m);

        this->setFieldSolver(std::make_shared<FieldSolver_t>(
            this->solver_m, &this->fcontainer_m->getRho(), &this->fcontainer_m->getE(),
            &this->fcontainer_m->getPhi()));

        this->fsolver_m->initSolver();
        
        /// ADA we need to be able to set a load balancer when not having a field solver
        this->setLoadBalancer(std::make_shared<LoadBalancer_t>(
            this->lbt_m, this->fcontainer_m, this->pcontainer_m, this->fsolver_m));

        *gmsg << "* Solver and Load Balancer set" << endl;
}

template<>
void PartBunch<double,3>::spaceChargeEFieldCheck() {

 Inform msg("EParticleStats");

 auto pE_view = this->pcontainer_m->E.getView();
 
 double avgE          = 0.0;
 double minEComponent = std::numeric_limits<double>::max();
 double maxEComponent = std::numeric_limits<double>::min();
 double minE          = std::numeric_limits<double>::max();
 double maxE          = std::numeric_limits<double>::min();
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
                         Kokkos::Sum<double>(avgE), Kokkos::Min<double>(minEComponent),
                         Kokkos::Max<double>(maxEComponent), Kokkos::Min<double>(minE),
                         Kokkos::Max<double>(maxE));

 MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &avgE, &avgE, 1, MPI_DOUBLE, MPI_SUM, 0, ippl::Comm->getCommunicator());
 MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &minEComponent, &minEComponent, 1, MPI_DOUBLE, MPI_MIN, 0, ippl::Comm->getCommunicator());
 MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &maxEComponent, &maxEComponent, 1, MPI_DOUBLE, MPI_MAX, 0, ippl::Comm->getCommunicator());
 MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &minE, &minE, 1, MPI_DOUBLE, MPI_MIN, 0, ippl::Comm->getCommunicator());
 MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &maxE, &maxE, 1, MPI_DOUBLE, MPI_MAX, 0, ippl::Comm->getCommunicator());
 
 avgE /= this->getTotalNum();
 msg << "avgENorm = " << avgE << endl;
 msg << "minEComponent = " << minEComponent << endl;
 msg << "maxEComponent = " << maxEComponent << endl;
 msg << "minE = " << minE << endl;
 msg << "maxE = " << maxE << endl;
}

template<>
void PartBunch<double,3>::calcBeamParameters() {
        // Usefull constants
        const int Dim = 3;

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
                Kokkos::Sum<double>(loc_centroid[i]), Kokkos::Sum<double>(loc_moment[i][0]),
                Kokkos::Sum<double>(loc_moment[i][1]), Kokkos::Sum<double>(loc_moment[i][2]),
                Kokkos::Sum<double>(loc_moment[i][3]), Kokkos::Sum<double>(loc_moment[i][4]),
                Kokkos::Sum<double>(loc_moment[i][5]));
            Kokkos::fence();
        }
        ippl::Comm->barrier();

        MPI_Allreduce(
            loc_moment, moment, 2 * Dim * 2 * Dim, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
        MPI_Allreduce(
            loc_centroid, centroid, 2 * Dim, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());

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
                Kokkos::Max<double>(rmax_loc[d]));

            Kokkos::parallel_reduce(
                "rel min", this->getLocalNum(),
                KOKKOS_LAMBDA(const int i, double& mm) {
                    double tmp_vel = Rview(i)[d];
                    mm             = tmp_vel < mm ? tmp_vel : mm;
                },
                Kokkos::Min<double>(rmin_loc[d]));
        }
        Kokkos::fence();
        MPI_Allreduce(rmax_loc, rmax, Dim, MPI_DOUBLE, MPI_MAX, ippl::Comm->getCommunicator());
        MPI_Allreduce(rmin_loc, rmin, Dim, MPI_DOUBLE, MPI_MIN, ippl::Comm->getCommunicator());

        // \todo can we do this nicer? 
        for (unsigned int i=0; i<Dim; i++) {
            rmax_m(i) = rmax[i];
            rmin_m(i) = rmin[i];
        }
        ippl::Comm->barrier();
    }

template<>
void PartBunch<double,3>::pre_run() {
    this->fcontainer_m->getRho() = 0.0;
    this->fsolver_m->runSolver();
    *gmsg << "* Solver warmup done" << endl;
}

template<>
Inform& PartBunch<double,3>::print(Inform& os) {
        // if (this->getLocalNum() != 0) {  // to suppress Nans
        Inform::FmtFlags_t ff = os.flags();
        const int Dim = 3;

        os << std::scientific;
        os << level1 << "\n";
        os << "* ************** B U N C H "
              "********************************************************* \n";
        os << "* PARTICLES       = " << this->getTotalNum() << "\n";
        os << "* CHARGE          = " << this->qi_m*this->getTotalNum() << " (Cb) \n";
        os << "* CORES           = " << ippl::Comm->size() << "\n";
        os << "* FIELD SOLVER    = " << solver_m << "\n";
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
        // }
        return os;
    }

template <>
void PartBunch<double,3>::bunchUpdate() {

    /* \brief
      1. calculates and set hr
      2. do repartitioning
     */


    auto *mesh = &this->fcontainer_m->getMesh();
    auto *FL   = &this->fcontainer_m->getFL();

    std::shared_ptr<ParticleContainer_t> pc = this->getParticleContainer();
    pc->computeMinMaxR();
    /* 
       This needs to go 
    auto Rview  = this->getParticleContainer()->R.getView();
    Kokkos::parallel_for(
                         "shift r to 0,0", ippl::getRangePolicy(Rview),
                         KOKKOS_LAMBDA(const int i) {
                             Rview(i) = Rview(i) - o ;
                         });
    */

    /// \brief assume o < 0.0?

    ippl::Vector<double, 3> o = pc->getMinR();
    ippl::Vector<double, 3> e = pc->getMaxR();
    ippl::Vector<double, 3> l = e - o; 
    hr_m = (1.0+this->OPALFieldSolver_m->getBoxIncr()/100.)*(l / this->nr_m);
    mesh->setMeshSpacing(hr_m);
    mesh->setOrigin(o-0.5*hr_m*this->OPALFieldSolver_m->getBoxIncr()/100.);
    
    pc->getLayout().updateLayout(*FL, *mesh);
    pc->update();

    this->getFieldContainer()->setRMin(o);
    this->getFieldContainer()->setRMax(e);
    this->getFieldContainer()->setHr(hr_m);
    
    /* old Tracker 
    this->calcBeamParameters();


    ippl::Vector<double, 3> o = this->get_origin();
    ippl::Vector<double, 3> e = this->get_maxExtent(); 
    ippl::Vector<double, 3> l = e - o; 

    hr_m = (1.0+this->OPALFieldSolver_m->getBoxIncr()/100.)*(l / this->nr_m);

    mesh->setMeshSpacing(hr_m);
    mesh->setOrigin(o-0.5*hr_m);

    this->getParticleContainer()->getLayout().updateLayout(*FL, *mesh);
    this->getParticleContainer()->update();
    */

    this->isFirstRepartition_m = true;
    this->loadbalancer_m->initializeORB(FL, mesh);

    // \fixme with the OPEN solver repartion does not work.
    // this->loadbalancer_m->repartition(FL, mesh, this->isFirstRepartition_m);
    this->updateMoments();
}

/*

    auto rhoView = rho.getView();
    using index_array_type_3d = typename ippl::RangePolicy<Dim>::index_array_type;
    const index_array_type_3d& args3d = {8,8,8};
    auto res3d = ippl::apply(rhoView,args3d);
    m << "Type of variable: " << typeid(res3d).name() << " rho[8,8,8]= " << res3d << endl;

    double qsum = 0.0;
    for (int i=0; i<nr_m[0]; i++) {
        for (int j=0; j<nr_m[1]; j++) {
            for (int k=0; k<nr_m[2]; k++) {
                const index_array_type_3d& args3d = {i,j,k};
                qsum += ippl::apply(rhoView,args3d);
            }
        }
    }

    m << "sum(rho_m)= " << qsum << endl;

    
    auto Qview =  this->pcontainer_m->Q.getView();
    using index_array_type_1d = typename ippl::RangePolicy<1>::index_array_type;
    const index_array_type_1d& args1d = {0}; 
    auto res1d = ippl::apply(Qview,args1d);
    m << "Type of variable: " << typeid(res1d).name() << " value " << res1d << endl;


    
 */



// ADA
template <>
void PartBunch<double,3>::computeSelfFields() {
    // Start binning and sorting
    std::shared_ptr<AdaptBins_t> bins = this->getBins();
    bins->doFullRebin(10); // rebin with 10 bins
    bins->sortContainerByBin(); 
    bins->print(); // For debugging...


    Inform m("computeSelfFields ");
    static IpplTimings::TimerRef SolveTimer = IpplTimings::getTimer("SolveTimer");
    IpplTimings::startTimer(SolveTimer);

    m << "Running binned solver routine." << endl;

    // Run solver for each bin
    this->Etmp_m = 0.0; // reset temporary field to zero
    for (binIndex_t i = 0; i < bins->getCurrentBinCount(); ++i) {
        // Scatter only for current bin index
        this->scatterCIC(i);

        // Run solver: obtains phi_m only for what was scattered in the previous step
        this->fsolver_m->runSolver();
        this->Etmp_m = this->Etmp_m + bins->LTrans(this->fcontainer_m->getE(), i);
    }

    // Gather built up temporary E field to the particles
    gather(this->pcontainer_m->E, this->Etmp_m, this->pcontainer_m->R);

    spaceChargeEFieldCheck();
    IpplTimings::stopTimer(SolveTimer);


    /**
     * Didn't change anything after this for the binning, might be nice for debug output.
     * However, not necessary once the binning is working. 
     * Or might need to be modified...
     */
    Field_t<3>& rho                          = this->fcontainer_m->getRho();
    auto rhoView                             = rho.getView();
    ippl::ParticleAttrib<double>& Q          = this->pcontainer_m->Q;
    typename Base::particle_position_type& R = this->pcontainer_m->R;

    rho = 0.0;

    Q           = Q*this->pcontainer_m->dt;
    this->qi_m  = this->qi_m * getdT();

    scatter(Q, rho, R);

    Q           = Q/this->pcontainer_m->dt;
    this->qi_m  = this->qi_m / getdT();
    rho         = rho/getdT();
    
    double gammaz = this->pcontainer_m->getMeanGammaZ();
    double scaleFactor = 1;
    // double scaleFactor = Physics::c * getdT();
    //and get meshspacings in real units [m]

    Vector_t<double, 3> hr_scaled = hr_m * scaleFactor;
    hr_scaled[2] *= gammaz;


    double tmp2 = 1 / hr_scaled[0] * 1 / hr_scaled[1] * 1 / hr_scaled[2];
    
    using index_array_type = typename ippl::RangePolicy<Dim>::index_array_type;

    //divide charge by a 'grid-cube' volume to get [C/m^3]
    rho = rho *tmp2;
        
    double Npoints = nr_m[0] * nr_m[1] * nr_m[2];

    double localDensity2 = 0.0, Density2=0.0;
    ippl::parallel_reduce(
                          "Density stats", ippl::getRangePolicy(rhoView),
                          KOKKOS_LAMBDA(const index_array_type& args, double& den) {
                              double val = ippl::apply(rhoView, args);
                              den  += Kokkos::pow(val, 2);
                          },
                          Kokkos::Sum<double>(localDensity2) );
    
    ippl::Comm->reduce(localDensity2, Density2, 1, std::plus<double>());
    
    rmsDensity_m = std::sqrt( (1.0 /Npoints) * Density2 / Physics::q_e / Physics::q_e );

    this->calcDebyeLength();

    m << "gammaz = "  << gammaz << endl;
    m << "hr_scaled = " << hr_scaled << endl;
    m << "coupling constant= " << getCouplingConstant() << endl;
    
    //this->fsolver_m->runSolver();    

    //gather(this->pcontainer_m->E, this->fcontainer_m->getE(), this->pcontainer_m->R);

    //spaceChargeEFieldCheck();

    //IpplTimings::stopTimer(SolveTimer);
}

template <>
void PartBunch<double,3>::scatterCIC(binIndex_t binIndex = -1) {
    /**
     * Scatters only particles in bin binIndex. Scatters all particles if binIndex=-1
     */
    Inform m("scatterCIC ");
    m << "Scattering binIndex = " << binIndex << " to grid." << endl;

    this->fcontainer_m->getRho() = 0.0;

    ippl::ParticleAttrib<double>* q          = &this->pcontainer_m->Q;
    typename Base::particle_position_type* R = &this->pcontainer_m->R;
    Field_t<3>* rho                          = &this->fcontainer_m->getRho();
    double Q                                 = this->qi_m * this->bins_m->getNPartInBin(binIndex, true);
    Vector_t<double, 3> rmin                 = rmin_m;
    Vector_t<double, 3> rmax                 = rmax_m;
    Vector_t<double, 3> hr                   = hr_m;

    if (binIndex == -1) {
        scatter(*q, *rho, *R);
    } else {
        scatter(*q, *rho, *R, this->bins_m->getBinIterationPolicy(binIndex), this->bins_m->getHashArray());
    }

    m << "gammz= " << this->pcontainer_m->getMeanP()[2] << endl;
    
    double relError = std::fabs((Q - (*rho).sum()) / Q);
    size_type TotalParticles = 0;
    size_type localParticles = this->pcontainer_m->getLocalNum();

    m << "computeSelfFields sum rho = " << (*rho).sum() << ", relError = " << relError << endl;
    
    ippl::Comm->reduce(localParticles, TotalParticles, 1, std::plus<size_type>());

    if (ippl::Comm->rank() == 0) {
        if (TotalParticles != totalP_m || relError > 1e-10) {
            m << "Time step: " << it_m << endl;
            m << "Total particles in the sim. " << totalP_m << " "
                << "after update: " << TotalParticles << endl;
            m << "Rel. error in charge conservation: " << relError << endl;
            ippl::Comm->abort();
        }
    }

    double cellVolume = std::reduce(hr.begin(), hr.end(), 1., std::multiplies<double>());
    (*rho)            = (*rho) / cellVolume;

    // double rhoNorm = norm(*rho);
    // rho = rho_e - rho_i (only if periodic BCs)
    if (this->fsolver_m->getStype() != "OPEN") {
        double size = 1;
        for (unsigned d = 0; d < 3; d++) {
            size *= rmax[d] - rmin[d];
        }
        *rho = *rho - (Q / size);
    }
}
