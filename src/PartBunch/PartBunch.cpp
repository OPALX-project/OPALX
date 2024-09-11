#include "PartBunch/PartBunch.hpp"
#include <boost/numeric/ublas/io.hpp>
#include "Utilities/Util.h"

#undef doDEBUG


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
void PartBunch<double,3>::spaceChargeEFieldCheck(Vector_t<double, 3> efScale) {

 Inform msg("EParticleStats");

 auto pE_view   = this->pcontainer_m->E.getView();
 auto fphi_view = this->fcontainer_m->getPhi().getView();


 double avgE          = 0.0;
 double avgphi        = 0.0;

 int myRank = ippl::Comm->rank();

 size_t nloc = this->getLocalNum();
 
 Kokkos::parallel_reduce(
                         "check e-field", nloc,
                         KOKKOS_LAMBDA(const size_t i, double& loc_avgE) {
                             double EX    = pE_view[i][0]*efScale[0];
                             double EY    = pE_view[i][1]*efScale[1];
                             double EZ    = pE_view[i][2]*efScale[2];

                             //                             if ((i<5) || (i>nloc-4))
                             //    std::cout << "E[" << i << "]= (" << EX << "," << EY << "," << EZ << ")" << std::endl;
                             
                             double ENorm = Kokkos::sqrt(EX*EX + EY*EY + EZ*EZ);
                             
                             loc_avgE += ENorm;

                         },
                         Kokkos::Sum<double>(avgE));
 
 MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &avgE, &avgE, 1, MPI_DOUBLE, MPI_SUM, 0, ippl::Comm->getCommunicator());
 avgE /= this->getTotalNum();
 msg << "avgENorm = " << avgE << endl;

 using mdrange_type             = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
 Kokkos::parallel_reduce(
                         "check phi", mdrange_type({0,0,0}, {fphi_view.extent(0),fphi_view.extent(1),fphi_view.extent(2)}),
                         KOKKOS_LAMBDA(const int i, const int j, const int k, double& loc_avgphi) {
                             double phi = fphi_view(i,j,k);
                             loc_avgphi += phi;
                         },
                         Kokkos::Sum<double>(avgphi));

 MPI_Reduce(myRank == 0 ? MPI_IN_PLACE : &avgphi, &avgphi, 1, MPI_DOUBLE, MPI_SUM, 0, ippl::Comm->getCommunicator());
 avgphi /= this->getTotalNum(); 
 msg << "avgphi = " << avgphi << endl;
 
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
void PartBunch<double,3>::bunchUpdate(ippl::Vector<double, 3> hr) {
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

    ippl::Vector<double, 3> o = pc->getMinR();
    ippl::Vector<double, 3> e = pc->getMaxR();
    ippl::Vector<double, 3> l = e - o;
    
    mesh->setMeshSpacing(hr);
    
    mesh->setOrigin(o-0.5*hr*this->OPALFieldSolver_m->getBoxIncr()/100.);
    
    pc->getLayout().updateLayout(*FL, *mesh);
    pc->update();

    this->getFieldContainer()->setRMin(o);
    this->getFieldContainer()->setRMax(e);
    this->getFieldContainer()->setHr(hr);
    
    this->isFirstRepartition_m = true;
    this->loadbalancer_m->initializeORB(FL, mesh);

    // \fixme with the OPEN solver repartion does not work.
    // this->loadbalancer_m->repartition(FL, mesh, this->isFirstRepartition_m);
    this->updateMoments();
    


}

template <>
void PartBunch<double,3>::bunchUpdate() {

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
    
    mesh->setMeshSpacing(hr_m);
    mesh->setOrigin(o-0.5*hr_m);

    this->getParticleContainer()->getLayout().updateLayout(*FL, *mesh);
    this->getParticleContainer()->update();

    this->isFirstRepartition_m = true;
    this->loadbalancer_m->initializeORB(FL, mesh);

    // \fixme with the OPEN solver repartion does not work.
    // this->loadbalancer_m->repartition(FL, mesh, this->isFirstRepartition_m);
    this->updateMoments();
}


template <>
void PartBunch<double,3>::computeSelfFields() {
    Inform m("computeSelfFields w CICScatter ");
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

    this->fcontainer_m->getRho() = 0.0;
    size_type TotalParticles     = 0;
    size_type localParticles     = this->pcontainer_m->getLocalNum();
    const double qtot            = this->qi_m * this->getTotalNum();
    
    ippl::ParticleAttrib<double>* Q          = &this->pcontainer_m->Q;
    typename Base::particle_position_type* R = &this->pcontainer_m->R;
    Field_t<3>* rho                          = &this->fcontainer_m->getRho();

    scatter(*Q, *rho, *R);

    double relError = std::fabs((qtot - (*rho).sum()) / qtot);
    
    ippl::Comm->reduce(localParticles, TotalParticles, 1, std::plus<size_type>());
    
    if ((ippl::Comm->rank() == 0) && (relError > 1.0E-13)) {
            m << "Time step: " << it_m
              << " total particles in the sim. " << totalP_m 
              << " missing : " << totalP_m-TotalParticles 
              << " rel. error in charge conservation: " << relError << endl;
    }

   /*

     scatterCIC end

    */

    this->fsolver_m->runSolver();    
    
    gather(this->pcontainer_m->E, this->fcontainer_m->getE(), this->pcontainer_m->R);

#ifdef doDEBUG
    double cellVolume = std::reduce(hr_m.begin(), hr_m.end(), 1., std::multiplies<double>());
    m << "cellVolume= " << cellVolume << endl;
    m << "Sum over E-field after gather = " << this->fcontainer_m->getE().sum() << endl;
#endif
    
    /*
      const double cc = getCouplingConstant();
      Vector_t<double, 3> efScale = Vector_t<double,3>(gammaz*cc/hr_scaled[0], gammaz*cc/hr_scaled[1], cc / gammaz / hr_scaled[2]);
      m << "efScale = " << efScale << endl;    
      spaceChargeEFieldCheck(efScale);
    */
    
    IpplTimings::stopTimer(SolveTimer);
}

template <>
void PartBunch<double,3>::scatterCIC() {

    Inform m("scatterCIC ");
    this->fcontainer_m->getRho() = 0.0;

        ippl::ParticleAttrib<double>* q          = &this->pcontainer_m->Q;
        typename Base::particle_position_type* R = &this->pcontainer_m->R;
        Field_t<3>* rho                          = &this->fcontainer_m->getRho();
        double Q                                 = this->qi_m * this->getTotalNum();
        Vector_t<double, 3> rmin                 = rmin_m;
        Vector_t<double, 3> rmax                 = rmax_m;
        Vector_t<double, 3> hr                   = hr_m;

        scatter(*q, *rho, *R);

        // m << "gammz= " << this->pcontainer_m->getMeanP()[2] << endl;
        
        double relError = std::fabs((Q - (*rho).sum()) / Q);
        size_type TotalParticles = 0;
        size_type localParticles = this->pcontainer_m->getLocalNum();

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
        
        m << "cellVolume= " << cellVolume << endl;

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
