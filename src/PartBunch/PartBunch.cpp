#include "PartBunch/PartBunch.hpp"
#include <boost/numeric/ublas/io.hpp>
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
            centroid_m[i] = 0.0;
            for (unsigned j = 0; j <= i; j++) {
                loc_moment[i][j] = 0.0;
                loc_moment[j][i] = 0.0;
                moments_m(i,j)   = 0.0;
                moments_m(j,i)   = 0.0;

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

        for (unsigned int i=0; i<2*Dim; i++) {
            centroid_m[i] = centroid[i];
            for (unsigned int j=0; j<2*Dim; j++) {
                moments_m(i,j) = moment[i][j];
            }
        }

        ippl::Comm->barrier();
    }

template<>
Inform& PartBunch<double,3>::print(Inform& os) {
        // if (this->getLocalNum() != 0) {  // to suppress Nans
        Inform::FmtFlags_t ff = os.flags();

        os << std::scientific;
        os << level1 << "\n";
        os << "* ************** B U N C H "
              "********************************************************* \n";
        os << "* PARTICLES       = " << this->getTotalNum() << "\n";
        os << "* CHARGE          = " << this->totalQ_m << "\n";
        os << "* CORES           = " << ippl::Comm->size() << "\n";
        os << "* FIELD SOLVER    = " << solver_m << "\n";
        os << "* INTEGRATOR      = " << integration_method_m << "\n";
        os << "* MIN R (origin)  = " << this->get_origin() << "\n";
        os << "* MAX R (max ext) = " << this->get_maxExtent() << "\n";
        os << "* RMS R           = " << this->get_rrms() << "\n";
        os << "* RMS P           = " << this->get_prms() << "\n";
        os << "* MESH SPACING    = " << this->fcontainer_m->getMesh().getMeshSpacing() << "\n";
        os << "* FIELD LAYOUT    = " << this->fcontainer_m->getFL() << "\n";
        os << "* "
              "********************************************************************************"
              "** "
           << endl;
        os.flags(ff);
        // }
        return os;
    }

template <>
void PartBunch<double,3>::resetFieldSolver() {

    auto *mesh = &this->fcontainer_m->getMesh();
    auto *FL = &this->fcontainer_m->getFL();

    // assume o < 0.0
    ippl::Vector<double, 3> o = this->get_origin();

    auto Rview  = this->getParticleContainer()->R.getView();
    Kokkos::parallel_for(
                         "shift r to 0,0", ippl::getRangePolicy(Rview),
                         KOKKOS_LAMBDA(const int i) {
                             Rview(i) = Rview(i) - o ;
                         });

    auto Qview  = this->getParticleContainer()->Q.getView();
    double q = totalQ_m / totalP_m;
    Kokkos::parallel_for(
                         "set q", ippl::getRangePolicy(Qview),
                         KOKKOS_LAMBDA(const int i) {
                             Qview(i) = q;
                         });

    ippl::Comm->barrier();
    this->calcBeamParameters();

    o = this->get_origin();

    Vector_t<double, 3> e = this->get_maxExtent(); 

    // assume e > o

    Vector_t<double, 3> l = e - o; 
    hr_m = 1.2*(l / this->nr_m);

    mesh->setMeshSpacing(hr_m);

    this->getParticleContainer()->update();        

    this->isFirstRepartition_m = true;
    this->loadbalancer_m->initializeORB(FL, mesh);
    this->loadbalancer_m->repartition(FL, mesh, this->isFirstRepartition_m);

}

template <>
void PartBunch<double,3>::computeSelfFields() {
        static IpplTimings::TimerRef SolveTimer = IpplTimings::getTimer("SolveTimer");    
        IpplTimings::startTimer(SolveTimer);
        this->par2grid();
        this->fsolver_m->runSolver();        
        // gather E / B field
        this->grid2par();
        IpplTimings::stopTimer(SolveTimer);
}

template <>
void PartBunch<double,3>::scatterCIC() {
        Inform m("scatterCIC ");

        this->fcontainer_m->getRho() = 0.0;

        ippl::ParticleAttrib<double>* q          = &this->pcontainer_m->Q;
        typename Base::particle_position_type* R = &this->pcontainer_m->R;
        Field_t<3>* rho                          = &this->fcontainer_m->getRho();
        double Q                                 = this->totalQ_m;
        Vector_t<double, 3> rmin                 = rmin_m;
        Vector_t<double, 3> rmax                 = rmax_m;
        Vector_t<double, 3> hr                   = hr_m;

        scatter(*q, *rho, *R);

        m << "final rho sum = " << (*rho).sum() << endl;

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
