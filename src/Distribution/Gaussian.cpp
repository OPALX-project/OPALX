#include "Distribution.h"
#include "SamplingBase.hpp"
#include "Gaussian.h"
#include <memory>
#include <cmath>

/**
 * @brief Constructs a Gaussian sampler.
 *
 * @param pc Shared pointer to the particle container.
 * @param fc Shared pointer to the field container.
 * @param opalDist Shared pointer to the distribution object.
 */
Gaussian::Gaussian(std::shared_ptr<ParticleContainer_t> &pc, 
                   std::shared_ptr<FieldContainer_t> &fc, 
                   std::shared_ptr<Distribution_t> &opalDist)
    : SamplingBase(pc, fc, opalDist) {
    samperTimer_m = IpplTimings::getTimer("SamplingTimer");
}

/**
 * @brief Generates particles following a Gaussian distribution.
 *
 * @param numberOfParticles The total number of particles to generate.
 * @param nr The number of grid points in each dimension (not used here).
 */
void Gaussian::generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) {
    extern Inform* gmsg;

    IpplTimings::startTimer(samperTimer_m);

    size_t randInit;

    if (Options::seed == -1) {
        randInit = 1234567;
        *gmsg << "* Seed = " << randInit << " on all ranks" << endl;
    } else {
        randInit = static_cast<size_t>(Options::seed + 100 * ippl::Comm->rank());
    }

    GeneratorPool rand_pool64(randInit);
    
    double mu[3], sd[3];

    Vector_t<double, 3> rmin = -opalDist_m->getCutoffR();
    Vector_t<double, 3> rmax = opalDist_m->getCutoffR();

    for (int i = 0; i < 3; i++) {
        mu[i] = 0.0;
        sd[i] = opalDist_m->getSigmaR()[i];
        rmin(i) = (rmin(i) + mu[i]) * opalDist_m->getSigmaR()[i];
        rmax(i) = (rmax(i) + mu[i]) * opalDist_m->getSigmaR()[i];
    }

    view_type &Rview = pc_m->R.getView();
    const double par[6] = {mu[0], sd[0], mu[1], sd[1], mu[2], sd[2]};

    using Dist_t = ippl::random::NormalDistribution<double, 3>;
    using sampling_t = ippl::random::InverseTransformSampling<double, 3, Kokkos::DefaultExecutionSpace, Dist_t>;
    Dist_t dist(par);

    MPI_Comm comm = MPI_COMM_WORLD;
    int nranks, rank;
    MPI_Comm_size(comm, &nranks);
    MPI_Comm_rank(comm, &rank);

    size_t nlocal = floor(numberOfParticles / nranks);
    size_t remaining = numberOfParticles - nlocal * nranks;

    if (remaining > 0 && rank == 0) {
        nlocal += remaining;
    }

    sampling_t sampling(dist, rmax, rmin, rmax, rmin, nlocal);
    nlocal = sampling.getLocalSamplesNum();
    pc_m->create(nlocal);
    sampling.generate(Rview, rand_pool64);

    double meanR[3] = {0.0, 0.0, 0.0}, loc_meanR[3] = {0.0, 0.0, 0.0};

    Kokkos::parallel_reduce(
        "calc moments of particle distr.", numberOfParticles,
        KOKKOS_LAMBDA(const int k, double &cent0, double &cent1, double &cent2) {
            cent0 += Rview(k)[0];
            cent1 += Rview(k)[1];
            cent2 += Rview(k)[2];
        },
        Kokkos::Sum<double>(loc_meanR[0]), Kokkos::Sum<double>(loc_meanR[1]), Kokkos::Sum<double>(loc_meanR[2])
    );
    Kokkos::fence();

    MPI_Allreduce(loc_meanR, meanR, 3, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
    ippl::Comm->barrier(); 

    for (int i = 0; i < 3; i++) {
        meanR[i] = meanR[i] / static_cast<double>(numberOfParticles);
    }

    Kokkos::parallel_for(
        numberOfParticles,
        KOKKOS_LAMBDA(const int k) {
            Rview(k)[0] -= meanR[0];
            Rview(k)[1] -= meanR[1];
            Rview(k)[2] -= meanR[2];
        }
    );
    Kokkos::fence();

    for (int i = 0; i < 3; i++) {
        mu[i] = 0.0;
        sd[i] = opalDist_m->getSigmaP()[i];
    }

    view_type &Pview = pc_m->P.getView();
    Kokkos::parallel_for(
        nlocal,
        ippl::random::randn<double, 3>(Pview, rand_pool64, mu, sd)
    );
    Kokkos::fence();

    double avrgpz = opalDist_m->getAvrgpz();
    Kokkos::parallel_for(
        nlocal,
        KOKKOS_LAMBDA(const int k) {
            Pview(k)[2] += avrgpz;
        }
    );
    Kokkos::fence();

    IpplTimings::stopTimer(samperTimer_m);
}

