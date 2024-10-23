#ifndef IPPL_GAUSSIAN_H
#define IPPL_GAUSSIAN_H

#include "Distribution.h"
#include "SamplingBase.hpp"
#include <memory>

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t = FieldContainer<double, 3>;
using Distribution_t = Distribution;
using GeneratorPool = typename Kokkos::Random_XorShift64_Pool<>;
using Dist_t = ippl::random::NormalDistribution<double, 3>;

class Gaussian : public SamplingBase {
public:
    Gaussian(std::shared_ptr<ParticleContainer_t> &pc, std::shared_ptr<FieldContainer_t> &fc, std::shared_ptr<Distribution_t> &opalDist)
        : SamplingBase(pc, fc, opalDist) {}

    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override {

        size_t randInit;
        if (Options::seed == -1) {
            randInit = 1234567;
            *gmsg << "* Seed = " << randInit << " on all ranks" << endl;
        }
        else
            randInit = (size_t)(Options::seed + 100 * ippl::Comm->rank());
                
        GeneratorPool rand_pool64(randInit);
        
        double mu[3], sd[3];

        // sample R
        Vector_t<double, 3> rmin = -opalDist_m->getCutoffR();
        Vector_t<double, 3> rmax =  opalDist_m->getCutoffR();
        for(int i=0; i<3; i++){
            mu[i] = 0.0;
       	    sd[i] = opalDist_m->getSigmaR()[i];
            rmin(i) = (rmin(i) + mu[i])*opalDist_m->getSigmaR()[i];
            rmax(i) = (rmax(i) + mu[i])*opalDist_m->getSigmaR()[i];
        }
        view_type &Rview = pc_m->R.getView();

        const double par[6] = {mu[0], sd[0], mu[1], sd[1], mu[2], sd[2]};
        using Dist_t = ippl::random::NormalDistribution<double, Dim>;
        using sampling_t = ippl::random::InverseTransformSampling<double, Dim, Kokkos::DefaultExecutionSpace, Dist_t>;
        Dist_t dist(par);

        MPI_Comm comm = MPI_COMM_WORLD;
        int nranks;
        int rank;
        MPI_Comm_size(comm, &nranks);
        MPI_Comm_rank(comm, &rank);

        size_type nlocal = floor(numberOfParticles/nranks);

        // if nlocal*nranks > numberOfParticles, put the remaining in rank 0
        size_t remaining = numberOfParticles - nlocal*nranks;
        if(remaining>0 && rank==0){
            nlocal += remaining;
        }
        sampling_t sampling(dist, rmax, rmin, rmax, rmin, nlocal);
        nlocal = sampling.getLocalSamplesNum();
        pc_m->create(nlocal);
        sampling.generate(Rview, rand_pool64);

        // sample P
        for(int i=0; i<3; i++){
            mu[i] = 0.0;
            sd[i] = opalDist_m->getSigmaP()[i];
        }
        view_type &Pview = pc_m->P.getView();
        Kokkos::parallel_for(
            nlocal, ippl::random::randn<double, 3>(Pview, rand_pool64, mu, sd)
        );
        Kokkos::fence();
        ippl::Comm->barrier();

        // correct the mean
        double avrgpz = opalDist_m->getAvrgpz();
        Kokkos::parallel_for(
            nlocal, KOKKOS_LAMBDA(
                    const int k) {
                    Pview(k)[2] += avrgpz;
            }
        );
        Kokkos::fence();
        ippl::Comm->barrier();
    }
};
#endif
