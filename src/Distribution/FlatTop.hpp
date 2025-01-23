#ifndef IPPL_FLAT_TOP_H
#define IPPL_FLAT_TOP_H

#include "Distribution.h"
#include "SamplingBase.hpp"
#include <Kokkos_Sort.hpp>
#include <memory>
#include <cmath>

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t = FieldContainer<double, 3>;
using Distribution_t = Distribution;
using GeneratorPool = typename Kokkos::Random_XorShift64_Pool<>;
using Dist_t = ippl::random::NormalDistribution<double, 3>;

class FlatTop : public SamplingBase {
public:
    FlatTop(std::shared_ptr<ParticleContainer_t> &pc, std::shared_ptr<FieldContainer_t> &fc, std::shared_ptr<Distribution_t> &opalDist)
        : SamplingBase(pc, fc, opalDist), rand_pool_m(determineRandInit()) {
            setParameters(opalDist);
        }

private:
    using size_type = ippl::detail::size_type;
    GeneratorPool rand_pool_m;
    double flattopTime_m;
    double normalizedFlankArea_m;
    double distArea_m;
    double sigmaTFall_m;
    double sigmaTRise_m;
    Vector_t<double, 3> cutoffR_m;
    double fallTime_m;
    double riseTime_m;
    bool emitting_m;
    size_type totalN_m;

    static size_t determineRandInit() {
        size_t randInit;
        if (Options::seed == -1) {
            randInit = 1234567;
            *gmsg << "* Seed = " << randInit << " on all ranks" << endl;
        } else {
            randInit = static_cast<size_t>(Options::seed + 100 * ippl::Comm->rank());
        }
        return randInit;
    }

    void setParameters(const std::shared_ptr<Distribution_t> &opalDist) {
        emitting_m = opalDist->emitting_m;
        // time span of fall is [0, riseTime, riseTime+flattopTime, fallTime+flattopTime+riseTime ]
        sigmaTFall_m = opalDist_m->getSigmaTFall();
        sigmaTRise_m = opalDist_m->getSigmaTRise();
        cutoffR_m = opalDist_m->getCutoffR();

        fallTime_m = sigmaTFall_m * cutoffR_m[2]; // fall is [0, fallTime]
        flattopTime_m = opalDist->getTPulseLengthFWHM()
            - std::sqrt(2.0 * std::log(2.0)) * (sigmaTRise_m + sigmaTFall_m);
        if (flattopTime_m < 0.0) {
            flattopTime_m = 0.0;
        }
        riseTime_m = sigmaTRise_m * cutoffR_m[2];

        // These expression are take from the old OPAL
        // I think normalizedFlankArea is int_0^{cutoff} exp(-(x/sigma)^2/2 ) / sigma
        // Instead of int_0^{cutoff} exp(-(x/sigma)^2/2 ) / sqrt(2*pi) / sigma, which is strange!
        // So the distribution of tails are exp(-(x/sigma)^2/2 ) and not Gaussian!
        normalizedFlankArea_m = 0.5 * std::sqrt(Physics::two_pi) * std::erf(cutoffR_m[2] / std::sqrt(2.0));
        distArea_m = flattopTime_m + (sigmaTRise_m + sigmaTFall_m) * normalizedFlankArea_m;
    }

public:
    void generateUniformDisk(size_type nlocal, size_t nNew) {

        GeneratorPool rand_pool = rand_pool_m;
        Vector_t<double, 3> rmin;
        Vector_t<double, 3> rmax;
        Vector_t<double, 3> hr;

        view_type Rview = pc_m->R.getView();
        view_type Pview = pc_m->P.getView();

        double pi = Physics::pi;
        Vector_t<double, 3> sigmaR = opalDist_m->getSigmaR();
        // Sample (Rx,Ry) on a unit ring, then scale with sigmaRx and sigmaRy, set Px=Py=0
        Kokkos::parallel_for(
               "unitDisk", Kokkos::RangePolicy<>(nlocal, nlocal+nNew), KOKKOS_LAMBDA(const size_t j) {
                auto generator = rand_pool.get_state();
                double r = Kokkos::sqrt( generator.drand(0., 1.) );
                double theta = 2.0 * pi * generator.drand(0., 1.);
                rand_pool.free_state(generator);

                Rview(j)[0] = r * Kokkos::cos(theta) * sigmaR[0];
                Rview(j)[1] = r * Kokkos::sin(theta) * sigmaR[1];
                Rview(j)[2]  = 0.0;
                Pview(j)[0] = 0.0;
                Pview(j)[1] = 0.0;
                Pview(j)[2] = 0.0;
        });
        Kokkos::fence();
    }

    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override {

        // initial allocation is similar for both emitting and non-emitting cases
        allocateParticles(numberOfParticles);

        if(emitting_m){
            // set nlocal to 0 for the very first time step, before sampling particles
            //pc_m->setLocalNum(0);
            
            /*
            Cannot use pc_m->setLocalNum(0) here, because there is no way to update totalNum_m
            in the ParticleBase class. Instead, one has to call destroy(), which is efficient, if
            all particles are destroyed. In that case it simply updates also the totalNum_m count
            and sets localNum_m to 0. 
            However, it expects a bool Kokkos::View that indicates which particles are to be deleted.
            It is never used when ALL particles are deleted, so we can just pass an empty one.

            Note: this ONLY works if ALL particles are to be deleted. This does not work if a few
            might be sampled and some space is to be reserved for the rest. In that case, we would éither need
            to implement a "reserve" function in ParticleBase/ParticleContainer or write a setter for totalNum_m.
            */ 
            Kokkos::View<bool*> tmp_invalid("tmp_invalid", 0);
            pc_m->destroy(tmp_invalid, pc_m->getLocalNum());
        }
    }

    double FlatTopProfile(double t){
        double t0;
        if(t<riseTime_m){
            t0 = riseTime_m;
            return exp( -pow((t-t0)/sigmaTRise_m,2) /2. );
            //  In the old opal, tails seem to be exp(-x^2/sigma^2/2) rather than Gaussian with normalizing factor.
        }
	else if( t>riseTime_m && t<riseTime_m + flattopTime_m){
            return 1.;
        }
	else if(t>riseTime_m + flattopTime_m && t < fallTime_m + flattopTime_m + riseTime_m){
            t0 = fallTime_m + flattopTime_m;
            return exp( -pow((t-t0)/sigmaTFall_m,2)/2. );
            //  In the old opal, tails seem to be exp(-x^2/sigma^2/2) rather than Gaussian with normalizing factor.
        }
	else
            return 0.;
    }

    size_t computeNlocal(size_t nglobal){
        MPI_Comm comm = MPI_COMM_WORLD;
        int nranks;
        int rank;
        MPI_Comm_size(comm, &nranks);
        MPI_Comm_rank(comm, &rank);

        size_type nlocal = floor(nglobal/nranks);

        size_t remaining = nglobal - nlocal*nranks;

        if(remaining>0 && rank==0){
            nlocal += remaining;
        }

        return nlocal;
    }

    double ingerateTrapezoidal(double x1, double x2, double y1, double y2){
        return 0.5 * (y1+y2) * fabs(x2-x1);
    }

    double countEnteringParticlesPerRank(double t0, double tf){
        double tArea = 0.0;

        tArea = ingerateTrapezoidal(t0, tf, FlatTopProfile(t0), FlatTopProfile(tf));

        size_type totalNew = floor(totalN_m * tArea / distArea_m);

        size_type nlocalNew = computeNlocal(totalNew);

        return nlocalNew;
    }

    void allocateParticles(size_t numberOfParticles){
        totalN_m = numberOfParticles;

        size_type nlocal;

        nlocal = computeNlocal(totalN_m);

        pc_m->create(nlocal);
    }

    void emitParticles(double t, double dt) override {
        // count number of new particles to be emitted
        size_type nNew = countEnteringParticlesPerRank(t, t + dt);
        std::cout << "New Particles = " << nNew << std::endl; 

        if(nNew > 0){
            // current number of particles per rank
            size_type nlocal = pc_m->getLocalNum();

            // extend particle container to accomodate new particles
            pc_m->create(nNew);

            // generate new particles on uniform disc
            *gmsg << "* generate particles on a disc" << endl;
            generateUniformDisk(nlocal, nNew);

            *gmsg << "* new particles emmitted" << endl;
        }
    }

    void testNumEmitParticles(size_type nsteps, double dt) override {
        size_type nNew;
        int rank, numRanks;
        double t = 0.0;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
        std::string filename = "timeNpart_" + std::to_string(rank) + ".txt";
        std::ofstream file(filename);

        for(size_type i=0; i<nsteps; i++){
            nNew = countEnteringParticlesPerRank(t, t + dt);
            file << t << "  " << nNew << "\n";
            t = t + dt;
        }
        file.close();
    }

    void testEmitParticles(size_type nsteps, double dt) override {
        double t = 0.0;

        for(size_type i=0; i<nsteps; i++){
            emitParticles(t, dt);
            t = t + dt;
        }
    }
};
#endif
