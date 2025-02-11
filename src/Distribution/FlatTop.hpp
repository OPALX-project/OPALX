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
            withDomainDecomp_m = false;
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
    bool withDomainDecomp_m;
    double emissionTime_m;
    Vector_t<double, 3> nr_m;
    Vector_t<double, 3> hr_m;

    void setWithDomainDecomp(bool withDomainDecomp) override {
        withDomainDecomp_m = withDomainDecomp;
    }

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

        emissionTime_m = fallTime_m + flattopTime_m + riseTime_m;
        opalDist_m->setTEmission(emissionTime_m);

        // These expression are take from the old OPAL
        // I think normalizedFlankArea is int_0^{cutoff} exp(-(x/sigma)^2/2 ) / sigma
        // Instead of int_0^{cutoff} exp(-(x/sigma)^2/2 ) / sqrt(2*pi) / sigma, which is strange!
        // So the distribution of tails are exp(-(x/sigma)^2/2 ) and not Gaussian!
        normalizedFlankArea_m = 0.5 * std::sqrt(Physics::two_pi) * std::erf(cutoffR_m[2] / std::sqrt(2.0));
        distArea_m = flattopTime_m + (sigmaTRise_m + sigmaTFall_m) * normalizedFlankArea_m;

        // make sure only z direction is decomposed
        fc_m->setDecomp({false, false, true});
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

    void setNr(Vector_t<double, 3> nr){
        nr_m = nr;
    }
    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override {

        setNr(nr);

        // initial allocation is similar for both emitting and non-emitting cases
        allocateParticles(numberOfParticles);

        if(emitting_m){
            // set nlocal to 0 for the very first time step, before sampling particles
            //pc_m->setLocalNum(0);

            Kokkos::View<bool*> tmp_invalid("tmp_invalid", 0);
            // \todo might be abuse of semantics: maybe think about new pc_m->setTotalNum or pc_m->updateTotal function instead?
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

    size_t computeNlocalUniformly(size_t nglobal){
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

    double integrateTrapezoidal(double x1, double x2, double y1, double y2){
        return 0.5 * (y1+y2) * fabs(x2-x1);
    }

    void initDomainDecomp(double BoxIncr) override {
        auto *mesh = &fc_m->getMesh();
        auto *FL   = &fc_m->getFL();
        Vector_t<double, 3> sigmaR = opalDist_m->getSigmaR();
        ippl::Vector<double, 3> o;
        ippl::Vector<double, 3> e;
        double tol = 1e-15; // enlarge grid by tol to avoid missing particles on boundaries
        o[0] = -sigmaR[0] - tol;
        e[0] =  sigmaR[0] + tol;
        o[1] = -sigmaR[1] - tol;
        e[1] =  sigmaR[1] + tol;
        o[2] = 0.0 - tol;
        e[2] = Physics::c * emissionTime_m + tol;

        ippl::Vector<double, 3> l = e - o;
        hr_m = (1.0+BoxIncr/100.)*(l / nr_m);
        mesh->setMeshSpacing(hr_m);
        mesh->setOrigin(o-0.5*hr_m*BoxIncr/100.);
        pc_m->getLayout().updateLayout(*FL, *mesh);
    }

    double countEnteringParticlesPerRank(double t0, double tf){
        size_type nlocalNew=0;
        double tArea = 0.0;
        tArea = integrateTrapezoidal(t0, tf, FlatTopProfile(t0), FlatTopProfile(tf));
        size_type totalNew = floor(totalN_m * tArea / distArea_m);
        nlocalNew = 0;

        if(totalNew>0){
          if(!withDomainDecomp_m){
            // the same number of particles per rank
            nlocalNew = computeNlocalUniformly(totalNew);
          }
          else{
            // select number of particles per rank using estimated domain decomposition at final emission time
            // find min/max of particle positions for [t,t+dt]
            Vector_t<double, 3> prmin, prmax;
            Vector_t<double, 3> sigmaR = opalDist_m->getSigmaR();
            prmin[0] = -sigmaR[0];
            prmax[0] =  sigmaR[0];
            prmin[1] = -sigmaR[1];
            prmax[1] =  sigmaR[1];
            prmin[2] = std::max(0.0, Physics::c * t0);
            prmax[2] = Physics::c*tf;

            double dx = prmax[0] - prmin[0];
            double dy = prmax[1] - prmin[1];
            double dz = prmax[2] - prmin[2];

            if (dx <= 0 || dy <= 0 || dz <= 0) {
                throw std::runtime_error("Invalid global particle volume: prmax must be greater than prmin.");
            }

            double globalpvolume = dx * dy * dz;

            // find min/max of subdomains for the current rank
            auto regions = pc_m->getLayout().getRegionLayout().gethLocalRegions();
            int rank = ippl::Comm->rank();

            if (rank < 0 || static_cast<size_t>(rank) >= regions.size()) {
                throw std::runtime_error("Invalid rank index in gethLocalRegions()");
            }

            Vector_t<double, 3> locrmin, locrmax;
            for (unsigned d = 0; d < Dim; ++d) {
               locrmax[d] = regions(rank)[d].max();
               locrmin[d] = regions(rank)[d].min();
            }

            if (prmax[0] >= locrmin[0] && prmin[0] <= locrmax[0] &&
                prmax[1] >= locrmin[1] && prmin[1] <= locrmax[1] &&
                prmax[2] >= locrmin[2] && prmin[2] <= locrmax[2]) {

                double x1 = std::max(prmin[0], locrmin[0]);
                double x2 = std::min(prmax[0], locrmax[0]);
                double y1 = std::max(prmin[1], locrmin[1]);
                double y2 = std::min(prmax[1], locrmax[1]);
                double z1 = std::max(prmin[2], locrmin[2]);
                double z2 = std::min(prmax[2], locrmax[2]);

                if (x2 >= x1 && y2 >= y1 && z2 >= z1) {
                    double locpvolume = (x2 - x1) * (y2 - y1) * (z2 - z1);
                    if (globalpvolume > 0) {
                        nlocalNew = static_cast<int>(totalNew * locpvolume / globalpvolume);
                    }
                    else{
                        nlocalNew = 0;
                    }
                } else {
                    nlocalNew = 0;
                }
            }
          }
        }
        return nlocalNew;
    }

    void allocateParticles(size_t numberOfParticles){
        totalN_m = numberOfParticles;

        size_type nlocal;

        nlocal = computeNlocalUniformly(totalN_m);

        pc_m->create(nlocal);
    }

    void emitParticles(double t, double dt) override {
        // count number of new particles to be emitted
        size_type nNew = countEnteringParticlesPerRank(t, t + dt);

        // current number of particles per rank
        size_type nlocal = pc_m->getLocalNum();

        // extend particle container to accomodate new particles
        pc_m->create(nNew);

        // generate new particles on uniform disc
        *gmsg << "* generate particles on a disc" << endl;
        generateUniformDisk(nlocal, nNew);

        *gmsg << "* new particles emmitted" << endl;
    }

    void testNumEmitParticles(size_type nsteps, double dt) override {
        size_type nNew;
        int rank, numRanks;
        double t = 0.0;
        double c = Physics::c;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

        std::string filename = "timeNpart_" + std::to_string(rank) + ".txt";
        std::ofstream file(filename);

        // Check if the file opened successfully
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        for(size_type i=0; i<nsteps; i++){
            nNew = countEnteringParticlesPerRank(t, t + dt);

            // current number of particles per rank
            size_type nlocal = pc_m->getLocalNum();

            // extend particle container to accomodate new particles
            pc_m->create(nNew);

            // generate new particles on uniform disc
            generateUniformDisk(nlocal, nNew);

            // write to a file
            auto rViewDevice  = pc_m->R.getView();
            auto rView = Kokkos::create_mirror_view(rViewDevice);
            Kokkos::deep_copy(rView,rViewDevice);

            for(size_type j=nlocal; j<nlocal+nNew; j++){
                file << t << " " << (emissionTime_m-t)*c << " " << rView(j)[0] << " " << rView(j)[1] << "\n";
            }
            file.flush();  // Ensure data is written to disk

            t = t + dt;
        }
        file.close();
        ippl::Comm->barrier();
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
