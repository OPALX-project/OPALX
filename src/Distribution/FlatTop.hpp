#ifndef IPPL_FLAT_TOP_H
#define IPPL_FLAT_TOP_H

#include "Distribution.h"
#include "SamplingBase.hpp"
#include <Kokkos_Sort.hpp>
//#include <Kokkos_NestedSort.hpp> 
#include <memory>

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t = FieldContainer<double, 3>;
using Distribution_t = Distribution;
using GeneratorPool = typename Kokkos::Random_XorShift64_Pool<>;
using Dist_t = ippl::random::NormalDistribution<double, 3>;

struct VectorComparator {
    KOKKOS_INLINE_FUNCTION
    bool operator()(const ippl::Vector<double, 1>& a, const ippl::Vector<double, 1>& b) const {
        // Define comparison based on specific needs, e.g., comparing the first element
        return a[0] < b[0];  // Customize comparison logic as needed
    }
};


class FlatTop : public SamplingBase {
public:
    FlatTop(std::shared_ptr<ParticleContainer_t> &pc, std::shared_ptr<FieldContainer_t> &fc, std::shared_ptr<Distribution_t> &opalDist)
        : SamplingBase(pc, fc, opalDist) {}

    void generateUniformDisk(size_t& numberOfParticles, Vector_t<double, 3> nr, GeneratorPool &rand_pool) {

        double mu[3], sd[3];
        Vector_t<double, 3> rmin;
        Vector_t<double, 3> rmax;
        Vector_t<double, 3> hr;

        view_type &Rview = pc_m->R.getView();
        view_type &Pview = pc_m->P.getView();

        // STEP1: sample (Rx,Ry) uniformly on a unit disk.
        // Here, we use Muller-Marsaglia method https://mathworld.wolfram.com/HyperspherePointPicking.html
        // first sample normal distribution
        for(int i=0; i<3; i++){
            mu[i] = 0.0;
            sd[i] = 1.0;
        }

        MPI_Comm comm = MPI_COMM_WORLD;
        int nranks;
        int rank;
        MPI_Comm_size(comm, &nranks);
        MPI_Comm_rank(comm, &rank);

        size_type nlocal = floor(numberOfParticles/nranks);

        size_t remaining = numberOfParticles - nlocal*nranks;

        if(remaining>0 && rank==0){
            nlocal += remaining;
        }

        pc_m->create(nlocal);

        auto rand_gen = ippl::random::randn<double, 3>(Rview, rand_pool, mu, sd);
        Kokkos::parallel_for(nlocal, KOKKOS_LAMBDA(const int i) {
            rand_gen(i);
        });

        Kokkos::fence();

        // then bring (Rx,Ry) to a unit ring, and then unit disk
        Kokkos::parallel_for(
               "unitDisk", nlocal, KOKKOS_LAMBDA(const size_t j) {
                double r = Kokkos::sqrt( Kokkos::pow( Rview(j)[0], 2 ) + Kokkos::pow( Rview(j)[1], 2 ) );
                // (Rx,Ry) to a unit ring
                Rview(j)[0] /= r;
                Rview(j)[1] /= r;
                Rview(j)[2]  = 0.0;
                // (Rx,Ry) to a unit disk
                auto generator = rand_pool.get_state();
                double scale = generator.drand(0., 1.);
                Rview(j)[0] *= scale;
                Rview(j)[1] *= scale;
                rand_pool.free_state(generator);
        });
        Kokkos::fence();

       // STEP2: scale Rx and Ry with sigmaRx and sigmaRy. Also, set Px=Py=0.
        Vector_t<double, 3> sigmaR = opalDist_m->getSigmaR();
        Kokkos::parallel_for(
               "scaleRxRy", nlocal, KOKKOS_LAMBDA(const size_t j) {
               Rview(j)[0] *= sigmaR[0];
               Rview(j)[1] *= sigmaR[1];
               Pview(j)[0] = 0.0;
               Pview(j)[1] = 0.0;
        });
	Kokkos::fence();
    }

    void generateLongFlattopT(size_t& nlocal, GeneratorPool &rand_pool){

        using size_type = ippl::detail::size_type;
        view_type &Rview = pc_m->R.getView();
        view_type &Pview = pc_m->P.getView();
        auto &tview = pc_m->t.getView();

        double flattopTime = opalDist_m->getTPulseLengthFWHM()
             - std::sqrt(2.0 * std::log(2.0)) * (opalDist_m->getSigmaTRise() + opalDist_m->getSigmaTFall());

        if (flattopTime < 0.0)
              flattopTime = 0.0;

        double normalizedFlankArea = 0.5 * std::sqrt(Physics::two_pi) * gsl_sf_erf(opalDist_m->getCutoffR()[2] / std::sqrt(2.0));
        double distArea = flattopTime
                + (opalDist_m->getSigmaTRise() + opalDist_m->getSigmaTFall()) * normalizedFlankArea;

        // Find number of particles in rise, fall and flat top.
        size_type numRise = nlocal * opalDist_m->getSigmaTRise() * normalizedFlankArea / distArea;
        size_type numFall = nlocal * opalDist_m->getSigmaTFall() * normalizedFlankArea / distArea;
        size_type numFlat = nlocal - numRise - numFall;

        // Generate particles in tail.
        double sigmaTFall = opalDist_m->getSigmaTFall();
        const double par[2] = {0.0, sigmaTFall};
        using Dist_t = ippl::random::NormalDistribution<double, 1>;
        using sampling_t = ippl::random::InverseTransformSampling<double, 1, Kokkos::DefaultExecutionSpace, Dist_t>;
        Dist_t dist(par);
        Vector<double, 1> tmin = 0.0;
        Vector<double, 1> tmax = opalDist_m->getSigmaTFall() * opalDist_m->getCutoffR()[2];

        using size_type = ippl::detail::size_type;
        sampling_t sampling(dist, tmax, tmin, tmax, tmin, numFall);
        sampling.generate(tview, rand_pool);

        Vector_t<double, 3> cutoffR = opalDist_m->getCutoffR();
        Kokkos::parallel_for(
               "falltime", numFall, KOKKOS_LAMBDA(const size_t j) {
               tview(j) = -tview(j) + sigmaTFall * cutoffR[2];
        });
	Kokkos::fence();

        // Generate particles in flat top.
        double modulationAmp = opalDist_m->getFTOSCAmplitude() / 100.0;
        if (modulationAmp > 1.0)
            modulationAmp = 1.0;
        double numModulationPeriods = opalDist_m->getFTOSCPeriods();

        double modulationPeriod = 0.0;
        if (numModulationPeriods != 0.0)
            modulationPeriod = flattopTime / numModulationPeriods;

        double two_pi = Physics::two_pi;
        Kokkos::parallel_for(
            "onflattop", Kokkos::RangePolicy<>(numFall, numFall+numFlat), KOKKOS_LAMBDA(const size_t j) {

            if (modulationAmp == 0.0 || numModulationPeriods == 0.0) {
                auto generator = rand_pool.get_state();
                double r = generator.drand(0., 1.);
                tview(j) = r * flattopTime;
                rand_pool.free_state(generator);
            }
            else{
                bool allow = false;
                double randNums[2] = {0.0, 0.0};
                double temp = 0.0;
                while (!allow) {
                   auto generator = rand_pool.get_state();
                   randNums[0]= generator.drand(0., 1.);
                   randNums[1]= generator.drand(0., 1.);
                   rand_pool.free_state(generator);

                   temp = randNums[0] * flattopTime;

                   double funcValue = (1.0 + modulationAmp
                                    * Kokkos::sin(two_pi * temp / modulationPeriod))
                                    / (1.0 + Kokkos::abs(modulationAmp));

                    allow = (randNums[1] <= funcValue);
                }
                tview(j) = temp;
                tview(j) += sigmaTFall * cutoffR[2];
            }
        });
	Kokkos::fence();

        // Generate particles in rise.
        const double par2[2] = {0.0, opalDist_m->getSigmaTRise() };
        using Dist_t = ippl::random::NormalDistribution<double, 1>;
        Dist_t dist2(par2);
        tmin = 0.0;
        tmax = opalDist_m->getSigmaTRise() * opalDist_m->getCutoffR()[2];

        sampling_t sampling2(dist2, tmax, tmin, numRise);

        Kokkos::View<ippl::Vector<double, 1>*> subview(&tview(numFall + numFlat), numRise);
        sampling2.generate(subview, rand_pool);

        Kokkos::parallel_for(
               "rise", numRise, KOKKOS_LAMBDA(const size_t j) {
                tview(j+numFall + numFlat) = subview(j) + sigmaTFall * cutoffR[2] + flattopTime;
        });

        // Set Pz=Rz=0
        Kokkos::parallel_for(
               "RzPz=0", nlocal, KOKKOS_LAMBDA(const size_t j) {
               Rview(j)[2] = 0.0;
               Pview(j)[2] = 0.0;
        });

        Kokkos::parallel_for("custom_sort", Kokkos::RangePolicy<>(0, nlocal - 1), KOKKOS_LAMBDA(const int i) {
            for (size_type j = i + 1; j < nlocal; ++j) {
                if (VectorComparator()(tview(j), tview(i))) {
                  auto temp = tview(i);
                  tview(i) = tview(j);
                  tview(j) = temp;
                }
            }
        });

    }

    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override {

        size_t randInit;
        if (Options::seed == -1) {
            randInit = 1234567;
            *gmsg << "* Seed = " << randInit << " on all ranks" << endl;
        }
	else
            randInit = (size_t)(Options::seed + 100 * ippl::Comm->rank());

        GeneratorPool rand_pool(randInit);

        *gmsg << "* generate particles on a disc" << endl;
        generateUniformDisk(numberOfParticles, nr, rand_pool);

        size_type nlocal = pc_m->getLocalNum();

        *gmsg << "* sample particle time" << endl;
        generateLongFlattopT(nlocal, rand_pool);

        *gmsg << "* Done with flat top particle generation" << endl;
    }
};
#endif
