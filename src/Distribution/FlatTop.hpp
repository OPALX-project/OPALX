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
        : SamplingBase(pc, fc, opalDist), rand_pool_m(determineRandInit()), setParameters(opalDist)) {}

private:
    GeneratorPool rand_pool_m;
    double flattopTime_m;
    double normalizedFlankArea_m;
    double distArea_m;
    double sigmaTFall_m;
    double sigmaTRise_m;
    Vector_t<double, 3> cutoffR_m;
    double fallTime_m;
    double flattopTime_m;
    double riseTime_m;
    bool emitting_m;

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

    static setParameters(const std::shared_ptr<Distribution_t> &opalDist) {
        emitting_m = opalDist->emitting_m;
        // time span of fall is [0, fallTime]
        sigmaTFall_m = opalDist_m->getSigmaTFall();
        cutoffR_m = opalDist_m->getCutoffR();
        fallTime_m = sigmaTFall_m * cutoffR_m[2]; // fall is [0, fallTime]

        // time of span flattop is [fallTime, fallTime+flattopTime]
        flattopTime_m = opalDist->getTPulseLengthFWHM()
            - std::sqrt(2.0 * std::log(2.0)) * (opalDist->getSigmaTRise() + opalDist->getSigmaTFall());
        if (flattopTime_m < 0.0) {
            flattopTime_m = 0.0;
        }

        // time span of rise is [fallTime+flattopTime, fallTime+flattopTime+riseTime]
        riseTime_m = opalDist_m->getSigmaTRise() * opalDist_m->getCutoffR()[2];

        normalizedFlankArea_m = 0.5 * std::sqrt(Physics::two_pi) * std::erf(opalDist_m->getCutoffR()[2] / std::sqrt(2.0));
        distArea_m = flattopTime_m
                + (opalDist_m->getSigmaTRise() + opalDist_m->getSigmaTFall()) * normalizedFlankArea_m;
    }

public:
    void generateUniformDisk(size_t& numberOfParticles) {

        GeneratorPool rand_pool = rand_pool_m;
        Vector_t<double, 3> rmin;
        Vector_t<double, 3> rmax;
        Vector_t<double, 3> hr;

        view_type &Rview = pc_m->R.getView();
        view_type &Pview = pc_m->P.getView();

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

        double pi = Physics::pi;
        Vector_t<double, 3> sigmaR = opalDist_m->getSigmaR();
        // Sample (Rx,Ry) on a unit ring, then scale with sigmaRx and sigmaRy, set Px=Py=0
        Kokkos::parallel_for(
               "unitDisk", nlocal, KOKKOS_LAMBDA(const size_t j) {
                auto generator = rand_pool.get_state();
                double r = Kokkos::sqrt( generator.drand(0., 1.) );
                double theta = 2.0 * pi * generator.drand(0., 1.);
                rand_pool.free_state(generator);

                Rview(j)[0] = r * Kokkos::cos(theta) * sigmaR[0];
                Rview(j)[1] = r * Kokkos::sin(theta) * sigmaR[1];
                Rview(j)[2]  = 0.0;
                Pview(j)[0] = 0.0;
                Pview(j)[1] = 0.0;
        });
        Kokkos::fence();
    }

    double FlatTopProfile(double t){
        double t0;
        if(t<fallTime_m){
            t0 = fallTime_m;
            return exp( -pow((t-t0)/fallTime_m,2) /2. ) / std::sqrt(Physics::two_pi) / sigmaTFall_m;
        }
        else if( t>fallTime_m && t<fallTime_m + flattopTime_m){
            return 1.
        }
        else if(t>fallTime_m + flattopTime_m && t < fallTime_m + flattopTime_m + riseTime_m){
            t0 = fallTime_m + flattopTime_m + riseTime_m;
            return exp( -pow((t-t0)/sigmaTRise_m,2)/2. ) / std::sqrt(Physics::two_pi) / sigmaTRise_m;
        }
        else
            return 0.;
    }

    double countEnteringParticles(double t0, double tf){
        double tmin, tmax;
        if(t0 < fallTime ){
            tmin = t0;
            tmax = min(tf, fallTime);
        }
        else if( tf > fallTime && t0 < fallTime + flattoptime){
            tmin = max(t0, fallTime); // in case t0<fallTime, tmin should be fallTime
            tmax = min(tf, fallTime+flattopTime); //  in case tf>fallTime+flattopTime, tmax should be fallTime+flattopTime
        }
        else if( tf > fallTime + flattopTime && t0 < fallTime + flattopTime + riseTime){
            tmin = max(t0, fallTime+flattopTime);
            tmax = min(tf, fallTime+flattopTime+riseTime);
        }
        double tArea = 0.5 * ( FlatTopProfile(tmin) + FlatTopProfile(tmax) ); // Trapezoidal rule
        return totalN_m * tArea / distArea;
    }

    void generateLongFlattopT(size_t& nlocal){

        GeneratorPool rand_pool = rand_pool_m;

        using size_type = ippl::detail::size_type;
        view_type &Rview = pc_m->R.getView();
        view_type &Pview = pc_m->P.getView();
        auto &tview = pc_m->t.getView();

        double flattopTime = opalDist_m->getTPulseLengthFWHM()
             - std::sqrt(2.0 * std::log(2.0)) * (opalDist_m->getSigmaTRise() + opalDist_m->getSigmaTFall());

        if (flattopTime < 0.0)
              flattopTime = 0.0;

        double normalizedFlankArea = 0.5 * std::sqrt(Physics::two_pi) * std::erf(opalDist_m->getCutoffR()[2] / std::sqrt(2.0));
        double distArea = flattopTime
                + (opalDist_m->getSigmaTRise() + opalDist_m->getSigmaTFall()) * normalizedFlankArea;

        // Find number of particles in rise, fall and flat top.
        size_type numRise = nlocal * opalDist_m->getSigmaTRise() * normalizedFlankArea / distArea;
        size_type numFall = nlocal * opalDist_m->getSigmaTFall() * normalizedFlankArea / distArea;
        size_type numFlat = nlocal - numRise - numFall;

        // Generate particles in tail.
        double sigmaTFall = opalDist_m->getSigmaTFall();
        Vector_t<double, 3> cutoffR = opalDist_m->getCutoffR();
        const double par[2] = {0.0, sigmaTFall};
        using Dist_t = ippl::random::NormalDistribution<double, 1>;
        using sampling_t = ippl::random::InverseTransformSampling<double, 1, Kokkos::DefaultExecutionSpace, Dist_t>;
        Dist_t dist(par);
        Vector<double, 1> tmin = 0.0;
        Vector<double, 1> tmax = sigmaTFall * cutoffR[2];

        sampling_t sampling(dist, tmax, tmin, tmax, tmin, numFall);
        sampling.generate(tview, rand_pool);

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
                rand_pool.free_state(generator);
                tview(j) = r * flattopTime;
            }
            else{
                bool allow = false;
                double randNums[2] = {0.0, 0.0};
                double temp = 0.0;
                while (!allow) {
                   auto generator = rand_pool.get_state();
                   randNums[0] = generator.drand(0., 1.);
                   randNums[1] = generator.drand(0., 1.);
                   rand_pool.free_state(generator);

                   temp = randNums[0] * flattopTime;

                   double funcValue = (1.0 + modulationAmp
                                    * Kokkos::sin(two_pi * temp / modulationPeriod))
                                    / (1.0 + Kokkos::abs(modulationAmp));

                    allow = (randNums[1] <= funcValue);
                }
                tview(j) = temp;
            }
            tview(j) += sigmaTFall * cutoffR[2];
        });
	Kokkos::fence();

        // Generate particles in rise.
        const double par2[2] = {0.0, opalDist_m->getSigmaTRise() };
        Dist_t dist2(par2);
        tmin = 0.0;
        tmax = opalDist_m->getSigmaTRise() * opalDist_m->getCutoffR()[2];

        sampling_t sampling2(dist2, tmax, tmin, tmax, tmin, numRise);

        auto subview = Kokkos::subview(tview, std::make_pair(numFall + numFlat, numFall + numFlat + numRise));
        sampling2.generate(subview, rand_pool);

        Kokkos::parallel_for(
               "rise", numRise, KOKKOS_LAMBDA(const size_t j) {
                subview(j) += + sigmaTFall * cutoffR[2] + flattopTime;
        });

        // Set Pz=Rz=0
        Kokkos::parallel_for(
               "RzPz=0", nlocal, KOKKOS_LAMBDA(const size_t j) {
               Rview(j)[2] = 0.0;
               Pview(j)[2] = 0.0;
        });

        // Assume tview is Kokkos::View<ippl::Vector<double, 1>*>
        Kokkos::View<double*> keys("keys", nlocal);

        // Extract the sortable keys (e.g., a[0])
        Kokkos::parallel_for("ExtractKeys", Kokkos::RangePolicy<>(0, nlocal), KOKKOS_LAMBDA(const int i) {
            keys(i) = tview(i)[0];
        });

        // Sort the keys, by default kokkos::sort sorts array in an ascending order
        Kokkos::sort(keys);

        // Reorder tview based on sorted keys in a descending order
        Kokkos::parallel_for("Reorder", Kokkos::RangePolicy<>(0, nlocal), KOKKOS_LAMBDA(const int i) {
            tview(i)[0] = keys(nlocal-i-1);
        });

        // sanity check
        bool sorted = true;
        Kokkos::parallel_reduce(
            "CheckSorted",
            tview.extent(0) - 1, // Check up to the second-to-last element
            KOKKOS_LAMBDA(const size_type i, bool &local_sorted) {
                if (tview(i)[0] < tview(i + 1)[0]) {
                    local_sorted = false;
                }
            },
            Kokkos::LAnd<bool>(sorted) // Logical AND reduction to ensure all checks pass
        );
        if(sorted == true)
            *gmsg << "* Sorting of particles w.r.t. emission time PAASED the sanity check!\n";
        else
            *gmsg << "* Sorting of particles w.r.t emission time FAILED the sanity check!\n";

    }

    void setEmissionTime() {
        opalDist_m->setTEmission( opalDist_m->getTPulseLengthFWHM() + (opalDist_m->getCutoffR()[2] - std::sqrt(2.0 * std::log(2.0))) * ( opalDist_m->getSigmaTRise() + opalDist_m->getSigmaTFall() ) );
    }

    void shiftBeam(size_t nlocal){

        double tEmission = opalDist_m->getTEmission();
        auto &tview = pc_m->t.getView();

        Kokkos::parallel_for("shift t", nlocal, KOKKOS_LAMBDA(const size_t j) {
               tview(j) -= tEmission;
        });

    }

    void generateParticles(size_t& numberOfParticles) override {
        if(!emitting_m){
            *gmsg << "* generate particles on a disc" << endl;
            generateUniformDisk(numberOfParticles);

            size_type nlocal = pc_m->getLocalNum();

            *gmsg << "* sample particle time" << endl;
            generateLongFlattopT(nlocal);

            *gmsg << "* shift beam to have negative time" << endl;

            setEmissionTime();
            shiftBeam(nlocal);
        }
        else{
            // TODO: somehow, we need to bring t and dt here
            n_new = countEnteringParticles(t, t + dt);

            *gmsg << "* generate particles on a disc" << endl;
            generateUniformDisk(n_new);// TODO: make sure that the ranges here are from nlocal to nlocal+n_new

        }
        *gmsg << "* Done with flat top particle generation" << endl;
/*
        auto tViewDevice  = pc_m->t.getView();
        auto tView = Kokkos::create_mirror_view(tViewDevice);
        Kokkos::deep_copy(tView,tViewDevice);

        int rank, numRanks;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
        std::string filename = "time_" + std::to_string(rank) + ".txt";
        std::ofstream file(filename);

        for (size_t i = 0; i < nlocal; ++i) {
           file << tView(i)[0] << "\n";
        }

        file.close();
*/
    }
};
#endif
