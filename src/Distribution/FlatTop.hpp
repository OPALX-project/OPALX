#ifndef IPPL_FLAT_TOP_H
#define IPPL_FLAT_TOP_H

#include "Distribution.h"
#include "SamplingBase.hpp"
#include <memory>

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t = FieldContainer<double, 3>;
using Distribution_t = Distribution;
using GeneratorPool = typename Kokkos::Random_XorShift64_Pool<>;
using Dist_t = ippl::random::NormalDistribution<double, 3>;

class FlatTop : public SamplingBase {
public:
    FlatTop(std::shared_ptr<ParticleContainer_t> &pc, std::shared_ptr<FieldContainer_t> &fc, std::shared_ptr<Distribution_t> &opalDist)
        : SamplingBase(pc, fc, opalDist) {}

    void generateUniformDisk(size_t& numberOfParticles, Vector_t<double, 3> nr) {
        GeneratorPool rand_pool64((size_t)(Options::seed + 100 * ippl::Comm->rank()));

        double mu[3], sd[3];
        Vector_t<double, 3> rmin;
        Vector_t<double, 3> rmax;
        Vector_t<double, 3> hr;

        view_type &Rview = pc_m->R.getView();
        view_type &Pview = pc_m->P.getView();

        auto& mesh = fc_m->getMesh();
        auto& FL = fc_m->getFL();

        // STEP1: sample (Rx,Ry) uniformly on a unit disk.
        // Here, we use Muller-Marsaglia method https://mathworld.wolfram.com/HyperspherePointPicking.html
        // first sample normal distribution
        for(int i=0; i<3; i++){
            mu[i] = 0.0;
            sd[i] = 1.0;
        }
        const double par[6] = {mu[0], sd[0], mu[1], sd[1], mu[2], sd[2]};
        using Dist_t = ippl::random::NormalDistribution<double, Dim>;
        using sampling_t = ippl::random::InverseTransformSampling<double, Dim, Kokkos::DefaultExecutionSpace, Dist_t>;
        Dist_t dist(par);

        hr = 6.0 / nr;
        rmin = -3.0;
        rmax = 3.0;
        mesh.setMeshSpacing(hr);
        mesh.setOrigin(rmin);
        pc_m->getLayout().updateLayout(FL, mesh);

        ippl::detail::RegionLayout<double, Dim, Mesh_t<Dim>> rlayout;
        rlayout = ippl::detail::RegionLayout<double, Dim, Mesh_t<Dim>>(FL, mesh);

        sampling_t sampling(dist, rmax, rmin, rlayout, numberOfParticles);

        size_type nlocal = sampling.getLocalSamplesNum();
        pc_m->create(nlocal);

        sampling.generate(Rview, rand_pool64);

        // then bring (Rx,Ry) to a unit ring, and then unit disk
        Kokkos::parallel_for(
               "unitDisk", pc_m->getLocalNum(), KOKKOS_LAMBDA(const size_t j) {
                double r = Kokkos::sqrt( Kokkos::pow( Rview(j)[0], 2 ) + Kokkos::pow( Rview(j)[1], 2 ) );
                // (Rx,Ry) to a unit ring
                Rview(j)[0] /= r;
                Rview(j)[1] /= r;
                Rview(j)[2]  = 0.0;
                // (Rx,Ry) to a unit disk
                auto generator = rand_pool64.get_state();
                double scale = generator.drand(0., 1.);
                Rview(j)[0] *= scale;
                Rview(j)[1] *= scale;
                rand_pool64.free_state(generator);
        });
        Kokkos::fence();

       // STEP2: scale Rx and Ry with sigmaRx and sigmaRy. Also, set Px=Py=0.
        Vector_t<double, 3> sigmaR = opalDist_m->getSigmaR();
        Kokkos::parallel_for(
               "scaleRxRy", pc_m->getLocalNum(), KOKKOS_LAMBDA(const size_t j) {
               Rview(j)[0] *= sigmaR[0];
               Rview(j)[1] *= sigmaR[1];
               Pview(j)[0] = 0.0;
               Pview(j)[1] = 0.0;
        });
	Kokkos::fence();

    }

    void generateLongFlattopT(size_t& numberOfParticles){

        using size_type = ippl::detail::size_type;
        GeneratorPool rand_pool64((size_t)(Options::seed + 100 * ippl::Comm->rank()));
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
        size_type numRise = numberOfParticles * opalDist_m->getSigmaTRise() * normalizedFlankArea / distArea;
        size_type numFall = numberOfParticles * opalDist_m->getSigmaTFall() * normalizedFlankArea / distArea;
        size_type numFlat = numberOfParticles - numRise - numFall;

        // Generate particles in tail.
        const double par[2] = {0.0, opalDist_m->getSigmaTFall() };
        using Dist_t = ippl::random::NormalDistribution<double, 1>;
        using sampling_t = ippl::random::InverseTransformSampling<double, 1, Kokkos::DefaultExecutionSpace, Dist_t>;
        Dist_t dist(par);
        Vector<double, 1> tmin = 0.0;
        Vector<double, 1> tmax = opalDist_m->getSigmaTFall() * opalDist_m->getCutoffR()[2];

        using size_type = ippl::detail::size_type;
        sampling_t sampling(dist, tmax, tmin, numFall);
        sampling.generate(tview, rand_pool64);

        double sigmaTFall = opalDist_m->getSigmaTFall();
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
                auto generator = rand_pool64.get_state();
                double r = generator.drand(0., 1.);
                tview(j) = r * flattopTime;
                rand_pool64.free_state(generator);
            }
            else{
                bool allow = false;
                double randNums[2] = {0.0, 0.0};
                double temp = 0.0;
                while (!allow) {
                   auto generator = rand_pool64.get_state();
                   randNums[0]= generator.drand(0., 1.);
                   randNums[1]= generator.drand(0., 1.);
                   rand_pool64.free_state(generator);

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
        sampling2.generate(subview, rand_pool64);

        Kokkos::parallel_for(
               "rise", numRise, KOKKOS_LAMBDA(const size_t j) {
                tview(j+numFall + numFlat) = subview(j) + sigmaTFall * cutoffR[2] + flattopTime;
        });

        // Set Pz=Rz=0
        Kokkos::parallel_for(
               "RzPz=0", numberOfParticles, KOKKOS_LAMBDA(const size_t j) {
               Rview(j)[2] = 0.0;
               Pview(j)[2] = 0.0;
        });
    }

    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override {

        GeneratorPool rand_pool64((size_t)(Options::seed + 100 * ippl::Comm->rank()));

        generateUniformDisk(numberOfParticles, nr);

        size_type nlocal = pc_m->getLocalNum();

        generateLongFlattopT(nlocal);
    }
};
#endif
