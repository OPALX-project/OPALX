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

    using Matrix_t = ippl::Vector< ippl::Vector<double, 6>, 6>;

    Vector_t<double, 3> muR_m, muP_m;
    Matrix_t cov_m;
    Matrix_t L_m;

    Vector_t<double, 3> rmin_m;
    Vector_t<double, 3> rmax_m;
    Vector_t<double, 3> pmin_m;
    Vector_t<double, 3> pmax_m;

    Vector_t<double, 3> normRmin_m;
    Vector_t<double, 3> normRmax_m;
    Vector_t<double, 3> normPmin_m;
    Vector_t<double, 3> normPmax_m;

    Vector_t<double, 6> min_m;
    Vector_t<double, 6> max_m;

    Vector_t<double, 6> normMin_m;
    Vector_t<double, 6> normMax_m;

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

        view_type &Rview = pc_m->R.getView();
        view_type &Pview = pc_m->P.getView();
        auto &Dtview = pc_m->dt.getView();

       // TODO: sample particle time according to the flat top profile
       // TODO: set Rz=Pz to zero


        double flattopTime = opalDist_m->getTPulseLengthFWHM()
             - std::sqrt(2.0 * std::log(2.0)) * (opalDist_m->getSigmaTRise() + opalDist_m->getSigmaTFall());

        if (flattopTime < 0.0)
              flattopTime = 0.0;

        double normalizedFlankArea = 0.5 * std::sqrt(Physics::two_pi) * gsl_sf_erf(opalDist_m->getCutoffR()[2] / std::sqrt(2.0));
        double distArea = flattopTime
                + (getSigmaTRise() + opalDist_m->getSigmaTFall()) * normalizedFlankArea;

        // Find number of particles in rise, fall and flat top.
        size_t numRise = numberOfParticles * opalDist_m->getSigmaTRise() * normalizedFlankArea / distArea;
        size_t numFall = numberOfParticles * opalDist_m->getSigmaTFall() * normalizedFlankArea / distArea;
        size_t numFlat = numberOfParticles - numRise - numFall;

        // Generate particles in tail.
        int saveProcessor = -1;
        const int myNode = Ippl::myNode();
        const int numNodes = Ippl::getNodes();
        const bool scalable = Attributes::getBool(itsAttr[Attrib::Distribution::SCALABLE]);
/*
        for (size_t partIndex = 0; partIndex < numFall; partIndex++) {

           double t = 0.0;
           double pz = 0.0;

          bool allow = false;
          while (!allow) {
              t = gsl_ran_gaussian_tail(randGen_m, 0, sigmaTFall_m);
              if (t <= sigmaTFall_m * cutoffR_m[2]) {
                  t = -t + sigmaTFall_m * cutoffR_m[2];
                  allow = true;
              }
          }

          // Save to each processor in turn.
          saveProcessor++;
          if (saveProcessor >= numNodes)
              saveProcessor = 0;

          if (scalable || myNode == saveProcessor) {
              tOrZDist_m.push_back(t);
          }
*/
    }

    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override {

        GeneratorPool rand_pool64((size_t)(Options::seed + 100 * ippl::Comm->rank()));

        generateUniformDisk(numberOfParticles, nr);

        generateLongFlattopT(numberOfParticles);
    }
};
#endif
