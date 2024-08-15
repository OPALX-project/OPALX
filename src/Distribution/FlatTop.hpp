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

    void generateUniformDisk(size_t& numberOfParticles) {
        GeneratorPool rand_pool64((size_t)(Options::seed + 100 * ippl::Comm->rank()));

        view_type &Rview = pc_m->R.getView();
        view_type &Pview = pc_m->P.getView();

       // TODO: sample (Rx,Ry) uniformly on a disk of radius.
       // TODO: then, scale Rx and Ry with sigmaR_m[0] and sigmaR_m[1]
       // TODO: set Px=Py=0

    }

    void generateLongFlattopT(size_t& numberOfParticles){

        view_type &Rview = pc_m->R.getView();
        view_type &Pview = pc_m->P.getView();
        auto &Dtview = pc_m->dt.getView();

       // TODO: sample particle time according to the flat top profile
       // TODO: set Rz=Pz to zero
    }

    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override {

        GeneratorPool rand_pool64((size_t)(Options::seed + 100 * ippl::Comm->rank()));

        generateUniformDisk(numberOfParticles);

        generateLongFlattopT(numberOfParticles);
    }
};
#endif
