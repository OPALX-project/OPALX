#ifndef IPPL_MULTI_VARIATE_GAUSSIAN_H
#define IPPL_MULTI_VARIATE_GAUSSIAN_H

#include "Distribution.h"
#include "SamplingBase.hpp"
#include <memory>

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t = FieldContainer<double, 3>;
using Distribution_t = Distribution;
using GeneratorPool = typename Kokkos::Random_XorShift64_Pool<>;
using Dist_t = ippl::random::NormalDistribution<double, 3>;

class MultiVariateGaussian : public SamplingBase {
public:
    MultiVariateGaussian(std::shared_ptr<ParticleContainer_t> &pc, std::shared_ptr<FieldContainer_t> &fc, std::shared_ptr<Distribution_t> &opalDist)
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

    void ComputeCholeskyFactorization(){
            for (unsigned int i = 0; i < 6; i++) {
               for (unsigned int j = 0; j < 6; j++) {
                   L_m[i][j] = 0.0;
               }
            }
            double sum = 0.0;
            for (unsigned int i = 0; i < 6; i++) {
               for (unsigned int j = 0; j <= i; j++) {
                  sum = 0.0;
                  for (unsigned int k = 0; k < j; k++) {
                      sum += L_m[i][k] * L_m[j][k];
                  }
                  if (j == i) {
                     L_m[j][j] = Kokkos::sqrt(cov_m[j][j] - sum);
                  }
                  else {
                     L_m[i][j] = (cov_m[i][j] - sum) / L_m[j][j];
                  }
               }
            }
    }

    void ComputeCenteredBounds(){
        rmin_m = -opalDist_m->getCutoffR();
        rmax_m =  opalDist_m->getCutoffR();
        pmin_m = -opalDist_m->getCutoffP();
        pmax_m =  opalDist_m->getCutoffP();

        for(int i=0; i<3; i++){
            rmin_m(i) *= opalDist_m->getSigmaR()[i];
            rmax_m(i) *= opalDist_m->getSigmaR()[i];
            pmin_m(i) *= opalDist_m->getSigmaP()[i];
            pmax_m(i) *= opalDist_m->getSigmaP()[i];

            min_m(i*2) = rmin_m(i);
            max_m(i*2) = rmax_m(i);
            min_m(i*2+1) = pmin_m(i);
            max_m(i*2+1) = pmax_m(i);
        }

        normMin_m = 0.0;
        normMax_m = 0.0;
        double sumMin, sumMax;
        for(int i=0; i<6; i++){
            sumMin = 0.0;
            sumMax = 0.0;
            for(int j=0; j<i; j++){
               sumMin += -L_m[i][j]*normMin_m(j);
               sumMax += -L_m[i][j]*normMax_m(j);
            }
            normMin_m(i) = (min_m(i)-sumMin)/L_m[i][i];
            normMax_m(i) = (max_m(i)-sumMax)/L_m[i][i];
        }

        for(int i=0; i<3; i++){
            normRmin_m(i) = min_m(2*i)/opalDist_m->getSigmaR()[i];
            normRmax_m(i) = max_m(2*i)/opalDist_m->getSigmaR()[i];
            normPmin_m(i) = min_m(2*i+1)/opalDist_m->getSigmaP()[i];
            normPmax_m(i) = max_m(2*i+1)/opalDist_m->getSigmaP()[i];

            rmin_m(i) /= opalDist_m->getSigmaR()[i];
            rmax_m(i) /= opalDist_m->getSigmaR()[i];
            pmin_m(i) /= opalDist_m->getSigmaP()[i];
            pmax_m(i) /= opalDist_m->getSigmaP()[i];
        }
    }

    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override {
        GeneratorPool rand_pool64((size_t)(Options::seed + 100 * ippl::Comm->rank()));

        // read the input covariance matrix from opal dist
        for (unsigned int i = 0; i < 6; i++) {
            for (unsigned int j = 0; j < 6; j++) {
                cov_m[i][j] = opalDist_m->correlationMatrix_m[i][j];
            }
        }

        // compute L using Cholesky factorization cov=L*LT
        ComputeCholeskyFactorization();

        // compute boundaries of random numbers
        ComputeCenteredBounds();

        Vector_t<double, 3> hr, hp;

        view_type &Rview = pc_m->R.getView();
        view_type &Pview = pc_m->P.getView();

        const double par[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
        using Dist_t = ippl::random::NormalDistribution<double, 3>;
        using sampling_t = ippl::random::InverseTransformSampling<double, 3, Kokkos::DefaultExecutionSpace, Dist_t>;
        Dist_t dist(par);

        auto& mesh = fc_m->getMesh();
        auto& FL = fc_m->getFL();

        hr = (rmax_m-rmin_m) / nr;
        mesh.setMeshSpacing(hr);
        mesh.setOrigin(rmin_m);
        pc_m->getLayout().updateLayout(FL, mesh);

        ippl::detail::RegionLayout<double, 3, Mesh_t<Dim>> rlayout;
        rlayout = ippl::detail::RegionLayout<double, 3, Mesh_t<Dim>>(FL, mesh);

        sampling_t samplingR(dist, normRmax_m, normRmin_m, rlayout, numberOfParticles);

        size_type nlocal = samplingR.getLocalSamplesNum();
        pc_m->create(nlocal);
        samplingR.generate(Rview, rand_pool64);


        auto meshP = fc_m->getMesh();
        auto FLP = fc_m->getFL();

        hp = (pmax_m-pmin_m) / nr;
        meshP.setMeshSpacing(hp);
        meshP.setOrigin(pmin_m);

        ippl::detail::RegionLayout<double, 3, Mesh_t<Dim>> playout;
        playout = ippl::detail::RegionLayout<double, 3, Mesh_t<Dim>>(FLP, meshP);

        sampling_t samplingP(dist, normPmax_m, normPmin_m, playout, numberOfParticles);
        //samplingP.setLocalSamplesNum( nlocal );
        samplingP.generate(Pview, rand_pool64);

        Matrix_t L;
        for (unsigned i = 0; i < 6; ++i){
            for (unsigned j = 0; j < 6; ++j){
                L[i][j] = L_m[i][j];
            }
        }

        Kokkos::parallel_for( nlocal,KOKKOS_LAMBDA(const int k) {
                    double vec_old[6], vec[6] = {0.0};

                    for (unsigned i = 0; i < 3; ++i){
                        vec_old[2*i] = Rview(k)[i];
                        vec_old[2*i+1] = Pview(k)[i];
                    }

                    for (unsigned i = 0; i < 6; ++i){
                        for (unsigned j = 0; j < i+1; ++j){
                            vec[i] += L[i][j]*vec_old[j];
                        }
                    }

                    for (unsigned i = 0; i < 3; ++i){
                       Rview(k)[i] = vec[2*i];
                       Pview(k)[i] = vec[2*i+1];
                   }
            }
	);
	Kokkos::fence();

        // zero mean of R
        double meanR[3], loc_meanR[3];
        for(int i=0; i<3; i++){
           meanR[i] = 0.0;
           loc_meanR[i] = 0.0;
        }

        Kokkos::parallel_reduce(
            "calc moments of particle distr.", nlocal,
            KOKKOS_LAMBDA(
                    const int k, double& cent0, double& cent1, double& cent2) {
                    cent0 += Rview(k)[0];
                    cent1 += Rview(k)[1];
                    cent2 += Rview(k)[2];
            },
            Kokkos::Sum<double>(loc_meanR[0]), Kokkos::Sum<double>(loc_meanR[1]), Kokkos::Sum<double>(loc_meanR[2]));
        Kokkos::fence();

        MPI_Allreduce(loc_meanR, meanR, 3, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
        ippl::Comm->barrier();

        for(int i=0; i<3; i++){
           meanR[i] = meanR[i]/(1.*numberOfParticles);
        }

        Kokkos::parallel_for(
                nlocal, KOKKOS_LAMBDA(
                    const int k) {
                    Rview(k)[0] -= meanR[0];
                    Rview(k)[1] -= meanR[1];
                    Rview(k)[2] -= meanR[2];
            }
        );
        Kokkos::fence();

        // zero mean of P
        double meanP[3], loc_meanP[3];
        for(int i=0; i<3; i++){
           meanP[i] = 0.0;
           loc_meanP[i] = 0.0;
        }
        Kokkos::parallel_reduce(
            "calc moments of particle distr.", nlocal,
            KOKKOS_LAMBDA(
                    const int k, double& cent0, double& cent1, double& cent2) {
                    cent0 += Pview(k)[0];
                    cent1 += Pview(k)[1];
                    cent2 += Pview(k)[2];
            },
	    Kokkos::Sum<double>(loc_meanP[0]), Kokkos::Sum<double>(loc_meanP[1]), Kokkos::Sum<double>(loc_meanP[2]));
        Kokkos::fence();

        MPI_Allreduce(loc_meanP, meanP, 3, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
        ippl::Comm->barrier();

        for(int i=0; i<3; i++){
           meanP[i] = meanP[i]/(1.*numberOfParticles);
        }

        Kokkos::parallel_for(
            nlocal, KOKKOS_LAMBDA(
                    const int k) {
                    Pview(k)[0] -= meanP[0];
                    Pview(k)[1] -= meanP[1];
                    Pview(k)[2] -= meanP[2];
            }
        );
        Kokkos::fence();

        // correct the means of R and P
        double avrgpz = opalDist_m->getAvrgpz();
        Kokkos::parallel_for(
            nlocal, KOKKOS_LAMBDA(
                    const int k) {
                    Pview(k)[2] += avrgpz;
            }
        );
        Kokkos::fence();
    }
};
#endif
