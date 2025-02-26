#include "MultiVariateGaussian.h"
#include <Kokkos_Core.hpp>
#include <mpi.h>

/**
 * @brief Constructs a MultiVariateGaussian sampler.
 */
MultiVariateGaussian::MultiVariateGaussian(std::shared_ptr<ParticleContainer_t> &pc, 
                                           std::shared_ptr<FieldContainer_t> &fc, 
                                           std::shared_ptr<Distribution_t> &opalDist)
    : SamplingBase(pc, fc, opalDist) {}

/**
 * @brief Computes the Cholesky decomposition of the covariance matrix.
 */
void MultiVariateGaussian::ComputeCholeskyFactorization() {
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
            } else {
                L_m[i][j] = (cov_m[i][j] - sum) / L_m[j][j];
            }
        }
    }
}

/**
 * @brief Computes normalized boundaries for the Gaussian sampling.
 */
void MultiVariateGaussian::ComputeCenteredBounds() {
    rmin_m = -opalDist_m->getCutoffR();
    rmax_m =  opalDist_m->getCutoffR();
    pmin_m = -opalDist_m->getCutoffP();
    pmax_m =  opalDist_m->getCutoffP();

    for (int i = 0; i < 3; i++) {
        rmin_m(i) *= opalDist_m->getSigmaR()[i];
        rmax_m(i) *= opalDist_m->getSigmaR()[i];
        pmin_m(i) *= opalDist_m->getSigmaP()[i];
        pmax_m(i) *= opalDist_m->getSigmaP()[i];

        min_m(i * 2) = rmin_m(i);
        max_m(i * 2) = rmax_m(i);
        min_m(i * 2 + 1) = pmin_m(i);
        max_m(i * 2 + 1) = pmax_m(i);
    }
}

/**
 * @brief Generates particles following a multivariate Gaussian distribution.
 */
void MultiVariateGaussian::generateParticles(size_t &numberOfParticles, Vector_t<double, 3> nr) {
    extern Inform* gmsg;
    size_t randInit;
    if (Options::seed == -1) {
        randInit = 1234567;
        *gmsg << "* Seed = " << randInit << " on all ranks" << endl;
    } else {
        randInit = (size_t)(Options::seed + 100 * ippl::Comm->rank());
    }

    GeneratorPool rand_pool64(randInit);

    // Initialize covariance matrix from the distribution.
    for (unsigned int i = 0; i < 6; i++) {
        for (unsigned int j = 0; j < 6; j++) {
            cov_m[i][j] = opalDist_m->correlationMatrix_m[i][j];
        }
    }

    ComputeCholeskyFactorization();
    ComputeCenteredBounds();

    view_type &Rview = pc_m->R.getView();
    view_type &Pview = pc_m->P.getView();

    const double par[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
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

    sampling_t sampling(dist, normRmax_m, normRmin_m, normRmax_m, normRmin_m, nlocal);
    pc_m->create(nlocal);
    sampling.generate(Rview, rand_pool64);

    sampling.updateBounds(normPmax_m, normPmin_m, normPmax_m, normPmin_m);
    sampling.generate(Pview, rand_pool64);

    Matrix_t L;
    for (unsigned int i = 0; i < 6; i++) {
        for (unsigned int j = 0; j < 6; j++) {
            L[i][j] = L_m[i][j];
        }
    }

    // Apply Cholesky transformation
    Kokkos::parallel_for(nlocal, KOKKOS_LAMBDA(const int k) {
        double vec_old[6], vec[6] = {0.0};
        for (unsigned i = 0; i < 3; ++i) {
            vec_old[2 * i] = Rview(k)[i];
            vec_old[2 * i + 1] = Pview(k)[i];
        }
        for (unsigned i = 0; i < 6; ++i) {
            for (unsigned j = 0; j < i + 1; ++j) {
                vec[i] += L[i][j] * vec_old[j];
            }
        }
        for (unsigned i = 0; i < 3; ++i) {
            Rview(k)[i] = vec[2 * i];
            Pview(k)[i] = vec[2 * i + 1];
        }
    });

    Kokkos::fence();
}
