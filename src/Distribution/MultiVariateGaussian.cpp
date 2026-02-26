#include "MultiVariateGaussian.h"
#include <mpi.h>
#include <Kokkos_Core.hpp>

using Matrix_t = ippl::Vector<ippl::Vector<double, 6>, 6>;

/**
 * @brief Constructs the MultiVariateGaussian class.
 * @param pc Shared pointer to the particle container.
 * @param fc Shared pointer to the field container.
 * @param opalDist Shared pointer to the distribution.
 */
MultiVariateGaussian::MultiVariateGaussian(
    std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
    std::shared_ptr<Distribution_t> opalDist
)
    : SamplingBase(pc, fc, opalDist) {
    // Initialize covariance matrix from the distribution.
    for (unsigned int i = 0; i < 6; i++) {
        for (unsigned int j = 0; j < 6; j++) {
            cov_m[i][j] = opalDist_m->correlationMatrix_m[i][j];
        }
    }

    setSigmaR(opalDist_m->getSigmaR());
    setSigmaP(opalDist_m->getSigmaP());
    setCutoffR(opalDist_m->getCutoffR());
    setCutoffP(opalDist_m->getCutoffP());

    meanR_m    = 0.0;
    meanP_m    = 0.0;
    meanP_m[2] = opalDist_m->getAvrgpz();

    samplerTimer_m = IpplTimings::getTimer("SamplingTimer");
    initRandomPool();
}

MultiVariateGaussian::MultiVariateGaussian(
    std::shared_ptr<ParticleContainer_t> pc, const Vector_t<double, 3>& meanR,
    const Vector_t<double, 3>& meanP, const Vector_t<double, 3>& sigmaR,
    const Vector_t<double, 3>& sigmaP, const Vector_t<double, 3>& cutoffR,
    const Vector_t<double, 3>& cutoffP, bool fixMeanR, bool fixMeanP
)
    : SamplingBase(pc) {
    // Initialize covariance matrix from the distribution.
    for (unsigned int i = 0; i < 6; i++) {
        for (unsigned int j = 0; j < 6; j++) {
            cov_m[i][j] = 0.0;
            if (i == j && i % 2 == 0) {
                cov_m[i][j] = sigmaR[i / 2] * sigmaR[i / 2];
            }
            if (i == j && i % 2 == 1) {
                cov_m[i][j] = sigmaP[i / 2] * sigmaP[i / 2];
            }
        }
    }
    setMeanR(meanR);
    setMeanP(meanP);
    setSigmaR(sigmaR);
    setSigmaP(sigmaP);
    setCutoffR(cutoffR);
    setCutoffP(cutoffP);
    setFixMeanR(fixMeanR);
    setFixMeanP(fixMeanP);

    samplerTimer_m = IpplTimings::getTimer("SamplingTimer");
    initRandomPool();
}

MultiVariateGaussian::MultiVariateGaussian(
    std::shared_ptr<ParticleContainer_t> pc, const Vector_t<double, 3>& meanR,
    const Vector_t<double, 3>& meanP, const Matrix_t& cov, const Vector_t<double, 3>& cutoffR,
    const Vector_t<double, 3>& cutoffP, bool fixMeanR, bool fixMeanP
)
    : SamplingBase(pc) {
    cov_m = cov;

    setMeanR(meanR);
    setMeanP(meanP);
    setSigmaR(
        ippl::Vector<double, 3>(
            Kokkos::sqrt(cov_m[0][0]), Kokkos::sqrt(cov_m[2][2]), Kokkos::sqrt(cov_m[4][4])
        )
    );

    setSigmaP(
        ippl::Vector<double, 3>(
            Kokkos::sqrt(cov_m[1][1]), Kokkos::sqrt(cov_m[3][3]), Kokkos::sqrt(cov_m[5][5])
        )
    );
    setCutoffR(cutoffR);
    setCutoffP(cutoffP);
    setFixMeanR(fixMeanR);
    setFixMeanP(fixMeanP);

    samplerTimer_m = IpplTimings::getTimer("SamplingTimer");
    initRandomPool();
}

/**
 * @brief Initializes the random number generator pool.
 */
void MultiVariateGaussian::initRandomPool() {
    extern Inform* gmsg;
    size_t randInit;

    if (Options::seed == -1) {
        randInit = 1234567;
        *gmsg << "* Seed = " << randInit << " on all ranks" << endl;
    } else {
        randInit = static_cast<size_t>(Options::seed + 100 * ippl::Comm->rank());
    }

    GeneratorPool rand_pool64(randInit);
    randPool_m = rand_pool64;
    return;
}

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
 * @brief Computes normalized boundaries for the multivariate Gaussian sampling.
 */
void MultiVariateGaussian::ComputeCenteredBounds() {
    rmin_m = -cutoffR_m;
    rmax_m = cutoffR_m;
    pmin_m = -cutoffP_m;
    pmax_m = cutoffP_m;

    for (int i = 0; i < 3; i++) {
        rmin_m(i) *= sigmaR_m(i);
        rmax_m(i) *= sigmaR_m(i);
        pmin_m(i) *= sigmaP_m(i);
        pmax_m(i) *= sigmaP_m(i);

        min_m(i * 2)     = rmin_m(i);
        max_m(i * 2)     = rmax_m(i);
        min_m(i * 2 + 1) = pmin_m(i);
        max_m(i * 2 + 1) = pmax_m(i);
    }

    normMin_m = 0.0;
    normMax_m = 0.0;
    double sumMin, sumMax;
    for (int i = 0; i < 6; i++) {
        sumMin = 0.0;
        sumMax = 0.0;
        for (int j = 0; j < i; j++) {
            sumMin += -L_m[i][j] * normMin_m(j);
            sumMax += -L_m[i][j] * normMax_m(j);
        }
        normMin_m(i) = (min_m(i) - sumMin) / L_m[i][i];
        normMax_m(i) = (max_m(i) - sumMax) / L_m[i][i];
    }

    for (int i = 0; i < 3; i++) {
        normRmin_m(i) = min_m(2 * i) / sigmaR_m(i);
        normRmax_m(i) = max_m(2 * i) / sigmaR_m(i);
        normPmin_m(i) = min_m(2 * i + 1) / sigmaP_m(i);
        normPmax_m(i) = max_m(2 * i + 1) / sigmaP_m(i);

        rmin_m(i) /= sigmaR_m(i);
        rmax_m(i) /= sigmaR_m(i);
        pmin_m(i) /= sigmaP_m(i);
        pmax_m(i) /= sigmaP_m(i);
    }
}

/**
 * @brief Generates particles following a multivariate Gaussian distribution.
 */
void MultiVariateGaussian::
    generateParticles(size_t& numberOfParticles, Vector_t<double, 3> /*nr*/) {
    IpplTimings::startTimer(samplerTimer_m);

    auto rand_pool64 = randPool_m;
    // compute L using Cholesky factorization cov=L*LT
    ComputeCholeskyFactorization();

    // compute boundaries of normal random numbers
    ComputeCenteredBounds();

    view_type& Rview = pc_m->R.getView();
    view_type& Pview = pc_m->P.getView();

    const double par[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    using Dist_t        = ippl::random::NormalDistribution<double, 3>;
    using sampling_t =
        ippl::random::InverseTransformSampling<double, 3, Kokkos::DefaultExecutionSpace, Dist_t>;
    Dist_t dist(par);

    MPI_Comm comm = MPI_COMM_WORLD;
    int nranks, rank;
    MPI_Comm_size(comm, &nranks);
    MPI_Comm_rank(comm, &rank);

    // if nlocal*nranks > numberOfParticles, put the remaining in rank 0
    size_t nlocal    = floor(numberOfParticles / nranks);
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
    Kokkos::parallel_for(
        nlocal, KOKKOS_LAMBDA(const int k) {
            double vec_old[6], vec[6] = {0.0};
            for (unsigned i = 0; i < 3; ++i) {
                vec_old[2 * i]     = Rview(k)[i];
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
        }
    );

    Kokkos::fence();

    // zero mean of R
    double meanR[3], loc_meanR[3];

    if (fixMeanR_m) {
        for (int i = 0; i < 3; i++) {
            meanR[i]     = 0.0;
            loc_meanR[i] = 0.0;
        }

        Kokkos::parallel_reduce(
            "calc moments of particle distr.", nlocal,
            KOKKOS_LAMBDA(const int k, double& cent0, double& cent1, double& cent2) {
                cent0 += Rview(k)[0];
                cent1 += Rview(k)[1];
                cent2 += Rview(k)[2];
            },
            Kokkos::Sum<double>(loc_meanR[0]), Kokkos::Sum<double>(loc_meanR[1]),
            Kokkos::Sum<double>(loc_meanR[2])
        );
        Kokkos::fence();

        MPI_Allreduce(loc_meanR, meanR, 3, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
        ippl::Comm->barrier();

        for (int i = 0; i < 3; i++) {
            meanR[i] = meanR[i] / (1. * numberOfParticles);
        }

        Kokkos::parallel_for(
            nlocal, KOKKOS_LAMBDA(const int k) {
                Rview(k)[0] -= meanR[0];
                Rview(k)[1] -= meanR[1];
                Rview(k)[2] -= meanR[2];
            }
        );
        Kokkos::fence();
    }

    // zero mean of P
    double meanP[3], loc_meanP[3];
    if (fixMeanP_m) {
        for (int i = 0; i < 3; i++) {
            meanP[i]     = 0.0;
            loc_meanP[i] = 0.0;
        }
        Kokkos::parallel_reduce(
            "calc moments of particle distr.", nlocal,
            KOKKOS_LAMBDA(const int k, double& cent0, double& cent1, double& cent2) {
                cent0 += Pview(k)[0];
                cent1 += Pview(k)[1];
                cent2 += Pview(k)[2];
            },
            Kokkos::Sum<double>(loc_meanP[0]), Kokkos::Sum<double>(loc_meanP[1]),
            Kokkos::Sum<double>(loc_meanP[2])
        );
        Kokkos::fence();

        MPI_Allreduce(loc_meanP, meanP, 3, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
        ippl::Comm->barrier();

        for (int i = 0; i < 3; i++) {
            meanP[i] = meanP[i] / (1. * numberOfParticles);
        }

        Kokkos::parallel_for(
            nlocal, KOKKOS_LAMBDA(const int k) {
                Pview(k)[0] -= meanP[0];
                Pview(k)[1] -= meanP[1];
                Pview(k)[2] -= meanP[2];
            }
        );
        Kokkos::fence();
    }

    // correct the means of R and P from input
    for (int i = 0; i < 3; i++) {
        meanR[i] = meanR_m[i];
        meanP[i] = meanP_m[i];
    }

    Kokkos::parallel_for(
        nlocal, KOKKOS_LAMBDA(const int k) {
            for (int i = 0; i < 3; i++) {
                Rview(k)[i] += meanR[i];
                Pview(k)[i] += meanP[i];
            }
        }
    );
    Kokkos::fence();

    IpplTimings::stopTimer(samplerTimer_m);
}
