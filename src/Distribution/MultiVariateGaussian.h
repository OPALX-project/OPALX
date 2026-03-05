#ifndef IPPL_MULTI_VARIATE_GAUSSIAN_H
#define IPPL_MULTI_VARIATE_GAUSSIAN_H

#include "Distribution.h"
#include "SamplingBase.hpp"
#include <Kokkos_Random.hpp>
#include "Ippl.h"
#include "Utilities/Options.h"
#include "OPALTypes.h"
#include <memory>
#include <cmath>

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t = FieldContainer<double, 3>;
using Distribution_t = Distribution;
using GeneratorPool = typename Kokkos::Random_XorShift64_Pool<>;
using Dist_t = ippl::random::NormalDistribution<double, 3>;
using Matrix_t = ippl::Vector<ippl::Vector<double, 6>, 6>;

/**
 * @class MultiVariateGaussian
 * @brief A particle generation method following multivariate Gaussian distribution.
 *
 * This class generates particles following a multivariate Gaussian distribution
 * using Cholesky factorization and inverse transformation sampling.
 *
 * Given covariance matrix cov_m = [ Cov(R0,R0), Cov(R0,P0), Cov(R0,R1), Cov(R0,P1), ...]
 * whose values are read from opalDist_m->correlationMatrix_m.
 * First, the Cholesky factorization is computed cov_m = L_m * L_m^T
 * Then, normally distribution particles R=P~N(0,I) are transformed to multivariate using L_m.
 */
class MultiVariateGaussian : public SamplingBase {
public:
    /**
     * @brief Constructor for MultiVariateGaussian.
     * 
     * @param pc Shared pointer to the particle container.
     * @param fc Shared pointer to the field container.
     * @param opalDist Shared pointer to the distribution.
     */
    MultiVariateGaussian(std::shared_ptr<ParticleContainer_t> pc, 
                         std::shared_ptr<FieldContainer_t> fc, 
                         std::shared_ptr<Distribution_t> opalDist);

    /**
     * @brief Constructor for MultiVariateGaussian with specified parameters.
     * @param meanR Mean position vector.
     * @param meanP Mean momentum vector.
     * @param sigmaR Standard deviation for position distribution.
     * @param sigmaP Standard deviation for momentum distribution.
     * @param cutoffR Cutoff position vector.
     * @param cutoffP Cutoff momentum vector.
     * @param fixMeanR Boolean flag to fix mean position.
     * @param fixMeanP Boolean flag to fix mean momentum.
     */
    MultiVariateGaussian(std::shared_ptr<ParticleContainer_t> pc,
                   const Vector_t<double, 3>& meanR,
                   const Vector_t<double, 3>& meanP,
                   const Vector_t<double, 3>& sigmaR,
                   const Vector_t<double, 3>& sigmaP,
                   const Vector_t<double, 3>& cutoffR,
                   const Vector_t<double, 3>& cutoffP,
                   bool fixMeanR=true,
                   bool fixMeanP=true);

    /**
     * @brief Constructs the MultiVariateGaussian class.
     * @param pc Shared pointer to the particle container.
     * @param meanR Mean position vector.
     * @param meanP Mean momentum vector.
     * @param cov Covariance matrix.
     * @param cutoffR Cutoff position vector.
     * @param cutoffP Cutoff momentum vector.
     * @param fixMeanR Boolean flag to fix mean position.
     * @param fixMeanP Boolean flag to fix mean momentum.
     */
    MultiVariateGaussian(std::shared_ptr<ParticleContainer_t> pc,
                   const Vector_t<double, 3>& meanR,
                   const Vector_t<double, 3>& meanP,
                   const Matrix_t &cov,
                   const Vector_t<double, 3>& cutoffR,
                   const Vector_t<double, 3>& cutoffP,
                   bool fixMeanR=true,
                   bool fixMeanP=true);
    /**
     * @brief Computes the Cholesky factorization of the covariance matrix.
     */
    void ComputeCholeskyFactorization();

    /**
     * @brief Computes centered bounds for the particle distribution.
     */
    void ComputeCenteredBounds();

    /**
     * @brief Generates particles based on the defined Gaussian distribution.
     * 
     * @param numberOfParticles Number of particles to generate.
     * @param nr Vector specifying additional sampling parameters.
     */
    void generateParticles(size_t &numberOfParticles, Vector_t<double, 3> nr) override;

    /**
     * @brief Timer for performance profiling.
     */
    IpplTimings::TimerRef samplerTimer_m;

    void setMeanR(const Vector_t<double, 3>& meanR) {
        meanR_m = meanR;
    }

    void setMeanP(const Vector_t<double, 3>& meanP) {
        meanP_m = meanP;
    }

    void setCutoffR(const Vector_t<double, 3>& cutoffR) {
        cutoffR_m = cutoffR;
    }

    void setCutoffP(const Vector_t<double, 3>& cutoffP) {
        cutoffP_m = cutoffP;
    }

    void setFixMeanR(bool fixMeanR) {
        fixMeanR_m = fixMeanR;
    }

    void setFixMeanP(bool fixMeanP) {
        fixMeanP_m = fixMeanP;
    }

    void setSigmaR(const Vector_t<double, 3>& sigmaR) {
        sigmaR_m = sigmaR;
    }
    void setSigmaP(const Vector_t<double, 3>& sigmaP) {
        sigmaP_m = sigmaP;
    }

    void setCovarianceMatrix(const Matrix_t &cov) {
        cov_m = cov;
    }

    void setL(const Matrix_t &L) {
        L_m = L;
    }

private:
    /*
    * @brief Mean vectors for position and momentum.
    */
    Vector_t<double, 3> meanR_m, meanP_m;

    /*
    * @brief Covariance matrix for the Gaussian distribution.
    */
    Matrix_t cov_m;

    /*
    * @brief Cholesky cov=LT*L factorized matrix for the Gaussian distribution.
    */
    Matrix_t L_m;

    /*
    * @brief Sampling bounds for the particle distribution.
    */
    Vector_t<double, 3> rmin_m, rmax_m, pmin_m, pmax_m;

    /*
    * @brief Normalized sampling bounds for the particle distribution.
    */
    Vector_t<double, 3> normRmin_m, normRmax_m, normPmin_m, normPmax_m;

    /**
     * @brief Min and Max bounds for all 6 dimensions (R0,P0,R1,P1,R2,P2).
     */
    Vector_t<double, 6> min_m, max_m, normMin_m, normMax_m;

    /**
     * @brief Cutoff multipliers for position and momentum distributions.
     */
    Vector_t<double, 3> cutoffR_m;
    Vector_t<double, 3> cutoffP_m;

    /**
     * @brief Initializes the random number generator pool.
     */
    void initRandomPool();

    /**
     * @brief Pool of random number generators for parallel sampling.
     */
    GeneratorPool randPool_m;

    /**
     * @brief Standard deviations for position and momentum distributions.
     */
    Vector_t<double, 3> sigmaR_m;
    Vector_t<double, 3> sigmaP_m;

    /**
     * @brief Flag to exactly fix the mean R and P of particles after sampling.
     */
    bool fixMeanR_m;
    bool fixMeanP_m;
};

#endif // IPPL_MULTI_VARIATE_GAUSSIAN_H

