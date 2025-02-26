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

/**
 * @brief Represents a multivariate Gaussian particle sampler.
 * 
 * This class generates particles following a multivariate Gaussian distribution 
 * using Cholesky factorization and transformation sampling.
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
    MultiVariateGaussian(std::shared_ptr<ParticleContainer_t> &pc, 
                         std::shared_ptr<FieldContainer_t> &fc, 
                         std::shared_ptr<Distribution_t> &opalDist);

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

private:
    using Matrix_t = ippl::Vector<ippl::Vector<double, 6>, 6>;

    Vector_t<double, 3> muR_m, muP_m;  ///< Mean vectors for position and momentum.
    Matrix_t cov_m;                    ///< Covariance matrix.
    Matrix_t L_m;                       ///< Cholesky factorized matrix.

    Vector_t<double, 3> rmin_m, rmax_m, pmin_m, pmax_m;  ///< Sampling bounds.
    Vector_t<double, 3> normRmin_m, normRmax_m, normPmin_m, normPmax_m;
    Vector_t<double, 6> min_m, max_m, normMin_m, normMax_m;
};

#endif // IPPL_MULTI_VARIATE_GAUSSIAN_H

