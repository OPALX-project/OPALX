#ifndef IPPL_GAUSSIAN_H
#define IPPL_GAUSSIAN_H

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
 * @class Gaussian
 * @brief A Gaussian-based particle sampler.
 *
 * This class generates particles following a Gaussian distribution
 * in both position and momentum space.
 */
class Gaussian : public SamplingBase {
public:
    /**
     * @brief Timer for performance profiling.
     */
    IpplTimings::TimerRef samperTimer_m;

    /**
     * @brief Constructor for the Gaussian sampler.
     * 
     * @param pc Shared pointer to the particle container.
     * @param fc Shared pointer to the field container.
     * @param opalDist Shared pointer to the distribution object.
     */
    Gaussian(std::shared_ptr<ParticleContainer_t> &pc, 
             std::shared_ptr<FieldContainer_t> &fc, 
             std::shared_ptr<Distribution_t> &opalDist);

    /**
     * @brief Generates particles with a Gaussian distribution.
     *
     * @param numberOfParticles The total number of particles to generate.
     * @param nr The number of grid cells in R (used in domain decomposition).
     */
    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;
};

#endif // IPPL_GAUSSIAN_H

