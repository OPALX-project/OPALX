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
 * @class Gaussian Distribution
 * @brief Generating particles following a Gaussian distribution.
 *
 * This function samples particle positions R and momenta P given following normal distributions
 * a normal distribution
 * R ~ N ( [0, 0, 0     ], sdR)
 * P ~ N ( [0, 0, avrgpz], sdP)
 * where sdR = [SigmaR[0] 0         0
                0         SigmaR[1] 0
                0         0         SigmaR[2]]
 *
 * and sdP =    [SigmaP[0] 0         0
                0          SigmaP[1] 0
                0          0         SigmaP[2]].
 *
 * Here, R is sampled in a bounded domains R \in [-CutoffR*SigmaR, CutoffR*SigmaR]^3
 * and corrected by translation to ensure mean = [0,0,0].
 *
 * @param numberOfParticles The total number of particles to generate.
 * @param nr The number of grid points in each dimension (not used here).
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
     * @param nr the number of grid cells in R (used in domain decomposition).
     */
    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;
};

#endif // IPPL_GAUSSIAN_H

