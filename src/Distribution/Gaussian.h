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
     * @brief Constructor for the Gaussian sampler without field container and distribution object.
     * @param pc Shared pointer to the particle container.
     * @param sigmaR Standard deviation for position distribution.
     * @param sigmaP Standard deviation for momentum distribution.
     * @param avrgpz Average momentum in the z-direction.
     * @param cutoffR Cutoff multiplier for position distribution.
     * @param fix_meanR Flag to exactly fix the mean position of particles after sampling.
     */
    Gaussian(std::shared_ptr<ParticleContainer_t>& pc,
            const Vector_t<double, 3>& sigmaR,
            const Vector_t<double, 3>& sigmaP,
            double avrgpz, const Vector_t<double, 3>& cutoffR,
            bool fix_meanR = true);
    /**
     * @brief Generates particles with a Gaussian distribution.
     *
     * @param numberOfParticles The total number of particles to generate.
     * @param nr the number of grid cells in R (used in domain decomposition).
     */
    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;

    void setSigmaR(const Vector_t<double, 3>& sigmaR) {
        sigmaR_m = sigmaR;
    }

    void setSigmaP(const Vector_t<double, 3>& sigmaP) {
        sigmaP_m = sigmaP;
    }

    void setAvrgpz(double avrgpz) {
        avrgpz_m = avrgpz;
    }

    void setCutoffR(const Vector_t<double, 3>& cutoffR) {
        cutoffR_m = cutoffR;
    }

    void getParameters(Vector_t<double, 3>& sigmaR,
                       Vector_t<double, 3>& sigmaP,
                       double& avrgpz,
                       Vector_t<double, 3>& cutoffR) const {
        sigmaR = sigmaR_m;
        sigmaP = sigmaP_m;
        avrgpz = avrgpz_m;
        cutoffR = cutoffR_m;
    }

    void setFixMeanR(bool fixMeanR) {
        fixMeanR_m = fixMeanR;
    }

    void getFixMeanR(bool& fixMeanR) const {
        fixMeanR = fixMeanR_m;
    }

private:
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
     * @brief Average momentum in the z-direction.
     */
    double avrgpz_m;

    /**
     * @brief Cutoff multiplier for position distribution.
     */
    Vector_t<double, 3> cutoffR_m;

    /**
     * @brief Flag to exactly fix the mean position of particles after sampling.
     */
    bool fixMeanR_m = true;
};

#endif // IPPL_GAUSSIAN_H

