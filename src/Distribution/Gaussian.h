#ifndef IPPL_GAUSSIAN_H
#define IPPL_GAUSSIAN_H

#include <Kokkos_Random.hpp>
#include <cmath>
#include <memory>
#include "Distribution.h"
#include "Ippl.h"
#include "OPALTypes.h"
#include "SamplingBase.hpp"
#include "Utilities/Options.h"

using GeneratorPool = typename Kokkos::Random_XorShift64_Pool<>;
using Dist_t        = ippl::random::NormalDistribution<double, 3>;

/**
 * @class Gaussian Distribution
 * @brief Generating particles following a Gaussian distribution.
 *
 * This sampler reproduces the phase-space generation of the GAUSS distribution
 * of (old) OPAL, see Distribution::generateGaussZ and Distribution::adjustPhaseSpace
 * there. Particle positions R and momenta P are drawn from
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
 * The cutoffs are applied collectively by rejection sampling on the
 * normalised (unit-sigma) coordinates, exactly as in OPAL:
 *   - transverse position:  (x / CutoffR[0])^2 + (y / CutoffR[1])^2 <= 1 (ellipse)
 *   - transverse momentum:  (px / CutoffP[0])^2 + (py / CutoffP[1])^2 <= 1 (ellipse)
 *   - longitudinal:         |z| <= CutoffR[2] and |pz| <= CutoffP[2]
 * Note that an independent per-axis truncation (a box) yields different
 * rms values and emittances than the elliptical cut (see issue #165).
 * All cutoffs must be positive; a very large value (e.g. the maximum double,
 * used by Distribution for a zero input cutoff, as in OPAL) disables the cut.
 *
 * If fixMeanR is set (default), the sampled phase space is adjusted as in
 * OPAL's adjustPhaseSpace: the means of x, y, z, px and py are shifted to
 * exactly zero, and pz is shifted such that the mean of
 *   pz + (px^2 + py^2) / (2 pz)
 * (a Taylor expansion of |p| for small transverse momenta) equals avrgpz,
 * i.e. the average beam energy is matched exactly.
 *
 * A zero spatial sigma (in metres) represents an exact point mass on that axis;
 * the corresponding coordinate is identically zero. The same holds for a zero
 * momentum sigma.
 *
 * Emission-source offsets R0/P0 are applied last and are not affected by the
 * mean/energy correction.
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
     * @param opalDist Borrowed distribution object.
     */
    Gaussian(std::shared_ptr<ParticleContainer_t> pc, Distribution_t* opalDist);

    /**
     * @brief Constructor for the Gaussian sampler without a distribution object.
     * @param pc Shared pointer to the particle container.
     * @param sigmaR Standard deviation for position distribution.
     * @param sigmaP Standard deviation for momentum distribution.
     * @param avrgpz Average momentum in the z-direction.
     * @param cutoffR Cutoff multiplier for position distribution (units of sigmaR).
     * @param fix_meanR Flag to exactly fix the mean position/momentum and the average
     *        energy of the particles after sampling (OPAL's adjustPhaseSpace).
     * @param cutoffP Cutoff multiplier for momentum distribution (units of sigmaP).
     */
    Gaussian(
            std::shared_ptr<ParticleContainer_t> pc, const Vector_t<double, 3>& sigmaR,
            const Vector_t<double, 3>& sigmaP, double avrgpz, const Vector_t<double, 3>& cutoffR,
            bool fix_meanR = true, const Vector_t<double, 3>& cutoffP = Vector_t<double, 3>(3.0));
    /**
     * @brief Generates particles with a Gaussian distribution.
     *
     * @param numberOfParticles The total number of particles to generate.
     * @param nr the number of grid cells in R (used in domain decomposition).
     */
    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;

    /**
     * @brief Time-stepped emission hook for one-shot delayed injection.
     *
     * For GAUSS, this emits once when [t, t+dt] crosses t0 for sources with t0 > 0.
     * Initial sampling at t0 == 0 is handled via generateParticles.
     */
    void emitParticles(double t, double dt) override;

    void setSigmaR(const Vector_t<double, 3>& sigmaR) { sigmaR_m = sigmaR; }

    void setSigmaP(const Vector_t<double, 3>& sigmaP) { sigmaP_m = sigmaP; }

    void setAvrgpz(double avrgpz) { avrgpz_m = avrgpz; }

    void setCutoffR(const Vector_t<double, 3>& cutoffR) { cutoffR_m = cutoffR; }

    void setCutoffP(const Vector_t<double, 3>& cutoffP) { cutoffP_m = cutoffP; }

    void getParameters(
            Vector_t<double, 3>& sigmaR, Vector_t<double, 3>& sigmaP, double& avrgpz,
            Vector_t<double, 3>& cutoffR) const {
        sigmaR  = sigmaR_m;
        sigmaP  = sigmaP_m;
        avrgpz  = avrgpz_m;
        cutoffR = cutoffR_m;
    }

    Vector_t<double, 3> getCutoffP() const { return cutoffP_m; }

    void setFixMeanR(bool fixMeanR) { fixMeanR_m = fixMeanR; }

    void getFixMeanR(bool& fixMeanR) const { fixMeanR = fixMeanR_m; }

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
     * @brief Cutoff multipliers for position and momentum distributions (units of sigma).
     */
    Vector_t<double, 3> cutoffR_m;
    Vector_t<double, 3> cutoffP_m = 3.0;

    /**
     * @brief Flag to exactly fix the mean position/momentum and the average energy
     *        of particles after sampling (equivalent of OPAL's adjustPhaseSpace).
     */
    bool fixMeanR_m = true;
};

#endif  // IPPL_GAUSSIAN_H
