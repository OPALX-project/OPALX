/**
 * @file FlatTop.h
 * @brief Defines the FlatTop class used for sampling emitting particles.
 */

#ifndef IPPL_FLAT_TOP_H
#define IPPL_FLAT_TOP_H

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
 * @class FlatTop
 * @brief Implements the sampling method for the flat-top distribution.
 *        x and y coordinates are uniformly distributed inside a circle
 *        and number of particles entering domain in [t, t+dt] follows flattop profile.
 *
 *        The FlatTop distribution is
 *        f(t)/Z = exp[ -((t-riseTime_m)/sigma)^2/2 ]                            t < riseTime
 *                 1.0                                                           riseTime < t < t<riseTime + flattopTime
 *                 exp[ -((t-(fallTime_m + flattopTime_m))/sigmaTFall_m)^2/2 ]   t>riseTime + flattopTime
 *         where Z is the normalizing factor.
 */
class FlatTop : public SamplingBase {
public:
    /**
     * @brief Constructor for FlatTop.
     * @param pc Shared pointer to ParticleContainer.
     * @param fc Shared pointer to FieldContainer.
     * @param opalDist Shared pointer to Distribution.
     */
    FlatTop(std::shared_ptr<ParticleContainer_t> &pc, std::shared_ptr<FieldContainer_t> &fc, std::shared_ptr<Distribution_t> &opalDist);

    /**
     * @brief Tests the number of emitted particles over a given number of steps.
     * @param nsteps Number of time steps to simulate.
     * @param dt Time step size.
     */
    void testNumEmitParticles(size_type nsteps, double dt) override;

    /**
     * @brief Tests particle emission over a given number of steps.
     * @param nsteps Number of time steps to simulate.
     * @param dt Time step size.
     */
    void testEmitParticles(size_type nsteps, double dt) override;

private:
    using size_type = ippl::detail::size_type;
    GeneratorPool rand_pool_m;  ///< Random number generator pool.
    double flattopTime_m;       ///< Time duration of when the time profile is flat.
    double normalizedFlankArea_m; ///< Normalized area of the distribution flanks.
    double distArea_m; ///< Total area of the flattop distribution.
    double sigmaTFall_m; ///< Standard deviation for fall time profile.
    double sigmaTRise_m; ///< Standard deviation for rise time profile.
    Vector_t<double, 3> cutoffR_m; ///< Cutoff radius.
    double fallTime_m; ///< Time duration for the fall phase.
    double riseTime_m; ///< Time duration for the rise phase.
    bool emitting_m; ///< Flag for particle emission status.
    size_type totalN_m; ///< Total number of particles.
    bool withDomainDecomp_m; ///< Flag for domain decomposition.
    double emissionTime_m; ///< Total emission time.
    Vector_t<double, 3> nr_m; ///< Number of grid points per direction.
    Vector_t<double, 3> hr_m; ///< Grid spacing.

    /**
     * @brief Sets whether to use domain decomposition.
     * @param withDomainDecomp Boolean flag for domain decomposition.
     */
    void setWithDomainDecomp(bool withDomainDecomp) override;

    /**
     * @brief Determines the random seed initialization.
     * @return The seed value.
     */
    static size_t determineRandInit();

    /**
     * @brief Sets distribution parameters.
     * @param opalDist Shared pointer to the distribution object.
     */
    void setParameters(const std::shared_ptr<Distribution_t> &opalDist);

public:
    /**
     * @brief Generates particles (x,y) uniformly on a disk distribution.
     * @param nlocal Number of local particles.
     * @param nNew Number of new particles to generate.
     */
    void generateUniformDisk(size_type nlocal, size_t nNew);

    /**
     * @brief Sets the number of grid points per direction.
     * @param nr Vector specifying the number of grid points.
     */
    void setNr(Vector_t<double, 3> nr);

    /**
     * @brief Generates particles with a given number and grid configuration.
     * @param numberOfParticles Number of particles to generate.
     * @param nr Number of grid points in each dimension.
     */
    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;

    /**
     * @brief Computes the flat-top profile value at a given time.
     * @param t Time value.
     * @return Profile value at time t.
     */
    double FlatTopProfile(double t);

    /**
     * @brief Computes the local number of particles uniformly distributed among ranks.
     * @param nglobal Total global number of particles.
     * @return Number of local particles.
     */
    size_t computeNlocalUniformly(size_t nglobal);

    /**
     * @brief Integrates using the trapezoidal rule.
     * @param x1 begining of interval [x1,x2].
     * @param x2 end of interval [x1,x2].
     * @param y1 value of f(x1).
     * @param y2 value of f(x2).
     * @return Integrated result.
     */
    double integrateTrapezoidal(double x1, double x2, double y1, double y2);

    /**
     * @brief Initializes the domain decomposition.
     * @param BoxIncr Box increment factor.
     */
    void initDomainDecomp(double BoxIncr) override;

    /**
     * @brief Counts the number of particles entering per rank in a given time interval.
     * @param t0 Start time.
     * @param tf End time.
     * @return Number of entering particles per rank.
     */
    double countEnteringParticlesPerRank(double t0, double tf);

    /**
     * @brief Allocates memory for a given number of particles.
     * @param numberOfParticles Number of particles to allocate.
     */
    void allocateParticles(size_t numberOfParticles);

    /**
     * @brief Emits new particles within a given time interval.
     * @param t Start time.
     * @param dt Time step.
     */
    void emitParticles(double t, double dt) override;
};

#endif // IPPL_FLAT_TOP_H

