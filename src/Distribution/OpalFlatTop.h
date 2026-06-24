/**
 * @file OpalFlatTop.h
 * @brief OPAL-like flat-top emission sampler with precomputed birth times.
 */

#ifndef IPPL_OPAL_FLAT_TOP_H
#define IPPL_OPAL_FLAT_TOP_H

#include <Kokkos_Random.hpp>

#include <algorithm>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "Distribution.h"
#include "Ippl.h"
#include "OPALTypes.h"
#include "SamplingBase.hpp"
#include "Utilities/Options.h"

using OpalFlatTopGeneratorPool = typename Kokkos::Random_XorShift64_Pool<>;

/**
 * @class OpalFlatTop
 * @brief Implements an old-OPAL-compatible flat-top emitter with precomputed birth times.
 *
 * x and y coordinates are sampled uniformly on an elliptical disk. The number
 * of particles emitted in each time interval follows the flat-top time profile:
 *
 * f(t)/Z = exp(-((t - riseTime_m) / sigmaTRise_m)^2 / 2),
 *          t < riseTime
 *          1.0,
 *          riseTime <= t <= riseTime + flattopTime
 *          exp(-((t - (fallTime_m + flattopTime_m)) / sigmaTFall_m)^2 / 2),
 *          t > riseTime + flattopTime
 *
 * where Z is the normalizing factor. Unlike FlatTop, this sampler builds a
 * sorted global inventory of particle birth times before tracking starts and
 * emits the matching slice during each tracker step.
 */
class OpalFlatTop : public SamplingBase {
public:
    /**
     * @brief Constructor for OpalFlatTop.
     *
     * @param pc Shared pointer to ParticleContainer.
     * @param fc Shared pointer to FieldContainer.
     * @param opalDist Borrowed Distribution.
     */
    OpalFlatTop(
            std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
            Distribution_t* opalDist);

    /**
     * @brief Constructor for OpalFlatTop with parameters passed directly.
     *
     * This is primarily intended for tests or callers that do not need a full
     * Distribution object.
     *
     * @param pc Shared pointer to ParticleContainer.
     * @param fc Shared pointer to FieldContainer.
     * @param emitting Flag indicating whether the distribution is emitting.
     * @param sigmaTFall Standard deviation of the falling flank in the time profile.
     * @param sigmaTRise Standard deviation of the rising flank in the time profile.
     * @param cutoff Cutoff multiplier for the spatial and longitudinal distributions.
     * @param tPulseLengthFWHM Time length of the pulse in FWHM.
     * @param sigmaR Semi-axis lengths for the transverse emission disk.
     * @param ftOscAmplitude Optional flat-top oscillation amplitude in percent.
     * @param ftOscPeriods Optional number of oscillation periods across the flat top.
     */
    OpalFlatTop(
            std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
            bool emitting, double sigmaTFall, double sigmaTRise, Vector_t<double, 3> cutoff,
            double tPulseLengthFWHM, Vector_t<double, 3> sigmaR, double ftOscAmplitude = 0.0,
            double ftOscPeriods = 0.0);

    /**
     * @brief Builds the global birth-time inventory for the requested particles.
     *
     * @param numberOfParticles Number of particles to generate.
     * @param nr Number of grid points in each direction.
     */
    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;

    /**
     * @brief Emits all particles whose birth times fall into the current step.
     *
     * @param t Start time of the current tracker step.
     * @param dt Time step size.
     */
    void emitParticles(double t, double dt) override;

    /**
     * @brief Returns whether all precomputed particles have been emitted.
     *
     * @return True if the inventory exists and the next global index is past its end.
     */
    bool isEmissionDone(double /*t*/) const override {
        return inventoryBuilt_m && nextGlobalIndex_m >= birthTimes_m.size();
    }

    /**
     * @brief Returns the global time shift needed to center old-OPAL pulse times.
     *
     * @return Non-negative shift from source time to the midpoint of the emission window.
     */
    double getGlobalTimeShift() const override {
        return std::max(0.0, 0.5 * emissionTime_m - t0_m);
    }

    /**
     * @brief Reports that this sampler provides an initial reference momentum.
     */
    bool hasInitialReferenceMomentum() const override { return true; }

    /**
     * @brief Returns the initial reference momentum used by the tracker.
     *
     * @return P0 for EMISSIONMODEL NONE, or the ASTRA half-sphere reference momentum.
     */
    Vector_t<double, 3> getInitialReferenceMomentum() const override;

    /**
     * @brief Returns the preferred emission time step.
     *
     * @return emissionTime divided by the configured number of emission steps.
     */
    double getEmissionTimeStep() const override {
        return emissionSteps_m > 0 ? emissionTime_m / static_cast<double>(emissionSteps_m) : 0.0;
    }

    /**
     * @brief Initializes a domain decomposition large enough for the emitted bunch.
     *
     * @param BoxIncr Percentage by which to enlarge the computed mesh box.
     */
    void initDomainDecomp(double BoxIncr) override;

    /**
     * @brief Sets whether to use domain decomposition during emission.
     *
     * @param withDomainDecomp Boolean flag for domain decomposition.
     */
    void setWithDomainDecomp(bool withDomainDecomp) override;

    /**
     * @brief Returns the total emission time.
     */
    double getEmissionTime() const { return emissionTime_m; }

    /**
     * @brief Returns the total area of the normalized flat-top distribution.
     */
    double getDistArea() const { return distArea_m; }

    /**
     * @brief Returns the number of precomputed birth times.
     */
    size_t getInventorySize() const { return birthTimes_m.size(); }

    /**
     * @brief Returns the next global inventory index to be emitted.
     */
    size_t getNextGlobalIndex() const { return nextGlobalIndex_m; }

    /**
     * @brief Returns the sorted global birth-time inventory.
     */
    const std::vector<double>& getBirthTimes() const { return birthTimes_m; }

    /**
     * @brief Replaces the birth-time inventory for tests.
     *
     * @param birthTimes Birth times to sort and install as the active inventory.
     */
    void setBirthTimesForTest(std::vector<double> birthTimes);

private:
    using size_type = ippl::detail::size_type;

    OpalFlatTopGeneratorPool rand_pool_m;  ///< Device random number generator pool.
    std::mt19937_64 host_rng_m;            ///< Host generator used to sample birth times.

    double flattopTime_m         = 0.0;  ///< Time duration of the flat profile section.
    double normalizedFlankArea_m = 0.0;  ///< Area of one normalized truncated Gaussian flank.
    double distArea_m            = 0.0;  ///< Total area of the flat-top distribution.
    double sigmaTFall_m          = 0.0;  ///< Standard deviation of the falling flank.
    double sigmaTRise_m          = 0.0;  ///< Standard deviation of the rising flank.
    Vector_t<double, 3> cutoffR_m;       ///< Cutoff multipliers for distribution support.
    double fallTime_m       = 0.0;       ///< Duration represented by the falling flank.
    double riseTime_m       = 0.0;       ///< Duration represented by the rising flank.
    bool emitting_m         = false;     ///< Flag for particle emission status.
    size_t totalN_m         = 0;         ///< Total number of particles in the inventory.
    bool withDomainDecomp_m = false;     ///< Flag for domain decomposition.
    double emissionTime_m   = 0.0;       ///< Total emission time.
    Vector_t<double, 3> nr_m;            ///< Number of grid points per direction.
    Vector_t<double, 3> hr_m;            ///< Grid spacing.
    Vector_t<double, 3> sigmaR_m;        ///< Semi-axis lengths of the transverse disk.
    double ftOscAmplitude_m = 0.0;       ///< Flat-top oscillation amplitude in percent.
    double ftOscPeriods_m   = 0.0;       ///< Number of oscillation periods across the flat top.
    size_t emissionSteps_m  = 100;       ///< Number of steps used to derive emission dt.

    std::vector<double> birthTimes_m;  ///< Sorted global particle birth times.
    size_t nextGlobalIndex_m = 0;      ///< First global birth-time index not emitted yet.
    bool inventoryBuilt_m    = false;  ///< True once birthTimes_m is ready for emission.

    /**
     * @brief Determines the device random seed initialization.
     *
     * @return Seed value for the per-rank Kokkos random pool.
     */
    static size_t determineRandInit();

    /**
     * @brief Determines the host random seed initialization.
     *
     * @return Seed value for host-side birth-time sampling.
     */
    static size_t determineHostSeed();

    /**
     * @brief Reads distribution parameters from the OPALX Distribution object.
     *
     * @param opalDist Borrowed Distribution object.
     */
    void setParameters(Distribution_t* opalDist);

    /**
     * @brief Sets derived flat-top profile and emission parameters.
     *
     * @param emitting Flag indicating whether the distribution is emitting.
     * @param sigmaTFall Standard deviation of the falling flank.
     * @param sigmaTRise Standard deviation of the rising flank.
     * @param cutoff Cutoff multiplier for distribution support.
     * @param tPulseLengthFWHM Time length of the pulse in FWHM.
     * @param sigmaR Semi-axis lengths for the transverse emission disk.
     * @param ftOscAmplitude Flat-top oscillation amplitude in percent.
     * @param ftOscPeriods Number of oscillation periods across the flat top.
     */
    void setInternalVariables(
            bool emitting, double sigmaTFall, double sigmaTRise, Vector_t<double, 3> cutoff,
            double tPulseLengthFWHM, Vector_t<double, 3> sigmaR, double ftOscAmplitude,
            double ftOscPeriods);

    /**
     * @brief Builds and sorts the global particle birth-time inventory.
     *
     * @param numberOfParticles Number of birth times to generate.
     */
    void buildBirthTimeInventory(size_t numberOfParticles);

    /**
     * @brief Samples an absolute value from a Gaussian truncated at a limit.
     *
     * @param sigma Standard deviation of the Gaussian.
     * @param limit Maximum accepted absolute value.
     * @return Sampled non-negative time offset.
     */
    double sampleTruncatedHalfGaussian(double sigma, double limit);

    /**
     * @brief Converts old-OPAL pulse time to tracker birth time.
     *
     * @param opalPulseTime Time coordinate in the old-OPAL pulse convention.
     * @return Birth time in the OPALX tracker convention.
     */
    double toBirthTime(double opalPulseTime) const;

    /**
     * @brief Computes the local contiguous subrange of a global emission batch.
     *
     * @param totalToEmit Number of particles emitted globally in the current step.
     * @return Pair of local offset into the global batch and local particle count.
     */
    std::pair<size_t, size_t> computeLocalEmitRange(size_t totalToEmit) const;

public:
    /**
     * @brief Generates the local particles assigned to this rank for one emission step.
     *
     * Particles are sampled on the transverse disk, assigned a per-particle
     * fractional step based on their birth time, and drifted for half of that
     * effective step.
     *
     * @param nlocalBefore Local particle count before allocation.
     * @param globalBegin First global birth-time index assigned to this rank.
     * @param nNew Number of local particles to generate.
     * @param tStart Start time of the current tracker step.
     * @param dt Time step size.
     */
    void generateLocalParticles(
            size_type nlocalBefore, size_t globalBegin, size_t nNew, double tStart, double dt);
};

#endif  // IPPL_OPAL_FLAT_TOP_H
