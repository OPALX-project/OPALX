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

class OpalFlatTop : public SamplingBase {
public:
    OpalFlatTop(
            std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
            Distribution_t* opalDist);

    OpalFlatTop(
            std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
            bool emitting, double sigmaTFall, double sigmaTRise, Vector_t<double, 3> cutoff,
            double tPulseLengthFWHM, Vector_t<double, 3> sigmaR,
            double ftOscAmplitude = 0.0, double ftOscPeriods = 0.0);

    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;
    void emitParticles(double t, double dt) override;
    bool isEmissionDone(double /*t*/) const override {
        return inventoryBuilt_m && nextGlobalIndex_m >= birthTimes_m.size();
    }
    double getGlobalTimeShift() const override {
        return std::max(0.0, 0.5 * emissionTime_m - t0_m);
    }

    void initDomainDecomp(double BoxIncr) override;
    void setWithDomainDecomp(bool withDomainDecomp) override;

    double getEmissionTime() const { return emissionTime_m; }
    double getDistArea() const { return distArea_m; }
    size_t getInventorySize() const { return birthTimes_m.size(); }
    size_t getNextGlobalIndex() const { return nextGlobalIndex_m; }
    const std::vector<double>& getBirthTimes() const { return birthTimes_m; }

    void setBirthTimesForTest(std::vector<double> birthTimes);

private:
    using size_type = ippl::detail::size_type;

    OpalFlatTopGeneratorPool rand_pool_m;
    std::mt19937_64 host_rng_m;

    double flattopTime_m         = 0.0;
    double normalizedFlankArea_m = 0.0;
    double distArea_m            = 0.0;
    double sigmaTFall_m          = 0.0;
    double sigmaTRise_m          = 0.0;
    Vector_t<double, 3> cutoffR_m;
    double fallTime_m       = 0.0;
    double riseTime_m       = 0.0;
    bool emitting_m         = false;
    size_t totalN_m         = 0;
    bool withDomainDecomp_m = false;
    double emissionTime_m   = 0.0;
    Vector_t<double, 3> nr_m;
    Vector_t<double, 3> hr_m;
    Vector_t<double, 3> sigmaR_m;
    double ftOscAmplitude_m = 0.0;
    double ftOscPeriods_m   = 0.0;

    std::vector<double> birthTimes_m;
    size_t nextGlobalIndex_m = 0;
    bool inventoryBuilt_m    = false;

    static size_t determineRandInit();
    static size_t determineHostSeed();

    void setParameters(Distribution_t* opalDist);
    void setInternalVariables(
            bool emitting, double sigmaTFall, double sigmaTRise, Vector_t<double, 3> cutoff,
            double tPulseLengthFWHM, Vector_t<double, 3> sigmaR, double ftOscAmplitude,
            double ftOscPeriods);

    void buildBirthTimeInventory(size_t numberOfParticles);
    double sampleTruncatedHalfGaussian(double sigma, double limit);
    double toBirthTime(double opalPulseTime) const;
    std::pair<size_t, size_t> computeLocalEmitRange(size_t totalToEmit) const;
    void generateLocalParticles(
            size_type nlocalBefore, size_t globalBegin, size_t nNew, double tStart, double dt);
};

#endif  // IPPL_OPAL_FLAT_TOP_H
