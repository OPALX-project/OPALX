#ifndef OPALX_SAMPLING_BASE_H
#define OPALX_SAMPLING_BASE_H

#include <memory>
#include "Distribution.h"
#include "Ippl.h"
#include "PartBunch/BunchStateHandler.h"

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t    = FieldContainer<double, 3>;
using Distribution_t      = Distribution;

using view_type = typename ippl::detail::ViewType<Vector_t<double, 3>, 1>::view_type;

class SamplingBase {
protected:
    std::shared_ptr<ParticleContainer_t> pc_m;
    std::shared_ptr<FieldContainer_t> fc_m;
    Distribution_t* opalDist_m;
    std::string samplingMethod_m;
    /// Emission source offset: position R0, momentum P0, start time t0 (applied in sample step).
    Vector_t<double, 3> R0_m    = 0.0;
    Vector_t<double, 3> P0_m    = 0.0;
    double t0_m                 = 0.0;
    std::string emissionModel_m = "NONE";

    /// For one-shot emitters (e.g. Gaussian at delayed t0): guard to avoid double sampling.
    bool hasEmittedOnce_m = false;

public:
    SamplingBase(
            std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
            Distribution_t* dist)
        : pc_m(pc), fc_m(fc), opalDist_m(dist) {}

    SamplingBase(std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc)
        : pc_m(pc), fc_m(fc) {}

    SamplingBase(std::shared_ptr<ParticleContainer_t> pc) : pc_m(pc) {}

    virtual ~SamplingBase() {}

    void setEmissionOffsets(
            ippl::Vector<double, 3> R0, ippl::Vector<double, 3> P0, double t0,
            const std::string& emissionModel = "NONE") {
        R0_m            = R0;
        P0_m            = P0;
        t0_m            = t0;
        emissionModel_m = emissionModel;
    }

    Vector_t<double, 3> getEmissionR0() const { return R0_m; }

    virtual void generateParticles(size_t& /*numberOfParticles*/, Vector_t<double, 3> /*nr*/) {}

    virtual void emitParticles(double /*t*/, double /*dt*/) {}

    /// @brief Whether this sampler has finished all emission (no more particles will be created).
    /// @param t Current simulation time (s).
    virtual bool isEmissionDone(double t) const {
        (void)t;
        return hasEmittedOnce_m;
    }

    // testNumEmitParticles is purely made for testing and should be removed
    virtual void testNumEmitParticles(size_t /*nsteps*/, double /*dt*/) {}

    // testEmitParticles is purely made for testing and should be removed
    virtual void testEmitParticles(size_t /*nsteps*/, double /*dt*/) {}

    virtual void initDomainDecomp(double /*BoxIncr*/) {}

    virtual void setWithDomainDecomp(bool /*withDomainDecomp*/) {}

    /**
     * @brief Computes the number of particles this rank should emit so that the global
     *        total equals totalToSample and no rank exceeds its capacity (space left).
     * @param totalToSample Global number of particles to emit this timestep.
     * @return Local number of particles this rank should create/emit.
     */
    size_t computeLocalEmitCount(size_t totalToSample) const;
};
#endif  // OPALX_SAMPLING_BASE_H
