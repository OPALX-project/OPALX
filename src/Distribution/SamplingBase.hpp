#ifndef IPPL_SAMPLING_BASE_H
#define IPPL_SAMPLING_BASE_H

#include "Distribution.h"
#include "Ippl.h"
#include <memory>

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t = FieldContainer<double, 3>;
using Distribution_t = Distribution;

class SamplingBase{
protected:
    std::shared_ptr<ParticleContainer_t> pc_m;
    std::shared_ptr<FieldContainer_t> fc_m;
    std::shared_ptr<Distribution_t> opalDist_m;
    std::string samplingMethod_m;
    /// Emission source offset: position R0, momentum P0, start time t0 (applied in sample step).
    ippl::Vector<double, 3> R0_m = {0.0, 0.0, 0.0};
    ippl::Vector<double, 3> P0_m = {0.0, 0.0, 0.0};
    double t0_m = 0.0;

public:
    
    SamplingBase(std::shared_ptr<ParticleContainer_t> pc,
        std::shared_ptr<FieldContainer_t> fc,
        std::shared_ptr<Distribution_t> dist
    )
        : pc_m(pc), fc_m(fc), opalDist_m(dist) {
    }

    SamplingBase(std::shared_ptr<ParticleContainer_t> pc,
        std::shared_ptr<FieldContainer_t> fc
    )
        : pc_m(pc), fc_m(fc) {
    }

    SamplingBase(std::shared_ptr<ParticleContainer_t> pc)
        : pc_m(pc) {
    }
    
    virtual ~SamplingBase() {}

    void setEmissionOffsets(ippl::Vector<double, 3> R0, ippl::Vector<double, 3> P0, double t0) {
        R0_m = R0;
        P0_m = P0;
        t0_m = t0;
    }

    virtual void generateParticles(size_t& /*numberOfParticles*/, Vector_t<double, 3> /*nr*/) {}

    virtual void emitParticles(double /*t*/, double /*dt*/) {}

    // testNumEmitParticles is purely made for testing and should be removed
    virtual void testNumEmitParticles(size_t /*nsteps*/, double /*dt*/) {}

    // testEmitParticles is purely made for testing and should be removed
    virtual void testEmitParticles(size_t /*nsteps*/, double /*dt*/) {}

    virtual void initDomainDecomp(double /*BoxIncr*/) {}

    virtual void setWithDomainDecomp(bool /*withDomainDecomp*/) {}
};
#endif

