#ifndef IPPL_SAMPLING_BASE_H
#define IPPL_SAMPLING_BASE_H

#include "Distribution.h"
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
    public:
        SamplingBase(std::shared_ptr<ParticleContainer_t> &pc, std::shared_ptr<FieldContainer_t> &fc, std::shared_ptr<Distribution_t> &dist)
           : pc_m(pc), fc_m(fc), opalDist_m(dist) {}

       virtual ~SamplingBase() {}

       virtual void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) {}

       virtual void emitParticles(double t, double dt) {}

       // testNumEmitParticles is purely made for testing and should be removed
       virtual void testNumEmitParticles(size_t nsteps, double dt) {}

       // testEmitParticles is purely made for testing and should be removed
       virtual void testEmitParticles(size_t nsteps, double dt) {}

       virtual void initDomainDecomp(Mesh_t<3> *mesh, FieldLayout_t<3> *FL, double BoxIncr) {}

       virtual void setWithDomainDecomp(bool withDomainDecomp) {}
};
#endif

