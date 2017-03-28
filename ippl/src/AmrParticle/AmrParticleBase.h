#ifndef AMR_PARTICLE_BASE_H
#define AMR_PARTICLE_BASE_H

#include "Particle/IpplParticleBase.h"

template<class PLayout>
class AmrParticleBase : public IpplParticleBase<PLayout> {

public:
    typedef typename PLayout::ParticlePos_t   ParticlePos_t;
    typedef typename PLayout::ParticleIndex_t ParticleIndex_t;
    typedef typename PLayout::SingleParticlePos_t SingleParticlePos_t;

    ParticleIndex_t level; // m_lev
    ParticleIndex_t grid;  // m_grid
    
public:
    void update();
    
    void update(const ParticleAttrib<char>& canSwap);
    
private:
    IpplTimings::TimerRef AssignDensityTimer_m;
    IpplTimings::TimerRef SortParticlesTimer_m;
    IpplTimings::TimerRef UpdateParticlesTimer_m;
};

#include "AmrParticleBase.hpp"

#endif
