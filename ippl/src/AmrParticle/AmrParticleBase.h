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
    AmrParticleBase();
    
    ~AmrParticleBase() {}
    
    // Update the particle object after a timestep.  This routine will change
    // our local, total, create particle counts properly.
    void update();
    
    
    // Update the particle object after a timestep.  This routine will change
    // our local, total, create particle counts properly.
    void update(const ParticleAttrib<char>& canSwap);
    
    // sort particles based on the grid and level that they belong to
    void sort();
    
    // sort the particles given a sortlist
    void sort(SortList_t &sortlist);
    
private:
    IpplTimings::TimerRef AssignDensityTimer_m;
    IpplTimings::TimerRef SortParticlesTimer_m;
    IpplTimings::TimerRef UpdateParticlesTimer_m;
};

#include "AmrParticleBase.hpp"

#endif
