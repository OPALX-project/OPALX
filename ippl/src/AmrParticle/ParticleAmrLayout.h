#ifndef PARTICLE_AMR_LAYOUT_H
#define PARTICLE_AMR_LAYOUT_H

#include "Particle/ParticleLayout.h"
// #include "AmrParticle/AmrParticleBase.h"

// #include <memory>

template<class T, unsigned Dim>
class ParticleAmrLayout : public ParticleLayout<T, Dim> {
    
public:
    // pair iterator definition ... this layout does not allow for pairlists
    typedef int pair_t;
    typedef pair_t* pair_iterator;
    typedef typename ParticleLayout<T, Dim>::SingleParticlePos_t
    SingleParticlePos_t;
    typedef typename ParticleLayout<T, Dim>::Index_t Index_t;
  
    // type of attributes this layout should use for position and ID
    typedef ParticleAttrib<SingleParticlePos_t> ParticlePos_t;
    typedef ParticleAttrib<Index_t>             ParticleIndex_t;
    
public:
/*    
    AmrLayout* getAmrLayout() {
        return layout_mp->get();
    }
    
    const AmrLayout* getAmrLayout() const {
        return layout_mp->get();
    }*/
    
//     void update(ParticleBase_t& PData, int minLevel = 0,
//                 const ParticleAttrib<char> canSwap = 0)
//     {
// //         layout_mp->update(PData, minLevel, canSwap);
//     }
    
// private:
//     std::unique_ptr<AmrLayout> layout_mp;
};

#endif
