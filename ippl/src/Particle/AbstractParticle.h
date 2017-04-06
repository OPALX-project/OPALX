#ifndef ABSTRACT_PARTICLE_H
#define ABSTRACT_PARTICLE_H

#include "Particle/ParticleLayout.h"
#include "Particle/ParticleAttrib.h"

template <class T, unsigned Dim>
class AbstractParticle {
    
public:
    typedef typename ParticleLayout<T, Dim>::SingleParticlePos_t
    SingleParticlePos_t;
    typedef typename ParticleLayout<T, Dim>::Index_t Index_t;
    typedef ParticleAttrib<SingleParticlePos_t> ParticlePos_t;
    typedef ParticleAttrib<Index_t>             ParticleIndex_t;
    typedef typename ParticleLayout<T, Dim>::UpdateFlags UpdateFlags;
    typedef typename ParticleLayout<T, Dim>::Position_t Position_t;
    
public:
    
    AbstractParticle() : R_p(0), ID_p(0) {}
    
//     AbstractParticle(ParticlePos_t& R,
//                      ParticleIndex_t& ID) : R_p(&R), ID_p(&ID)
//                      {
//                          std::cout << "AbstractParticle()" << std::endl;
//                     }
    
    virtual void addAttribute(ParticleAttribBase& pa) = 0;
    
    virtual size_t getTotalNum() const = 0;
    virtual size_t getLocalNum() const = 0;
    virtual size_t getDestroyNum() const = 0;
    virtual size_t getGhostNum() const = 0;
    virtual void setTotalNum(size_t n) = 0;
    virtual void setLocalNum(size_t n) = 0;
    
    virtual void update() = 0;
    virtual void update(const ParticleAttrib<char>& canSwap) = 0;
    
    virtual void createWithID(unsigned id) = 0;
    virtual void create(size_t) = 0;
    virtual void globalCreate(size_t np) = 0;
    
    virtual void destroy(size_t, size_t, bool = false) = 0;
    
public:
    ParticlePos_t* R_p;
    ParticleIndex_t* ID_p;
};

#endif
