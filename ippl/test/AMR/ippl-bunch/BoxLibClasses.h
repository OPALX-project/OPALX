#ifndef BOXLIB_CLASSES_H
#define BOXLIB_CLASSES_H

#include "AmrParticle/AmrParticleBase.h"
#include "AmrParticle/ParticleAmrLayout.h"

template <class T, unsigned Dim >
class BoxLibLayout : public ParticleAmrLayout<T, Dim>
{
public:
    typedef typename ParticleAmrLayout<T, Dim>::pair_t pair_t;
    typedef typename ParticleAmrLayout<T, Dim>::pair_iterator pair_iterator;
    typedef typename ParticleAmrLayout<T, Dim>::SingleParticlePos_t SingleParticlePos_t;
    typedef typename ParticleAmrLayout<T, Dim>::Index_t Index_t;
    
    typedef double AmrField_t;
    typedef double AmrFieldContainer_t;
    typedef typename ParticleAmrLayout<T, Dim>::ParticlePos_t ParticlePos_t;
    typedef ParticleAttrib<Index_t> ParticleIndex_t;
    
    BoxLibLayout() { }
    
    
    void update(AmrParticleBase< BoxLibLayout<T,Dim> >& PData, int lev_min = 0,
                const ParticleAttrib<char>* canSwap = 0)
    {
        std::cout << "BoxLibLayout::update()" << std::endl;
    }
    
    void update(IpplParticleBase< BoxLibLayout<T,Dim> >& PData, 
                const ParticleAttrib<char>* canSwap=0)
    {
        std::cout << "IpplBase update" << std::endl;
        //TODO: exit since we need AmrParticleBase with grids and levels for particles for this layout
        //if IpplParticleBase is used something went wrong
    }
};

template <class PLayout>
class BoxLibParticle : public AmrParticleBase<PLayout>
{
public:
    
    BoxLibParticle(PLayout *layout) : AmrParticleBase<PLayout>(layout) {
        this->initializeAmr();
    }
    
    BoxLibParticle() {
        this->initializeAmr();
    }
    
};

#endif