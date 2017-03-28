#ifndef BOXLIB_LAYOUT_H
#define BOXLIB_LAYOUT_H

#include "AmrParticle/ParticleAmrLayout.h"
#include <MultiFab.H>

class BoxLibLayout : public ParGDB {
    
public:
    typedef MultiFab AmrField_t;
    typedef ParticleAmrLayout::ParticleBase_t ParticleBase_t;
    
public:
    
    void update(ParticleBase_t& PData, const ParticleAttrib<char> canSwap = 0);
    
    template <unsigned Dim, class PT>
    void scatter(AmrField_t& f,
                 const ParticleAttrib< Vektor<PT, Dim> >& pp,
                 ParticleAttrib<CacheData>& cache) const;
};

#endif
