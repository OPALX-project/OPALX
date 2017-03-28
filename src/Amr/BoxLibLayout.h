#ifndef BOXLIB_LAYOUT_H
#define BOXLIB_LAYOUT_H

#include "AmrParticle/ParticleAmrLayout.h"

#include "Amr/AmrBoxLib.h"

template<class T, unsigned Dim>
class BoxLibLayout : public ParGDB {
    
public:
    typedef AmrBoxLib::AmrField_t AmrField_t;
    typedef typename ParticleAmrLayout<T, Dim, BoxLibLayout>::ParticleBase_t
    ParticleBase_t;
    
public:
    
    void update(ParticleBase_t& PData, const ParticleAttrib<char> canSwap = 0);
    
    template <class PT, class CacheData>
    void scatter(AmrField_t& f, const ParticleAttrib<Vektor<PT, Dim>& pp,
                 ParticleAttrib<CacheData>& cache,
                 int lbase, int lfine) const;
    
    template <unsigned Dim, class PT, class CacheData>
    void gather(AmrField_t& f, const ParticleAttrib<Vektor<PT, Dim>& pp,
                ParticleAttrib<CacheData>& cache,
                int lbase, int lfine) const
    
private:
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    //get the cell of the particle
    static IntVect Index (ParticleBase_t& p, 
                          const unsigned int ip,
                          const Geometry& geom);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    //get the cell of the particle
    static IntVect Index (SingleParticlePos_t &R, const Geometry& geom);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Checks/sets a particles location on levels lev_min and higher.
    // Returns false if the particle does not exist on that level.
    static bool Where_m(ParticleBase_t& prt,
                        const unsigned int ip,
                        int lev_min = 0, int finest_level = -1);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Checks/sets whether the particle has crossed a periodic boundary in such a way
    // that it is on levels lev_min and higher.
    static bool PeriodicWhere_m(ParticleBase_& prt,
                                const unsigned int ip,
                                int lev_min = 0, int finest_level = -1);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Checks/sets whether a particle is within its grid (including grow cells).
    static bool RestrictedWhere_m(ParticleBase_t& p,
                                  const unsigned int ip,
                                  int ngrow);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Returns true if the particle was shifted.
    static bool PeriodicShift_m(SingleParticlePos_t R);
};

#include "BoxLibLayout.hpp"

#endif
