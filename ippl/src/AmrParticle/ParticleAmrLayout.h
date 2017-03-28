#ifndef PARTICLE_AMR_LAYOUT_H
#define PARTICLE_AMR_LAYOUT_H

#include "Particle/ParticleLayout.h"
#include "AmrParticle/AmrParticleBase.h"

#include <memory>

template<class T, unsigned Dim, class AmrLayout>
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
    
    typedef AmrParticleBase<ParticleAmrLayout<T, Dim, AmrLayout> > ParticleBase_t;
    typedef AmrLayout<T, Dim>::AmrField_t AmrField_t;
    
public:
    
    AmrLayout* getAmrLayout() {
        return layout_mp->get();
    }
    
    const AmrLayout* getAmrLayout() const {
        return layout_mp->get();
    }
    
    void update(ParticleBase_t& PData, int minLevel = 0,
                const ParticleAttrib<char> canSwap = 0)
    {
//         unsigned N = Ippl::getNodes();
//         unsigned myN = Ippl::myNode();
    
//         std::size_t LocalNum = PData.getLocalNum();
//         std::size_t DestroyNum = PData.getDestroyNum();
//         std::size_t TotalNum;
//         int node;
        
        layout_mp->update(PData, minLevel, canSwap);
    }
    
    /*!
     * Put data from particle onto the grid(s)
     * @param f is the field data on the grid
     * @param pp is the particle attribute that is interpolated to the grid
     * @param cache
     * @param lbase is the start level
     * @param lfine is the end level
     */
    template <unsigned Dim, class PT, class CacheData>
    void scatter(AmrField_t& f, const ParticleAttrib<Vektor<PT, Dim>& pp,
                 ParticleAttrib<CacheData>& cache,
                 int lbase = 0, int lfine = -1) const
    {
        layout_mp->scatter(f, pp, cache, lbase, lfine);
    }
    
    /*!
     * Get data from grid(s) to the particles
     * @param f is the field data on the grid
     * @param pp is the particle attribute that is interpolated from the grid(s)
     * @param cache
     * @param lbase is the start level
     * @param lfine is the end level
     */
    template <unsigned Dim, class PT, class CacheData>
    void gather(AmrField_t& f, const ParticleAttrib<Vektor<PT, Dim>& pp,
                ParticleAttrib<CacheData>& cache,
                int lbase = 0, int lfine = -1) const
    {
        layout_mp->gather(f, pp, cache, lbase, lfine);
    }
    
private:
    std::unique_ptr<AmrLayout> layout_mp;
};

#endif