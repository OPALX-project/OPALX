#ifndef AMR_PARTICLE_BASE_H
#define AMR_PARTICLE_BASE_H

#include "Ippl.h"

#include "Particle/IpplParticleBase.h"

#include "AmrParticleLevelCounter.h"

/* The derived classes need to extend the base class by subsequent methods.
 * 
 * scatter the data from the given attribute onto the given Field, using
 * the given Position attribute
 * 
 * template <class FT, unsigned Dim, class PT>
 * void scatter(const ParticleAttrib<FT>& attrib, AmrField_t& f,
 *              const ParticleAttrib<Vektor<PT, Dim> >& pp,
 *              int lbase = 0, int lfine = -1) const;
 * 
 * 
 * gather the data from the given Field into the given attribute, using
 * the given Position attribute
 * 
 * template <class FT, unsigned Dim, class PT>
 * void gather(ParticleAttrib<FT>& attrib, const AmrField_t& f,
 *             const ParticleAttrib<Vektor<PT, Dim> >& pp,
 *             int lbase = 0, int lfine = -1) const;
 */


template<class PLayout>
class AmrParticleBase : public IpplParticleBase<PLayout> {

public:
    typedef typename PLayout::ParticlePos_t         ParticlePos_t;
    typedef typename PLayout::ParticleIndex_t       ParticleIndex_t;
    typedef typename PLayout::SingleParticlePos_t   SingleParticlePos_t;
    typedef typename PLayout::AmrField_t            AmrField_t;
    typedef typename PLayout::AmrFieldContainer_t   AmrFieldContainer_t;
    
    typedef long                                    SortListIndex_t;
    typedef std::vector<SortListIndex_t>            SortList_t;
    typedef std::vector<ParticleAttribBase *>       attrib_container_t;

    ParticleIndex_t Level; // m_lev
    ParticleIndex_t Grid;  // m_grid
    
    typedef AmrParticleLevelCounter<size_t, size_t> ParticleLevelCounter_t;
    
public:
    
    AmrParticleBase();
    
    AmrParticleBase(PLayout* layout) : IpplParticleBase<PLayout>(layout) { }
    
    ~AmrParticleBase() {}
    
    //initialize AmrParticleBase class - add level and grid variables to attribute list
    void initializeAmr() {
        this->addAttribute(Level);
        this->addAttribute(Grid);
    }
    
    const ParticleLevelCounter_t& getLocalNumPerLevel() const;
    
    void setLocalNumPerLevel(const ParticleLevelCounter_t& LocalNumPerLevel);
    
    // Update the particle object after a timestep.  This routine will change
    // our local, total, create particle counts properly.
    void update();
    
    /*!
     * There's is NO check performed if lev_min <= lev_max and
     * lev_min >= 0.
     * @param lev_min is the start level to update
     * @param lev_max is the last level to update
     */
    void update(int lev_min, int lev_max);
    
    // Update the particle object after a timestep.  This routine will change
    // our local, total, create particle counts properly.
    void update(const ParticleAttrib<char>& canSwap);
    
    // sort particles based on the grid and level that they belong to
    void sort();
    
    // sort the particles given a sortlist
    void sort(SortList_t &sortlist);
    
    PLayout& getAmrLayout() { return this->getLayout(); }
    const PLayout& getAmrLayout() const { return this->getLayout(); }
    
protected:
    IpplTimings::TimerRef UpdateParticlesTimer_m;
    IpplTimings::TimerRef SortParticlesTimer_m;
    
private:
    ParticleLevelCounter_t LocalNumPerLevel_m;
};

#include "AmrParticleBase.hpp"

#endif
