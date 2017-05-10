#ifndef BOXLIB_LAYOUT_H
#define BOXLIB_LAYOUT_H

#include "AmrParticle/ParticleAmrLayout.h"
#include "AmrParticle/AmrParticleBase.h"

#include "Amr/AmrDefs.h"

#include <ParGDB.H>

template<class T, unsigned Dim>
class BoxLibLayout : public ParticleAmrLayout<T, Dim>,
                     public ParGDB
{
    
public:
    typedef typename ParticleAmrLayout<T, Dim>::pair_t pair_t;
    typedef typename ParticleAmrLayout<T, Dim>::pair_iterator pair_iterator;
    typedef typename ParticleAmrLayout<T, Dim>::SingleParticlePos_t SingleParticlePos_t;
    typedef typename ParticleAmrLayout<T, Dim>::Index_t Index_t;
    
    typedef amr::AmrField_t AmrField_t;
    typedef amr::AmrFieldContainer_t AmrFieldContainer_t;
    typedef typename ParticleAmrLayout<T, Dim>::ParticlePos_t ParticlePos_t;
    typedef ParticleAttrib<Index_t> ParticleIndex_t;

public:
    
    /*!
     * Initializes default Geometry, DistributionMapping and BoxArray.
     */
    BoxLibLayout();
    
    BoxLibLayout(int nGridPoints, int maxGridSize, double lower, double upper);
    
    BoxLibLayout(const Geometry &geom,
                 const DistributionMapping &dmap,
                 const BoxArray &ba);
    
    BoxLibLayout(const Array<Geometry> &geom,
                 const Array<DistributionMapping> &dmap,
                 const Array<BoxArray> &ba,
                 const Array<int> &rr);
    
    void update(IpplParticleBase< BoxLibLayout<T,Dim> >& PData, const ParticleAttrib<char>* canSwap = 0);
    
    void update(AmrParticleBase< BoxLibLayout<T,Dim> >& PData, int lev_min = 0,
                const ParticleAttrib<char>* canSwap = 0);
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    //get the cell of the particle
    /*static*/ IntVect Index (AmrParticleBase< BoxLibLayout<T,Dim> >& p, 
                          const unsigned int ip,
                          const Geometry& geom);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    //get the cell of the particle
    /*static*/ IntVect Index (SingleParticlePos_t &R, const Geometry& geom);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Checks/sets a particles location on levels lev_min and higher.
    // Returns false if the particle does not exist on that level.
    /*static*/ bool Where_m(AmrParticleBase< BoxLibLayout<T,Dim> >& prt,
//                         const PLayout* layout_p,
                        const unsigned int ip,
                        int lev_min = 0, int finest_level = -1);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Checks/sets whether the particle has crossed a periodic boundary in such a way
    // that it is on levels lev_min and higher.
    /*static*/ bool PeriodicWhere_m(AmrParticleBase< BoxLibLayout<T,Dim> >& prt,
//                                 const PLayout* layout_p
                                const unsigned int ip,
                                int lev_min = 0, int finest_level = -1);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Checks/sets whether a particle is within its grid (including grow cells).
    /*static*/ bool RestrictedWhere_m(AmrParticleBase< BoxLibLayout<T,Dim> >& p,
//                                   const PLayout* layout_p
                                  const unsigned int ip,
                                  int ngrow);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Returns true if the particle was shifted.
    /*static*/ bool PeriodicShift_m(SingleParticlePos_t R/*, const PLayout* layout_p*/);
    
    
    void initDefaultBox(int nGridPoints, int maxGridSize,
                        double lower, double upper);
    
};

#include "BoxLibLayout.hpp"

#endif
