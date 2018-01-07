#ifndef AMR_PERIODIC_BOUNDARY_H
#define AMR_PERIODIC_BOUNDARY_H

#include "AmrBoundary.h"

/*!
 * Dirichlet boundary is on faces of physical domain the boundary
 * value would be at different locations depending on the level.
 */
template <class AmrMultiGridLevel>
class AmrPeriodicBoundary : public AmrBoundary<AmrMultiGridLevel> {

public:
    typedef typename AmrMultiGridLevel::umap_t      umap_t;
    typedef typename AmrMultiGridLevel::lo_t        lo_t;
    typedef typename AmrMultiGridLevel::scalar_t    scalar_t;
    typedef amr::AmrIntVect_t                       AmrIntVect_t;

public:
    
    AmrPeriodicBoundary() : AmrBoundary<AmrMultiGridLevel>(1) { }
    
    void apply(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               AmrMultiGridLevel* mglevel,
               const lo_t* nr);
};


template <class AmrMultiGridLevel>
void AmrPeriodicBoundary<AmrMultiGridLevel>::apply(const AmrIntVect_t& iv,
						   const lo_t& dir,
                                                   umap_t& map,
                                                   const scalar_t& value,
                                                   AmrMultiGridLevel* mglevel,
                                                   const lo_t* nr)
{
    // find interior neighbour cell on opposite site
    AmrIntVect_t niv = iv;
    niv[dir] = ( iv[dir] == -1 ) ? nr[dir] - 1 : 0;
    
    map[mglevel->serialize(niv)] += value;
}


#endif
