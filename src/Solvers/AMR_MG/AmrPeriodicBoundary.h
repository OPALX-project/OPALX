#ifndef AMR_PERIODIC_BOUNDARY_H
#define AMR_PERIODIC_BOUNDARY_H

#include "AmrBoundary.h"

/*!
 * Dirichlet boundary is on faces of physical domain the boundary
 * value would be at different locations depending on the level.
 */
template <class Level>
class AmrPeriodicBoundary : public AmrBoundary<Level> {

public:
    typedef typename Level::umap_t      umap_t;
    typedef typename Level::lo_t        lo_t;
    typedef typename Level::go_t        go_t;
    typedef typename Level::scalar_t    scalar_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;

public:
    
    AmrPeriodicBoundary() : AmrBoundary<Level>(1) { }
    
    void apply(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               Level* mglevel,
               const go_t* nr);
};


template <class Level>
void AmrPeriodicBoundary<Level>::apply(const AmrIntVect_t& iv,
                                       const lo_t& dir,
                                       umap_t& map,
                                       const scalar_t& value,
                                       Level* mglevel,
                                       const go_t* nr)
{
    // find interior neighbour cell on opposite site
    AmrIntVect_t niv = iv;
    niv[dir] = ( iv[dir] == -1 ) ? nr[dir] - 1 : 0;
    
    map[mglevel->serialize(niv)] += value;
}


#endif
