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
    
    AmrPeriodicBoundary() : AmrBoundary<AmrMultiGridLevel>(1) { }
    
    void apply(const AmrIntVect_t& iv,
               typename AmrMultiGridLevel::umap_t& map,
               const typename AmrMultiGridLevel::scalar_t& value,
               AmrMultiGridLevel* mglevel,
               const int* nr);
};


template <class AmrMultiGridLevel>
void AmrPeriodicBoundary<AmrMultiGridLevel>::apply(const AmrIntVect_t& iv,
                                                   typename AmrMultiGridLevel::umap_t& map,
                                                   const typename AmrMultiGridLevel::scalar_t& value,
                                                   AmrMultiGridLevel* mglevel,
                                                   const int* nr)
{
    // find interior neighbour cell on opposite site
    AmrIntVect_t niv;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        niv[d] = ( iv[d] == -1 ) ? nr[d] - 1 : 0;
    }
    
    map[mglevel->serialize(niv)] += value;
}

#endif
