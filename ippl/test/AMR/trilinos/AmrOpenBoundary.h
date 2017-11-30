#ifndef AMR_OPEN_BOUNDARY_H
#define AMR_OPEN_BOUNDARY_H

#include "AmrBoundary.h"

template <class AmrMultiGridLevel>
class AmrOpenBoundary : public AmrBoundary<AmrMultiGridLevel> {

public:
    typedef typename AmrMultiGridLevel::umap_t umap_t;
    typedef typename AmrMultiGridLevel::lo_t lo_t;
    typedef typename AmrMultiGridLevel::scalar_t scalar_t;
    
public:
    
    AmrOpenBoundary() : AmrBoundary<AmrMultiGridLevel>(2) { }
    
    void apply(const AmrIntVect_t& iv,
               umap_t& map,
               const scalar_t& value,
               AmrMultiGridLevel* mglevel,
               const lo_t* nr);
};


template <class AmrMultiGridLevel>
void AmrOpenBoundary<AmrMultiGridLevel>::apply(const AmrIntVect_t& iv,
					       umap_t& map,
					       const scalar_t& value,
					       AmrMultiGridLevel* mglevel,
					       const lo_t* nr)
{
    /* there should be only one boundary at a time, i.e.
     * either x-, y- or z-direction
     */
    
    /* depending on boundary we need forward
     * or backward difference for the gradient
     */
    
    // find interior neighbour cells
    AmrIntVect_t niv = iv;
    AmrIntVect_t n2iv = iv; // next interior cell
    
    lo_t d = 0;
    for ( ; d < AMREX_SPACEDIM; ++d) {
        
        if ( niv[d] == -1 ) {
            // lower boundary --> forward difference
            niv[d] = 0;
            n2iv[d] = 1;
            break;
            
        } else if ( niv[d] == nr[d] ) {
            // upper boundary --> backward difference
            niv[d] = nr[d] - 1;
            n2iv[d] = nr[d] - 2;
            break;
        }
    }
    
    // cell size in direction
    scalar_t h = 1.0 / scalar_t(nr[d]);
    scalar_t r = 0.358;
    
    // 1st order
    map[mglevel->serialize(niv)] += 2.0 * r / (2.0 * r + h) * value;
}

#endif
