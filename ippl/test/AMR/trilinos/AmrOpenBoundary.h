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
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               AmrMultiGridLevel* mglevel,
               const lo_t* nr);
    
};


template <class AmrMultiGridLevel>
void AmrOpenBoundary<AmrMultiGridLevel>::apply(const AmrIntVect_t& iv,
					       const lo_t& dir,
					       umap_t& map,
					       const scalar_t& value,
					       AmrMultiGridLevel* mglevel,
					       const lo_t* nr)
{
    /* depending on boundary we need forward
     * or backward difference for the gradient
     */

    
    // find interior neighbour cells
    AmrIntVect_t niv = iv;
    AmrIntVect_t n2iv = iv;
    
    if ( niv[dir] == -1 ) {
	// lower boundary --> forward difference
	niv[dir]  = 0;
	n2iv[dir] = 1;
    } else {
	// upper boundary --> backward difference
	niv[dir]  = nr[dir] - 1;
	n2iv[dir] = nr[dir] - 2;
    }
    
    // celll size in direction
    scalar_t h = 1.0 / scalar_t(nr[dir]);
    scalar_t r = 0.358;
    
    // 1st order
    map[mglevel->serialize(niv)] += 2.0 * r / (2.0 * r + h) * value;
}

#endif
