#ifndef AMR_DIRICHLET_BOUNDARY_H
#define AMR_DIRICHLET_BOUNDARY_H

#include "AmrBoundary.h"

/*!
 * Dirichlet boundary is on faces of physical domain the boundary
 * value would be at different locations depending on the level.
 */
template <class AmrMultiGridLevel>
class AmrDirichletBoundary : public AmrBoundary<AmrMultiGridLevel> {
    
public:
    typedef typename AmrMultiGridLevel::umap_t umap_t;
    typedef typename AmrMultiGridLevel::lo_t lo_t;
    typedef typename AmrMultiGridLevel::scalar_t scalar_t;
    
public:
    
    AmrDirichletBoundary() : AmrBoundary<AmrMultiGridLevel>(1) { }
    
    void apply(const AmrIntVect_t& iv,
	       const lo_t& dir,
	       umap_t& map,
	       const scalar_t& value,
	       AmrMultiGridLevel* mglevel,
	       const lo_t* nr);
};


template <class AmrMultiGridLevel>
void AmrDirichletBoundary<AmrMultiGridLevel>::apply(const AmrIntVect_t& iv,
						    const lo_t& dir,
						    umap_t& map,
						    const scalar_t& value,
						    AmrMultiGridLevel* mglevel,
						    const lo_t* nr)
{
    // find interior neighbour cell
    AmrIntVect_t niv = iv;
    niv[dir] = (iv[dir] == -1) ? iv[dir] + 1 : iv[dir] - 1;
    map[mglevel->serialize(niv)] -= value;
}




#endif
