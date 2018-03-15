#ifndef AMR_OPEN_BOUNDARY_H
#define AMR_OPEN_BOUNDARY_H

#include "AmrBoundary.h"

template <class Level>
class AmrOpenBoundary : public AmrBoundary<Level> {

public:
    typedef typename Level::umap_t      umap_t;
    typedef typename Level::lo_t        lo_t;
    typedef typename Level::scalar_t    scalar_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    
public:
    
    AmrOpenBoundary() : AmrBoundary<Level>(2) { }
    
    void apply(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               Level* mglevel,
               const lo_t* nr);
    
};


template <class Level>
void AmrOpenBoundary<Level>::apply(const AmrIntVect_t& iv,
                                   const lo_t& dir,
                                   umap_t& map,
                                   const scalar_t& value,
                                   Level* mglevel,
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
    
    // cell size in direction
    scalar_t h = mglevel->cellSize(dir);
    scalar_t r = 1.475625 - 0.5 * h; //0.358;
    
    // 1st order
//     map[mglevel->serialize(niv)] += 2.0 * r / (2.0 * r + h) * value;
    map[mglevel->serialize(niv)] -= 2.0 * h / r * value;
    map[mglevel->serialize(n2iv)] += value;
//    map[mglevel->serialize(niv)] += (1.0 - h / r) * value;
}

#endif
