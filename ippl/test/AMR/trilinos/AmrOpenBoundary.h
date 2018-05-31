#ifndef AMR_OPEN_BOUNDARY_H
#define AMR_OPEN_BOUNDARY_H

#include "AmrBoundary.h"

template <class Level>
class AmrOpenBoundary : public AmrBoundary<Level> {

public:
    typedef typename Level::umap_t      umap_t;
    typedef typename Level::lo_t        lo_t;
    typedef typename Level::go_t        go_t;
    typedef typename Level::scalar_t    scalar_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    
public:
    
    AmrOpenBoundary() : AmrBoundary<Level>(4) { }
    
    void apply(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               Level* mglevel,
               const go_t* nr);
    
    
private:
    
    void gradient_m(const AmrIntVect_t& iv,
                    const lo_t& dir,
                    umap_t& map,
                    const scalar_t& value,
                    Level* mglevel,
                    const go_t* nr);
    
    
    scalar_t coordinate_m(const AmrIntVect_t& iv,
                          const lo_t& dir,
                          Level* mglevel,
                          const go_t* nr);
    
    /*
     * Robin boundary condition (i.e. Dirichlet + Neumann)
     */
    void robin_m(const AmrIntVect_t& iv,
                 const lo_t& dir,
                 umap_t& map,
                 const scalar_t& value,
                 Level* mglevel,
                 const go_t* nr);
    
    /*
     * Asymptotic boundary condition first order (ABC1)
     * BUG Not working in 3 dimensions
     */
    void abc1_m(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               Level* mglevel,
               const go_t* nr);
    
    /*
     * Asymptotic boundary condition second order (ABC2)
     * FIXME Not working
     */
    void abc2_m(const AmrIntVect_t& iv,
                const lo_t& dir,
                umap_t& map,
                const scalar_t& value,
                Level* mglevel,
                const go_t* nr);
};


template <class Level>
void AmrOpenBoundary<Level>::apply(const AmrIntVect_t& iv,
                                   const lo_t& dir,
                                   umap_t& map,
                                   const scalar_t& value,
                                   Level* mglevel,
                                   const go_t* nr)
{
    this->robin_m(iv, dir, map, value, mglevel, nr);
//     this->abc1_m(iv, dir, map, value, mglevel, nr);
    
//     this->abc2_m(iv, dir, map, value, mglevel, nr);
}


template <class Level>
void AmrOpenBoundary<Level>::robin_m(const AmrIntVect_t& iv,
                                     const lo_t& dir,
                                     umap_t& map,
                                     const scalar_t& value,
                                     Level* mglevel,
                                     const go_t* nr)
{
    AmrIntVect_t niv = iv;
    AmrIntVect_t n2iv = iv;
    
    scalar_t sign = 1.0;
    
    if ( iv[dir] == -1 ) {
        // lower boundary
        niv[dir] = 0;
    } else {
        // upper boundary        
        niv[dir] = nr[dir] - 1;
    }
    
    // correct cell size
    scalar_t h = mglevel->cellSize(dir);
    
    // coordinate of inner cell
    scalar_t c = this->coordinate_m(niv, dir, mglevel, nr);
    
    // artificial distance
    scalar_t d = 3.0 * h + std::abs(c);
    
    map[mglevel->serialize(niv)] += (1.0 - h / d) * value;
    
    std::cout << dir << " "
              << niv << " "
              << h << " "
              << d << " "
              << (1.0 - h / d)
              << std::endl;
    
//     map[mglevel->serialize(niv)]  += 2.0 * sign * h / (d) * value;
//     map[mglevel->serialize(n2iv)] += value;
    
    // 1st order
    //scalar_t r = 1.475625 - 0.5 * h; //0.358;
// //     map[mglevel->serialize(niv)] += 2.0 * r / (2.0 * r + h) * value;
//     map[mglevel->serialize(niv)] -= 2.0 * h / r * value;
//     map[mglevel->serialize(n2iv)] += value;
// //    map[mglevel->serialize(niv)] += (1.0 - h / r) * value;
}


template <class Level>
void AmrOpenBoundary<Level>::abc1_m(const AmrIntVect_t& iv,
                                    const lo_t& dir,
                                    umap_t& map,
                                    const scalar_t& value,
                                    Level* mglevel,
                                    const go_t* nr)
{
    /* depending on boundary we need forward
     * or backward difference for the gradient
     */

    
    // find interior neighbour cells
    AmrIntVect_t niv = iv;
    AmrIntVect_t n2iv = iv;
    
    
    if ( niv[dir] == -1 ) {
        niv[dir]  = 0;
    } else {
        niv[dir]  = nr[dir] - 1;
    }
    
    // coordinate of inner cell
    scalar_t c1 = this->coordinate_m(niv, dir, mglevel, nr);
    
    // cell size in direction
    scalar_t h1 = mglevel->cellSize(dir);
    
    scalar_t prefac = h1 * value;
    
    map[mglevel->serialize(niv)] += (1.0 - prefac / c1);
    
    // gradients other directions
    for (lo_t i = 1; i < 3; ++i) {
        lo_t d = (dir + i) % 3;
        
        scalar_t c2 = this->coordinate_m(niv, d, mglevel, nr);
        
        scalar_t fac = - prefac * c2 / c1;
        
        this->gradient_m(niv, d, map, fac, mglevel, nr);
        
        
        std::cout << dir << " " << niv << " " << d << " " << c1 << " " << c2 << std::endl;
    }
}


template <class Level>
void AmrOpenBoundary<Level>::abc2_m(const AmrIntVect_t& iv,
                                    const lo_t& dir,
                                    umap_t& map,
                                    const scalar_t& value,
                                    Level* mglevel,
                                    const go_t* nr)
{
    // find interior neighbour cells
    AmrIntVect_t niv = iv;
    AmrIntVect_t n2iv = iv;
    
    if ( niv[dir] == -1 ) {
        // lower boundary // --> forward difference
        niv[dir]  = 0;
        n2iv[dir] = 1;
        
    } else {
        // upper boundary // --> backward difference
        niv[dir]  = nr[dir] - 1;
        n2iv[dir] = nr[dir] - 2;
    }
    
    scalar_t coord1 = 10.0 * this->coordinate_m(niv, dir, mglevel, nr);
    
    // cell size in direction
    scalar_t h = mglevel->cellSize(dir);
    
    // normal
    scalar_t prefac = 1.0 / ( 1.0 + coord1 / (2.0 * h) );
    
    map[mglevel->serialize(n2iv)] += prefac * value;
    map[mglevel->serialize(niv)]  += prefac * coord1 / (2.0 * h) * value;
    map[mglevel->serialize(n2iv)] -= prefac * coord1 / (2.0 * h) * value;
    
    
    map[mglevel->serialize(niv)] -= h / coord1;
    
    // tangential
    
    scalar_t coord3 = 1.0;
    
    for (lo_t i = 1; i < 3; ++i) {
        lo_t d = (dir + i) % 3;
        
        
        scalar_t hd = mglevel->cellSize(d);
        
        scalar_t coord2 = this->coordinate_m(niv, d, mglevel, nr);
        
        coord3 *= coord2;
        
        scalar_t scale = -0.5 * h * coord2 * coord2 / coord1 * value;
        
        AmrIntVect_t liv = niv;
        liv[d] += 1;
        
        AmrIntVect_t riv = niv;
        riv[d] -= 1;
        
        map[mglevel->serialize(liv)] += scale / ( hd * hd );
        map[mglevel->serialize(riv)] += scale / ( hd * hd );
        map[mglevel->serialize(niv)] += 2.0 * scale / ( hd * hd );
        
    }
    
    // mixed derivative
    coord3 = h * coord3 / coord1;
    
    
}


template <class Level>
void AmrOpenBoundary<Level>::gradient_m(const AmrIntVect_t& iv,
                                        const lo_t& dir,
                                        umap_t& map,
                                        const scalar_t& value,
                                        Level* mglevel,
                                        const go_t* nr)
{
    
    scalar_t denom = mglevel->cellSize(dir);
    
    if ( iv[dir] == 0 ) {
        // forward differencing
        
        // i + 1
        AmrIntVect_t l1iv = iv;
        l1iv[dir] += 1;
        
        map[mglevel->serialize(iv)]   -= denom * value;
        map[mglevel->serialize(l1iv)] += denom * value;
        
        
    } else if ( iv[dir] == nr[dir] - 1 ) {
        // backward differencing
        AmrIntVect_t l1iv = iv;
        
        // i - 1
        l1iv[dir] -= 1;
        map[mglevel->serialize(iv)]   -= denom * value;
        map[mglevel->serialize(l1iv)] += denom * value;
        
    } else {
        // central difference
        AmrIntVect_t uiv = iv;
        AmrIntVect_t liv = iv;
        
        // i + 1
        uiv[dir] += 1;
        
        // i - 1
        liv[dir] -= 1;
        
        map[mglevel->serialize(uiv)] += 0.5 * denom * value;
        map[mglevel->serialize(liv)] -= 0.5 * denom * value;
    }
}


template <class Level>
typename AmrOpenBoundary<Level>::scalar_t
AmrOpenBoundary<Level>::coordinate_m(const AmrIntVect_t& iv,
                                     const lo_t& dir,
                                     Level* mglevel,
                                     const go_t* nr)
{
    /*
     * subtract / add half cell length in order to get
     * cell centered position
     */
    scalar_t h = mglevel->cellSize(dir);
    
    scalar_t lower = mglevel->geom.ProbLo(dir) + 0.5 * h;
    scalar_t upper = mglevel->geom.ProbHi(dir) - 0.5 * h;
    
    scalar_t m = (upper - lower) / ( nr[dir] - 1 );
    
    return m * iv[dir] + lower;
}

#endif
