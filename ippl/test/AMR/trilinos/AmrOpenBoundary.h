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
    
    AmrOpenBoundary(int inc) : AmrBoundary<Level>(4), inc_m(inc) { }
    
    void apply(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               Level* mglevel,
               const go_t* nr);
    
    
private:
    
    int inc_m;
    
    void gradient_m(const AmrIntVect_t& iv,
                    const lo_t& dir,
                    umap_t& map,
                    const scalar_t& value,
                    Level* mglevel,
                    const go_t* nr);
    
    void laplace_m(const AmrIntVect_t& iv,
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
//     this->robin_m(iv, dir, map, value, mglevel, nr);
    this->abc2_m(iv, dir, map, value, mglevel, nr);
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
    scalar_t d = inc_m * h + std::abs(c);
    
    map[mglevel->serialize(niv)] += (1.0 - h / d) * value;
    
    std::cout << dir << " "
              << niv << " "
              << h << " "
              << d << " "
              << (1.0 - h / d)
              << std::endl;
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
    
    map[mglevel->serialize(niv)] += ((1.0 - inc_m / 100.0) * value - prefac / c1);
    
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
    AmrIntVect_t n1iv = iv;
    
    if ( n1iv[dir] == -1 ) {
        n1iv[dir]  = 0;
    } else {
        n1iv[dir]  = nr[dir] - 1;
    }
    
    
    /* check if ABC-2 is possible (due to mixed derivative),
     * otherwise do ABC-1
     */
    bool isPossible = true;
    for (lo_t i = 1; i < 3; ++i) {
        lo_t d = (dir + i) % 3;
        
        if ( n1iv[d] == 0 || n1iv[d] == nr[d] - 1 )
            isPossible = false;
    }
    
    if ( !isPossible ) {
        // do ABC-1 instead
        this->abc1_m(iv, dir, map, value, mglevel, nr);
        
        return;
    }
    
    // do ABC-2
    AmrIntVect_t n2iv = n1iv;
    
    if ( iv[dir] == -1 ) {
        n2iv[dir] = 1;
    } else {
        n2iv[dir] = nr[dir] - 2;
    }
    
    scalar_t h1 = mglevel->cellSize(dir);
    
    scalar_t c1 = this->coordinate_m(n1iv, dir, mglevel, nr);
    
    scalar_t prefac = 1.0 / ( 1.0 + 0.5 * c1 / h1 );
    
    map[mglevel->serialize(n2iv)] += prefac * value * (1.0 - 0.5 * c1 / h1);
    map[mglevel->serialize(n1iv)] += prefac * value * c1 / h1;
    
    prefac = 2.0 * h1 * value;
    
    map[mglevel->serialize(n1iv)] -= 0.5 * prefac / c1;
    
    for (lo_t i = 1; i < 3; ++i) {
        lo_t d = (dir + i) % 3;
        
        scalar_t c2 = this->coordinate_m(n1iv, d, mglevel, nr);
        
        scalar_t fac = 0.25 * prefac * c2 * c2;
        
        this->laplace_m(n1iv, d, map, fac, mglevel, nr);
    }
    
    // mixed derivative
    lo_t d2 = (dir + 1) % 3;
    lo_t d3 = (dir + 2) % 3;
    
    scalar_t c2 = this->coordinate_m(n1iv, d2, mglevel, nr);
    scalar_t c3 = this->coordinate_m(n1iv, d3, mglevel, nr);
    
    scalar_t h2 = mglevel->cellSize(d2);
    scalar_t h3 = mglevel->cellSize(d3);
    
    scalar_t constant = 0.25 / ( h2 * h3 ) * 0.5 * c2 * c3 / c1 * prefac;
    
    scalar_t sign = 1.0;
    
    for (int i = -1; i < 2; i += 2) {
        for (int j = -1; j < 2; j += 2) {
            AmrIntVect_t miv = n1iv;
            
            miv[d2] += i;
            miv[d3] += j;
            
            if ( i == j )
                sign = 1.0;
            else
                sign = -1.0;
            
            map[mglevel->serialize(miv)] += sign * constant;
        }
    }
}


template <class Level>
void AmrOpenBoundary<Level>::laplace_m(const AmrIntVect_t& iv,
                                       const lo_t& dir,
                                       umap_t& map,
                                       const scalar_t& value,
                                       Level* mglevel,
                                       const go_t* nr)
{
    AmrIntVect_t lower = iv;
    lower[dir] -= 1;
    
    AmrIntVect_t upper = iv;
    upper[dir] += 1;
    
    scalar_t h = mglevel->cellSize(dir);
    
    scalar_t invh2 = 1.0 / (h * h);
    
    map[mglevel->serialize(lower)] += value * invh2;
    map[mglevel->serialize(upper)] += value * invh2;
    map[mglevel->serialize(iv)]    += - 2.0 * value * invh2;
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
        // i + 1
        AmrIntVect_t l1iv = iv;
        l1iv[dir] += 1;
        
        map[mglevel->serialize(iv)]   += denom * value;
        map[mglevel->serialize(l1iv)] -= denom * value;
        
        
    } else if ( iv[dir] == nr[dir] - 1 ) {
        // backward differencing
        AmrIntVect_t l1iv = iv;
        
        // i - 1
        l1iv[dir] -= 1;
        map[mglevel->serialize(iv)]   += denom * value;
        map[mglevel->serialize(l1iv)] -= denom * value;
        
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
