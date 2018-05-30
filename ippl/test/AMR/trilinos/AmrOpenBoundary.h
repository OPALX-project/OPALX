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
    this->abc1_m(iv, dir, map, value, mglevel, nr);
    
//     this->abc2_m(iv, dir, map, value, mglevel, nr);
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
    
    // cell size in direction
    scalar_t h = mglevel->cellSize(dir);
    
    scalar_t coord1 = this->coordinate_m(niv, dir, mglevel, nr);
    
    scalar_t sign = 1.0;
    
    if ( niv[dir] == -1 ) {
        // lower boundary // --> forward difference
        
        niv[dir] = 1;
        map[mglevel->serialize(niv)] += value;
        
        sign = -1.0;
        
//         niv[dir]  = 0;
//         map[mglevel->serialize(niv)]  += (1.0 + h / coord1) * value;
        
    } else {
        // upper boundary // --> backward difference
        niv[dir] = nr[dir] - 2;
        map[mglevel->serialize(niv)] += value;
        
        
//         niv[dir]  = nr[dir] - 1;
//         map[mglevel->serialize(niv)]  += (1.0 - h / coord1) * value;
    }
    
    //scalar_t r = 1.475625 - 0.5 * h; //0.358;
    
    /*
     * 
     */
//     map[mglevel->serialize(n2iv)] += value;
//     map[mglevel->serialize(niv)] -= 2.0 / n * value;
    
//     std::cout << dir << " " << iv << " " << niv << std::endl;
//     

    
    
//     map[mglevel->serialize(n2iv)] += value;
    
    // gradients other directions
    for (lo_t i = 1; i < 3; ++i) {
        lo_t d = (dir + i) % 3;
        
        
        scalar_t hd = mglevel->cellSize(d);
        
        scalar_t coord2 = this->coordinate_m(niv, d, mglevel, nr);
        
        scalar_t scale = - sign * 2.0 * h * coord2 / coord1 * value;
        
        std::cout << dir << " " << niv << " " << d << " " << coord1 << " " << coord2 << " " << scale << std::endl;
        
        this->gradient_m(niv, d, map, scale, mglevel, nr);
    }
    
    
    // 1st order
// //     map[mglevel->serialize(niv)] += 2.0 * r / (2.0 * r + h) * value;
//     map[mglevel->serialize(niv)] -= 2.0 * h / r * value;
//     map[mglevel->serialize(n2iv)] += value;
// //    map[mglevel->serialize(niv)] += (1.0 - h / r) * value;
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
    
    scalar_t coord1 = this->coordinate_m(niv, dir, mglevel, nr);
    
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
    
    scalar_t denom = 0.5 / mglevel->cellSize(dir);
    
    if ( iv[dir] == 0 ) {
        // forward differencing
        
        // i + 1
        AmrIntVect_t l1iv = iv;
        l1iv[dir] += 1;
        
        
        // i + 2
        AmrIntVect_t l2iv = l1iv;
        l2iv[dir] += 1;
        
        map[mglevel->serialize(iv)] -= 3.0 * denom * value;
        map[mglevel->serialize(l1iv)] += 4.0 * denom * value;
        map[mglevel->serialize(l2iv)] -= denom * value;
        
        
    } else if ( iv[dir] == nr[dir] - 1 ) {
        // backward differencing
        AmrIntVect_t l1iv = iv;
        
        // i - 1
        l1iv[dir] -= 1;
        
        // i - 2
        AmrIntVect_t l2iv = l1iv;
        l2iv[dir] -= 1;
        
        map[mglevel->serialize(iv)] += 3.0 * denom * value;
        map[mglevel->serialize(l1iv)] -= 4.0 * denom * value;
        map[mglevel->serialize(l2iv)] += denom * value;
        
    } else {
        // central difference
        AmrIntVect_t uiv = iv;
        AmrIntVect_t liv = iv;
        
        // i + 1
        uiv[dir] += 1;
        
//         // i + 2
//         AmrIntVect_t u2iv = uiv;
//         u2iv[dir] += 2;
        
        // i - 1
        liv[dir] -= 1;
        
//         AmrIntVect_t l2iv = liv;
//         l2iv[dir] -= 1;
        
//         scalar_t frac = value / (12.0 * mglevel->cellSize(dir) );
        
//         map[mglevel->serialize(u2iv)] -= frac;
//         map[mglevel->serialize(uiv)]  += 8.0 / frac;
//         map[mglevel->serialize(liv)]  -= 8.0 / frac;
//         map[mglevel->serialize(l2iv)] += frac;
        
        map[mglevel->serialize(uiv)] += denom * value;
        map[mglevel->serialize(liv)] -= denom * value;
    }
}


template <class Level>
typename AmrOpenBoundary<Level>::scalar_t
AmrOpenBoundary<Level>::coordinate_m(const AmrIntVect_t& iv,
                                     const lo_t& dir,
                                     Level* mglevel,
                                     const go_t* nr)
{
    scalar_t lower = mglevel->geom.ProbLo(dir);
    scalar_t upper = mglevel->geom.ProbHi(dir);
    
    scalar_t m = (upper - lower) / ( nr[dir] - 1 );
    
    return m * iv[dir] + lower;
}

#endif
