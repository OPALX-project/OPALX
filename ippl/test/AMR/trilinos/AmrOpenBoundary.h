#ifndef AMR_OPEN_BOUNDARY_H
#define AMR_OPEN_BOUNDARY_H

#include "AmrBoundary.h"

#include "Utilities/OpalException.h"

template <class Level>
class AmrOpenBoundary : public AmrBoundary<Level> {

public:
    typedef typename Level::umap_t      umap_t;
    typedef typename Level::lo_t        lo_t;
    typedef typename Level::go_t        go_t;
    typedef typename Level::scalar_t    scalar_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    
public:
    
//     AmrOpenBoundary(ABC order) : AmrBoundary<Level>(4), order_m(order) { }
    AmrOpenBoundary(int order) : AmrBoundary<Level>(4), order_m(order) { }
    
    void apply(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               Level* mglevel,
               const go_t* nr);
    
    
    enum ABC {
        Zeroth,
        First,
        Second,
        Third
    };
    
    
private:
    
//     ABC order_m;
    int order_m;
    
    scalar_t coordinate_m(const AmrIntVect_t& iv,
                          const lo_t& dir,
                          Level* mglevel,
                          const go_t* nr);
    
    /*!
     * Asymptotic boundary condition 0th order (ABC0)
     */
    void abc0_m(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               Level* mglevel,
               const go_t* nr);
    
    
    /*!
     * Asymptotic boundary condition 1st order (ABC1)
     */
    void abc1_m(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               Level* mglevel,
               const go_t* nr);
    
    /*!
     * Asymptotic boundary condition 2nd order (ABC2)
     */
    void abc2_m(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               Level* mglevel,
               const go_t* nr);
    
    
    /*!
     * Asymptotic boundary condition 3rd order (ABC3)
     */
    void abc3_m(const AmrIntVect_t& iv,
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
    switch ( order_m ) {
        case ABC::Third:
        {
            this->abc3_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        case ABC::Second:
        {
            this->abc2_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        case ABC::First:
        {
            this->abc1_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        case ABC::Zeroth:
        {
            this->abc0_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        default:
        {
            throw OpalException("AmrOpenBoundary<Level>::apply()",
                                "Order not supported.");
        }
    }
}


template <class Level>
void AmrOpenBoundary<Level>::abc0_m(const AmrIntVect_t& iv,
                                    const lo_t& dir,
                                    umap_t& map,
                                    const scalar_t& value,
                                    Level* mglevel,
                                    const go_t* nr)
{
    // ABC0 == Dirichlet BC
    AmrIntVect_t niv = iv;
    
    if ( iv[dir] == -1 ) {
        // lower boundary
        niv[dir] = 0;
        
    } else {
        // upper boundary        
        niv[dir] = nr[dir] - 1;
    }
    
    map[mglevel->serialize(niv)] -= value;
}


template <class Level>
void AmrOpenBoundary<Level>::abc1_m(const AmrIntVect_t& iv,
                                    const lo_t& dir,
                                    umap_t& map,
                                    const scalar_t& value,
                                    Level* mglevel,
                                    const go_t* nr)
{
    // ABC1 == Robin BC (i.e. Dirichlet BC + Neumann BC)
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
    scalar_t x = this->coordinate_m(niv, 0, mglevel, nr);
    scalar_t y = this->coordinate_m(niv, 1, mglevel, nr);
    scalar_t z = this->coordinate_m(niv, 2, mglevel, nr);
    
    scalar_t r = std::sqrt( x * x + y * y + z * z );
    
    // artificial distance
    scalar_t d = nr[dir] * h + r;
    
    map[mglevel->serialize(niv)] += (1.0 - h / d) * value;
}


template <class Level>
void AmrOpenBoundary<Level>::abc2_m(const AmrIntVect_t& iv,
                                    const lo_t& dir,
                                    umap_t& map,
                                    const scalar_t& value,
                                    Level* mglevel,
                                    const go_t* nr)
{
    AmrIntVect_t niv = iv;
    AmrIntVect_t n2iv = iv;
    
    if ( iv[dir] == -1 ) {
        // lower boundary
        niv[dir] = 0;
        n2iv[dir] = 1;
        
    } else {
        // upper boundary        
        niv[dir] = nr[dir] - 1;
        n2iv[dir] = nr[dir] - 2;
    }
    
    // correct cell size
    scalar_t h = mglevel->cellSize(dir);
    
    // coordinate of inner cell
    scalar_t x = this->coordinate_m(niv, 0, mglevel, nr);
    scalar_t y = this->coordinate_m(niv, 1, mglevel, nr);
    scalar_t z = this->coordinate_m(niv, 2, mglevel, nr);
    
    scalar_t r = std::sqrt(x * x + y * y + z * z);
    
    scalar_t d = nr[dir] * h + r;
    
    map[mglevel->serialize(niv)] += (-2.0 * h * h / ( d * d + 4.0 * h * d) +
                                      2.0 * d / ( 4.0 * h + d) + 4.0 * h / (4.0 * h + d)) * value;
    map[mglevel->serialize(n2iv)] += - d / ( 4.0 * h + d ) * value;
}


template <class Level>
void AmrOpenBoundary<Level>::abc3_m(const AmrIntVect_t& iv,
                                    const lo_t& dir,
                                    umap_t& map,
                                    const scalar_t& value,
                                    Level* mglevel,
                                    const go_t* nr)
{
    //TODO
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
