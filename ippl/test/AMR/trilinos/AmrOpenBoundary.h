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
    
//     void gradient_m(const AmrIntVect_t& iv,
//                     const lo_t& dir,
//                     umap_t& map,
//                     const scalar_t& value,
//                     Level* mglevel,
//                     const go_t* nr);
//     
//     void laplace_m(const AmrIntVect_t& iv,
//                    const lo_t& dir,
//                    umap_t& map,
//                    const scalar_t& value,
//                    Level* mglevel,
//                    const go_t* nr);
//     void mixed_m(const AmrIntVect_t& iv,
//                  const lo_t& dir,
//                  umap_t& map,
//                  const scalar_t& value,
//                  Level* mglevel,
//                  const go_t* nr);
    
    
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
        case ABC::Zeroth:
        {
            this->abc0_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        case ABC::First:
        {
            this->abc1_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        case ABC::Second:
        {
            this->abc2_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        case ABC::Third:
        {
            this->abc3_m(iv, dir, map, value, mglevel, nr);
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


/*
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
        niv[dir]  = 0;
        n2iv[dir] = 1;
        
    } else {
        niv[dir]  = nr[dir] - 1;
        n2iv[dir] = nr[dir] - 2;
    }
    
    scalar_t x = this->coordinate_m(niv, 0, mglevel, nr);
    scalar_t y = this->coordinate_m(niv, 1, mglevel, nr);
    scalar_t z = this->coordinate_m(niv, 2, mglevel, nr);
    
    scalar_t x2 = x * x;
    scalar_t y2 = y * y;
    scalar_t z2 = z * z;
    
    scalar_t rho2    = x2 + y2;
    scalar_t invrho2 = 1.0 / rho2;
    scalar_t r2      = rho2 + z2;
    scalar_t invr2   = 1.0 / r2;
    
    scalar_t h = mglevel->cellSize(dir);
    
    scalar_t prefac = 2.0 * h * value;
    
//     std::cout << x << " "
//               << y << " "
//               << z << " "
//               << r2 << " "
//               << rho2 << " "
//               << h << std::endl;
    
    switch ( dir ) {
        
        case 0:
        {
            // x-direction
            scalar_t alpha = - x * invr2;
            scalar_t beta  = - 0.5 * z * ( 3.0 * x2 * x + 4.0 * x * y2 ) * invr2 * invrho2;
            scalar_t gamma = 0.5 * (x2 * x + x * y2) * invr2;
            scalar_t zeta  = 0.5 * invr2 * invrho2 * invrho2 * (4.0 * x2 * x * y * z2 +
                                                                3.0 * x * y2 * y * z2 -
                                                                x * y * rho2 * (rho2 + 2.0 * r2)
                                                               );
            scalar_t eta   = 0.5 * invr2 * invrho2 * ( z2 * y2 * x + x2 * x * r2 );
            scalar_t xi    = - x * y * z * invr2;
            
            // phi
            map[mglevel->serialize(niv)] += alpha * prefac;
            
            map[mglevel->serialize(n2iv)] += value;
            
            // phi_z
            this->gradient_m(niv, 2, map, beta * prefac, mglevel, nr);
            
            // phi_zz
            this->laplace_m(niv, 2, map, gamma * prefac, mglevel, nr);
            
            // phi_y
            this->gradient_m(niv, 1, map, zeta * prefac, mglevel, nr);
            
            // phi_yy
            this->laplace_m(niv, 1, map, eta * prefac, mglevel, nr);
            
            // phi_yz
            this->mixed_m(niv, 0, map, xi * prefac, mglevel, nr);
            
            break;
        }
        
        case 1:
        {
            // y-direction
            scalar_t alpha = -y * invr2;
            scalar_t beta  = -0.5 * z * (3.0 * y * y2 + 4 * y * x2) * invr2 * invrho2;
            scalar_t gamma = 0.5 * invr2 * (y * y2 + y * x2);
            scalar_t zeta  = 0.5 * invr2 * invrho2 * invrho2 * (4.0 * y * y2 * x * z2 +
                                                                3.0 * y * x * x2 * z2 -
                                                                x * y * rho2 * (rho2 + 2.0 * r2)
                                                                );
            scalar_t eta   = 0.5 * invr2 * invrho2 * ( y * x2 * z2 + y * y2 * r2 );
            scalar_t xi    = - x * y * z * invr2;
            
            // phi
            map[mglevel->serialize(niv)] += alpha * prefac;
            
            map[mglevel->serialize(n2iv)] += value;
            
            // phi_z
            this->gradient_m(niv, 2, map, beta * prefac, mglevel, nr);
            
            // phi_zz
            this->laplace_m(niv, 2, map, gamma * prefac, mglevel, nr);
            
            // phi_x
            this->gradient_m(niv, 0, map, zeta * prefac, mglevel, nr);
            
            // phi_xx
            this->laplace_m(niv, 0, map, eta * prefac, mglevel, nr);
            
            // phi_xz
            this->mixed_m(niv, 1, map, xi * prefac, mglevel, nr);
            
            break;
        }
        
        case 2:
        {
            // z-direction
            scalar_t alpha = - z * invr2;
            scalar_t beta  = - 1.5 * invr2 * x * z;
            scalar_t gamma = 0.5 * invr2 * invrho2 * ( x2 * z * z2 + z * y2 * r2 );
            scalar_t zeta  = - 1.5 * invr2 * y * z;
            scalar_t eta   = 0.5 * invr2 * invrho2 * ( y2 * z * z2 + z * x2 * r2 );
            scalar_t xi    = - invr2 * x * y * z;
            
            // phi
            map[mglevel->serialize(niv)] += alpha * prefac;
            
            map[mglevel->serialize(n2iv)] += value;
            
            // phi_x
            this->gradient_m(niv, 0, map, beta * prefac, mglevel, nr);
            
            // phi_xx
            this->laplace_m(niv, 0, map, gamma * prefac, mglevel, nr);
            
            // phi_y
            this->gradient_m(niv, 1, map, zeta * prefac, mglevel, nr);
            
            // phi_yy
            this->laplace_m(niv, 1, map, eta * prefac, mglevel, nr);
            
            // phi_xy
            this->mixed_m(niv, 2, map, xi * prefac, mglevel, nr);
            
            break;
        }
        default:
            throw OpalException("AmrOpenBoundary::abc2_m()",
                                "Infeasible direction");
            break;
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
    bool isLeftValid  = ( iv[dir] != 0 );
    bool isRightValid = ( iv[dir] != nr[dir] - 1 );
    
    scalar_t invh = mglevel->invCellSize(dir);
    scalar_t invh2 = invh * invh;
    
    if ( isLeftValid && isRightValid ) {
        // central stencil
        AmrIntVect_t left = iv;
        left[dir] -= 1;
        
        AmrIntVect_t right = iv;
        right[dir] += 1;
        
        map[mglevel->serialize(left)] += value * invh2;
        map[mglevel->serialize(right)] += value * invh2;
        map[mglevel->serialize(iv)]    += - 2.0 * value * invh2;
    }
    else if ( isLeftValid ) {
        // backward stencil
        
        map[mglevel->serialize(iv)] += 2.0 * invh2 * value;
        
        AmrIntVect_t left = iv;
        left[dir] -= 1;
        map[mglevel->serialize(left)] += -5.0 * invh2 * value;
        
        left[dir] -= 1;
        map[mglevel->serialize(left)] += 4.0 * invh2 * value;
        
        left[dir] -= 1;
        map[mglevel->serialize(left)] += - invh2 * value;
    }
    else if ( isRightValid ) {
        // forward stencil
        
        map[mglevel->serialize(iv)] += 2.0 * invh2 * value;
        
        AmrIntVect_t right = iv;
        right[dir] += 1;
        map[mglevel->serialize(right)] += -5.0 * invh2 * value;
        
        right[dir] += 1;
        map[mglevel->serialize(right)] += 4.0 * invh2 * value;
        
        right[dir] += 1;
        map[mglevel->serialize(right)] += - invh2 * value;
        
    } else {
        throw OpalException("AmrOpenBoundary::laplace_m()",
                            "Laplace discretization not possible.");
    }
}


template <class Level>
void AmrOpenBoundary<Level>::mixed_m(const AmrIntVect_t& iv,
                                     const lo_t& dir,
                                     umap_t& map,
                                     const scalar_t& value,
                                     Level* mglevel,
                                     const go_t* nr)
{
    lo_t d1 = (dir + 1) % 3;
    lo_t d2 = (dir + 2) % 3;
    
    
    scalar_t invh1 = mglevel->invCellSize(d1);
    scalar_t invh2 = mglevel->invCellSize(d2);
    
    
    bool isLeft1Valid  = ( iv[d1] != 0 );
    bool isRight1Valid = ( iv[d1] != nr[d1] - 1 );
    
    bool isLeft2Valid  = ( iv[d2] != 0 );
    bool isRight2Valid = ( iv[d2] != nr[d2] - 1 );
    
    if ( isLeft1Valid && isRight1Valid &&
         isLeft2Valid && isRight2Valid )
    {
        // ll = lower-left
        AmrIntVect_t ll = iv;
        ll[d1] -= 1;
        ll[d2] -= 1;
        map[mglevel->serialize(ll)] += 0.25 * invh1 * invh2 * value;
        
        // ur = upper-right
        AmrIntVect_t ur = iv;
        ur[d1] += 1;
        ur[d2] += 1;
        map[mglevel->serialize(ur)] += 0.25 * invh1 * invh2 * value;
        
        // ul = upper-left
        AmrIntVect_t ul = iv;
        ul[d1] -= 1;
        ul[d2] += 1;
        map[mglevel->serialize(ul)] += - 0.25 * invh1 * invh2 * value;
        
        // lr = lower-right
        AmrIntVect_t lr = iv;
        lr[d1] += 1;
        lr[d2] -= 1;
        map[mglevel->serialize(lr)] += - 0.25 * invh1 * invh2 * value;
    }
    else if ( isLeft1Valid && isLeft2Valid ) {
        // corner upper right
        
        map[mglevel->serialize(iv)] += invh1 * invh2 * value;
        
        // ll = lower-left
        AmrIntVect_t ll = iv;
        ll[d1] -= 1;
        ll[d2] -= 1;
        map[mglevel->serialize(ll)] += invh1 * invh2 * value;
        
        // ul = upper-left
        AmrIntVect_t ul = iv;
        ul[d1] -= 1;
        map[mglevel->serialize(ul)] += - invh1 * invh2 * value;
        
        // lr = lower-right
        AmrIntVect_t lr = iv;
        lr[d2] -= 1;
        map[mglevel->serialize(lr)] += - invh1 * invh2 * value;
    }
    else if ( isLeft1Valid && isRight2Valid ) {
        // corner lower right
        
        map[mglevel->serialize(iv)] += invh1 * invh2 * value;
        
        // ll = lower-left
        AmrIntVect_t ll = iv;
        ll[d1] -= 1;
        map[mglevel->serialize(ll)] -= invh1 * invh2 * value;
        
        // ul = upper-left
        AmrIntVect_t ul = iv;
        ul[d1] -= 1;
        ul[d2] += 1;
        map[mglevel->serialize(ul)] += invh1 * invh2 * value;
        
        // lr = upper-right
        AmrIntVect_t lr = iv;
        lr[d2] += 1;
        map[mglevel->serialize(lr)] += - invh1 * invh2 * value;
    }
    else if ( isRight1Valid && isRight2Valid ) {
        // corner lower left
        
        map[mglevel->serialize(iv)] += invh1 * invh2 * value;
        
        // ul = upper-left
        AmrIntVect_t ul = iv;
        ul[d2] += 1;
        map[mglevel->serialize(ul)] -= invh1 * invh2 * value;
        
        // ur = upper-right
        AmrIntVect_t ur = iv;
        ur[d1] += 1;
        ur[d2] += 1;
        map[mglevel->serialize(ur)] += invh1 * invh2 * value;
        
        // lr = lower-right
        AmrIntVect_t lr = iv;
        lr[d1] += 1;
        map[mglevel->serialize(lr)] += - invh1 * invh2 * value;
    }
    else if ( isRight1Valid && isLeft2Valid ) {
        // corner upper left
        
        map[mglevel->serialize(iv)] += invh1 * invh2 * value;
        
        // ur = upper-right
        AmrIntVect_t ur = iv;
        ur[d1] += 1;
        map[mglevel->serialize(ur)] -= invh1 * invh2 * value;
        
        // lr = lower-right
        AmrIntVect_t lr = iv;
        lr[d1] += 1;
        lr[d2] -= 1;
        map[mglevel->serialize(lr)] += invh1 * invh2 * value;
        
        // ll = lower-left
        AmrIntVect_t ll = iv;
        ll[d2] -= 1;
        map[mglevel->serialize(lr)] += - invh1 * invh2 * value;
    }
    else {
        throw OpalException("AmrOpenBoundary::mixed_m()",
                            "Error in mixed 2nd derivative");
    }
}



template <class Level>
void AmrOpenBoundary<Level>::gradient_m(const AmrIntVect_t& iv,
                                        const lo_t& dir,
                                        umap_t& map,
                                        const scalar_t& value,
                                        Level* mglevel,
                                        const go_t* nr)
{
    
    scalar_t denom = mglevel->invCellSize(dir);
    
    if ( iv[dir] == 0 ) {
        // i + 1
        AmrIntVect_t l1iv = iv;
        l1iv[dir] += 1;
        
        map[mglevel->serialize(iv)]   -= denom * value;
        map[mglevel->serialize(l1iv)] += denom * value;
        
        
    } else if ( iv[dir] == nr[dir] - 1 ) {
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
*/

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
