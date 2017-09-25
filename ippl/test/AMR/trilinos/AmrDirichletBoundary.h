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
    void apply(const AmrIntVect_t& iv,
               typename AmrMultiGridLevel::indices_t& indices,
               typename AmrMultiGridLevel::coefficients_t& values,
               const double& value,
               AmrMultiGridLevel* mglevel,
               const int* nr);
};


template <class AmrMultiGridLevel>
void AmrDirichletBoundary<AmrMultiGridLevel>::apply(const AmrIntVect_t& iv,
                                                    typename AmrMultiGridLevel::indices_t& indices,
                                                    typename AmrMultiGridLevel::coefficients_t& values,
                                                    const double& value,
                                                    AmrMultiGridLevel* mglevel,
                                                    const int* nr)
{
    // find interior neighbour cell
    AmrIntVect_t niv;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        if ( iv[i] > -1 && iv[i] < nr[i] )
            niv[i] = iv[i];
        else
            niv[i] = (iv[i] == -1) ? iv[i] + 1 : iv[i] - 1;
    }
    
    indices.push_back( mglevel->serialize(niv) );
    values.push_back( -value );
}

#endif
