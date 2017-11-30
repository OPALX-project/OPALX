#ifndef AMR_BOUNDARY_H
#define AMR_BOUNDARY_H

#include <vector>

#include "AmrDefs.h"

template <class AmrMultiGridLevel>
class AmrBoundary {
    
public:
    
    AmrBoundary(int nPoints) : nPoints_m(nPoints) { };
    
    
    /*!
     * Check if we are on the physical boundary
     * @param iv cell to check
     * @param nr is the number of grid points
     */
    bool isBoundary(const AmrIntVect_t& iv, const int* nr) const {
#if AMREX_SPACEDIM == 3
    return ( iv[0] < 0 || iv[0] >= nr[0] ||
             iv[1] < 0 || iv[1] >= nr[1] ||
             iv[2] < 0 || iv[2] >= nr[2] );
#else
    return ( iv[0] < 0 || iv[0] >= nr[0] ||
             iv[1] < 0 || iv[1] >= nr[1] );
#endif
    }
    
    /*!
     * @param iv is the cell where we want to have the boundary value
     * @param map with indices global matrix indices and matrix values
     * @param value matrix entry (coefficients)
     * @param value of matrix entry that is supposed for index
     * @param nr is the number of grid points
     */
    virtual void apply(const AmrIntVect_t& iv,
                       typename AmrMultiGridLevel::umap_t& map,
                       const typename AmrMultiGridLevel::scalar_t& value,
                       AmrMultiGridLevel* mglevel,
                       const int* nr) = 0;
    
    
    const int& getNumberOfPoints() const {
        return nPoints_m;
    }
    
private:
    const int nPoints_m;    ///< Number of points used for boundary
};

#endif
