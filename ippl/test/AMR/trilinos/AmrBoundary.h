#ifndef AMR_BOUNDARY_H
#define AMR_BOUNDARY_H

#include <vector>

#include "AmrDefs.h"

template <class AmrMultiGridLevel>
class AmrBoundary {

public:
    typedef typename AmrMultiGridLevel::umap_t umap_t;
    typedef typename AmrMultiGridLevel::lo_t lo_t;
    typedef typename AmrMultiGridLevel::scalar_t scalar_t;
    
public:

    AmrBoundary(lo_t nPoints) : nPoints_m(nPoints) { };
    
    
    /*!
     * Check if we are on the physical boundary (all directions)
     * @param iv cell to check
     * @param nr is the number of grid points
     */
    bool isBoundary(const AmrIntVect_t& iv, const lo_t* nr) const {
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
     * Check if we are on the physical boundary (certain direction)
     * @param iv cell to check
     * @param nr is the number of grid points
     */
    bool isBoundary(const AmrIntVect_t& iv,
		    const lo_t& dir,
		    const lo_t* nr) const {
	return ( iv[dir] < 0 || iv[dir] >= nr[0] );
    }
    
    /*!
     * Checks all directions if physical / mesh boundary.
     * @param iv is the cell where we want to have the boundary value
     * @param map with indices global matrix indices and matrix values
     * @param value matrix entry (coefficients)
     * @param value of matrix entry that is supposed for index
     * @param nr is the number of grid points
     */
    void apply(const AmrIntVect_t& iv,
	       umap_t& map,
	       const scalar_t& value,
	       AmrMultiGridLevel* mglevel,
	       const lo_t* nr);
    
    /*!
     * Apply boundary in a certain direction.
     * @param iv is the cell where we want to have the boundary value
     * @param dir direction of physical / mesh boundary
     * @param map with indices global matrix indices and matrix values
     * @param value matrix entry (coefficients)
     * @param value of matrix entry that is supposed for index
     * @param nr is the number of grid points
     */
    virtual void apply(const AmrIntVect_t& iv,
		       const lo_t& dir,
                       umap_t& map,
                       const scalar_t& value,
                       AmrMultiGridLevel* mglevel,
                       const lo_t* nr) = 0;
    
    
    const lo_t& getNumberOfPoints() const {
        return nPoints_m;
    }
    
private:
    const lo_t nPoints_m;    ///< Number of points used for boundary
};


template <class AmrMultiGridLevel>
void AmrBoundary<AmrMultiGridLevel>::apply(const AmrIntVect_t& iv,
					   umap_t& map,
					   const scalar_t& value,
					   AmrMultiGridLevel* mglevel,
					   const lo_t* nr)
{
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
	if ( this->isBoundary(iv, d, nr) )
	    this->apply(iv, d, map, value, mglevel, nr);
    }
}

#endif
