#ifndef AMR_INTERPOLATER_H
#define AMR_INTERPOLATER_H

#include <AMReX_BoxArray.H>

///< Abstract base class for all coarse to fine cell interpolaters
template <class AmrMultiGridLevel>
class AmrInterpolater {
    
public:
    
    /*!
     * @param nPoints is the number of interpolation points used
     */
    AmrInterpolater(int nPoints) : nPoints_m(nPoints) { }
    
    /*!
     * Number of cell points used for interpolation.
     */
    const int& getNumberOfPoints() const {
        return nPoints_m;
    }
    
    /*!
     * Get the stencil to interpolate a value from coarse to fine level
     * @param iv is the fine cell where we want to have the interpolated value
     * @param map with global matrix indices of coarse level cells and
     * matrix entries of coarse level cells (coefficients)
     * @param scale to apply to matrix values
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary values at physical domain, e.g. Dirichlet, open BC
     */
    virtual void stencil(const AmrIntVect_t& iv,
                         typename AmrMultiGridLevel::umap_t& map,
			 const typename AmrMultiGridLevel::scalar_t& scale,
                         AmrMultiGridLevel* mglevel) = 0;
    
    /*!
     * Coarse-Fine-Interface
     * Get stencil of coarse side
     * @param iv is the coarse cell at the interface (center cell of Laplacian)
     * @param map with global matrix indices of coarse level cells and
     * matrix entries of coarse level cells (coefficients)
     * @param scale of matrix values
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param ba contains all coarse cells that got refined
     * @param riv is the fine cell at the interface
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary values at physical domain, e.g. Dirichlet, open BC
     */
    virtual void coarse(const AmrIntVect_t& iv,
                        typename AmrMultiGridLevel::umap_t& map,
			const typename AmrMultiGridLevel::scalar_t& scale,
                        int dir, int shift, const amrex::BoxArray& ba,
                        const AmrIntVect_t& riv,
                        AmrMultiGridLevel* mglevel) = 0;
    
    /*!
     * Coarse-Fine-Interface
     * Get stencil of fine side
     * @param iv is the fine ghost cell at the interface (on coarse cell that is not
     * refined)
     * @param map with global matrix indices of fine level cells and
     * matrix entries of fine level cells (coefficients)
     * @param scale of matrix values
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param ba contains all coarse cells that got refined
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary avlues at physical domain, e.g. Dirichlet, open BC
     */
    virtual void fine(const AmrIntVect_t& iv,
                      typename AmrMultiGridLevel::umap_t& map,
                      const typename AmrMultiGridLevel::scalar_t& scale,
                      int dir, int shift, const amrex::BoxArray& ba,
                      AmrMultiGridLevel* mglevel) = 0;
    
protected:
    const int nPoints_m;    ///< Number of points used for interpolation
};

#endif
