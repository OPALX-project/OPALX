#ifndef AMR_INTERPOLATER_H
#define AMR_INTERPOLATER_H

template <class AmrMultiGridLevel>
class AmrInterpolater {
    
public:
    
    /*!
     * Number of cell points used for interpolation.
     */
    virtual const int& getNumberOfPoints() const = 0;
    
    /*!
     * Get the stencil to interpolate a value from coarse to fine level
     * @param iv is the fine cell where we want to have the interpolated value
     * @param indices global matrix indices of coarse level cells
     * @param values matrix entries of coarse level cells (coefficients)
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary values at physical domain, e.g. Dirichlet, open BC
     */
    virtual void stencil(const AmrIntVect_t& iv,
                         typename AmrMultiGridLevel::indices_t& indices,
                         typename AmrMultiGridLevel::coefficients_t& values,
                         AmrMultiGridLevel* mglevel) = 0;
    
    /*!
     * Coarse-Fine-Interface
     * Get stencil of coarse side
     * @param iv is the coarse cell at the interface (center cell of Laplacian)
     * @param indices global matrix indices of coarse level cells
     * @param values matrix entries of coarse level cells (coefficients)
     * @param direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary values at physical domain, e.g. Dirichlet, open BC
     */
    virtual void coarse(const AmrIntVect_t& iv,
                        typename AmrMultiGridLevel::indices_t& indices,
                        typename AmrMultiGridLevel::coefficients_t& values,
                        int direction, int shift, AmrMultiGridLevel* mglevel) = 0;
    
    /*!
     * Coarse-Fine-Interface
     * Get stencil of fine side
     * @param iv is the fine ghost cell at the interface (on coarse cell that is not
     * refined)
     * @param indices global matrix indices of fine level cells
     * @param values matrix entries of fine level cells (coefficients)
     * @param direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary avlues at physical domain, e.g. Dirichlet, open BC
     */
    virtual void fine(const AmrIntVect_t& iv,
                      typename AmrMultiGridLevel::indices_t& indices,
                      typename AmrMultiGridLevel::coefficients_t& values,
                      int direction, int shift, AmrMultiGridLevel* mglevel) = 0;
};

#endif
