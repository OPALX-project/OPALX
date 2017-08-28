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
     * @param numEntries in matrix row (increase when adding index)
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary values at physical domain, e.g. Dirichlet, open BC
     */
    virtual void stencil(const AmrIntVect_t& iv,
                         typename AmrMultiGridLevel::indices_t& indices,
                         typename AmrMultiGridLevel::coefficients_t& values,
                         int& numEntries,
                         AmrMultiGridLevel* mglevel) = 0;
};

#endif
