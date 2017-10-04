#ifndef AMR_LAGRANGE_INTERPOLATER_H
#define AMR_LAGRANGE_INTERPOLATER_H

#include "AmrInterpolater.h"

template <class AmrMultiGridLevel>
class AmrLagrangeInterpolater : public AmrInterpolater<AmrMultiGridLevel>
{
public:
    
    enum Order {
        LINEAR = 1,
        QUADRATIC
    };
    
    AmrLagrangeInterpolater(Order order);
    
    void stencil(const AmrIntVect_t& iv,
                 typename AmrMultiGridLevel::indices_t& indices,
                 typename AmrMultiGridLevel::coefficients_t& values,
                 AmrMultiGridLevel* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                typename AmrMultiGridLevel::indices_t& indices,
                typename AmrMultiGridLevel::coefficients_t& values,
                int dir, int shift, AmrMultiGridLevel* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              typename AmrMultiGridLevel::indices_t& indices,
              typename AmrMultiGridLevel::coefficients_t& values,
              int dir, int shift, AmrMultiGridLevel* mglevel);
    
private:
    
    /*!
     * First order interpolation on fine cell interface side
     * 
     * @param iv is the fine ghost cell at the interface (on coarse cell that is not
     * refined)
     * @param indices global matrix indices of fine level cells
     * @param values matrix entries of fine level cells (coefficients)
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary avlues at physical domain, e.g. Dirichlet, open BC
     */
    void fineLinear_m(const AmrIntVect_t& iv,
                  typename AmrMultiGridLevel::indices_t& indices,
                  typename AmrMultiGridLevel::coefficients_t& values,
                  int dir, int shift, AmrMultiGridLevel* mglevel);
    
    /*!
     * Second order interpolation on fine cell interface side
     * 
     * @param iv is the fine ghost cell at the interface (on coarse cell that is not
     * refined)
     * @param indices global matrix indices of fine level cells
     * @param values matrix entries of fine level cells (coefficients)
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary avlues at physical domain, e.g. Dirichlet, open BC
     */
    void fineQuadratic_m(const AmrIntVect_t& iv,
                     typename AmrMultiGridLevel::indices_t& indices,
                     typename AmrMultiGridLevel::coefficients_t& values,
                     int dir, int shift, AmrMultiGridLevel* mglevel);
    
    /*
     * First oder interpolation on coarse cell interface side
     * @param iv is the coarse cell at the interface (center cell of Laplacian)
     * @param indices global matrix indices of coarse level cells
     * @param values matrix entries of coarse level cells (coefficients)
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary values at physical domain, e.g. Dirichlet, open BC
     */
    void crseLinear_m(const AmrIntVect_t& iv,
                      typename AmrMultiGridLevel::indices_t& indices,
                      typename AmrMultiGridLevel::coefficients_t& values,
                      int dir, int shift, AmrMultiGridLevel* mglevel);
    
    /*!
     * Second order interpolation on coarse cell interface side
     * 
     * @param iv is the coarse cell at the interface (center cell of Laplacian)
     * @param indices global matrix indices of coarse level cells
     * @param values matrix entries of coarse level cells (coefficients)
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary values at physical domain, e.g. Dirichlet, open BC
     */
    void crseQuadratic_m(const AmrIntVect_t& iv,
                         typename AmrMultiGridLevel::indices_t& indices,
                         typename AmrMultiGridLevel::coefficients_t& values,
                         int dir, int shift, AmrMultiGridLevel* mglevel);
};

#include "AmrLagrangeInterpolater.hpp"

#endif
