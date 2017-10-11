#ifndef AMR_LAGRANGE_INTERPOLATER_H
#define AMR_LAGRANGE_INTERPOLATER_H

#include "AmrInterpolater.h"

#if BL_SPACEDIM == 3
    #include <bitset>
    #include <iterator>
    #include <array>
#endif

template <class AmrMultiGridLevel>
class AmrLagrangeInterpolater : public AmrInterpolater<AmrMultiGridLevel>
{
public:
    
    enum Order {
        LINEAR = 1,
        QUADRATIC
    };
    
#if BL_SPACEDIM == 3
    typedef std::bitset<25> bits_t; ///< for checking the neighbour cells
    typedef std::array<unsigned int long, 9> qpattern_t;    ///< quadratic pattern
    typedef std::array<unsigned int long, 4> lpattern_t     ///< linear pattern
#endif
    
public:
    
    AmrLagrangeInterpolater(Order order);
    
    void stencil(const AmrIntVect_t& iv,
                 typename AmrMultiGridLevel::indices_t& indices,
                 typename AmrMultiGridLevel::coefficients_t& values,
                 AmrMultiGridLevel* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                typename AmrMultiGridLevel::indices_t& indices,
                typename AmrMultiGridLevel::coefficients_t& values,
                int dir, int shift, const amrex::BoxArray& ba,
                const AmrIntVect_t& riv,
                AmrMultiGridLevel* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              typename AmrMultiGridLevel::indices_t& indices,
              typename AmrMultiGridLevel::coefficients_t& values,
              int dir, int shift, const amrex::BoxArray& ba,
              AmrMultiGridLevel* mglevel);
    
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
     * @param ba contains all coarse cells that got refined
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary avlues at physical domain, e.g. Dirichlet, open BC
     */
    void fineLinear_m(const AmrIntVect_t& iv,
                  typename AmrMultiGridLevel::indices_t& indices,
                  typename AmrMultiGridLevel::coefficients_t& values,
                  int dir, int shift, const amrex::BoxArray& ba,
                  AmrMultiGridLevel* mglevel);
    
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
     * @param ba contains all coarse cells that got refined
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary avlues at physical domain, e.g. Dirichlet, open BC
     */
    void fineQuadratic_m(const AmrIntVect_t& iv,
                     typename AmrMultiGridLevel::indices_t& indices,
                     typename AmrMultiGridLevel::coefficients_t& values,
                     int dir, int shift, const amrex::BoxArray& ba,
                     AmrMultiGridLevel* mglevel);
    
    /*!
     * First oder interpolation on coarse cell interface side
     * @param iv is the coarse cell at the interface (center cell of Laplacian)
     * @param indices global matrix indices of coarse level cells
     * @param values matrix entries of coarse level cells (coefficients)
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param ba contains all coarse cells that got refined
     * @param riv is the fine cell at the interface
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary values at physical domain, e.g. Dirichlet, open BC
     */
    void crseLinear_m(const AmrIntVect_t& iv,
                      typename AmrMultiGridLevel::indices_t& indices,
                      typename AmrMultiGridLevel::coefficients_t& values,
                      int dir, int shift, const amrex::BoxArray& ba,
                      const AmrIntVect_t& riv,
                      AmrMultiGridLevel* mglevel);
    
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
     * @param ba contains all coarse cells that got refined
     * @param riv is the fine cell at the interface
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary values at physical domain, e.g. Dirichlet, open BC
     */
    void crseQuadratic_m(const AmrIntVect_t& iv,
                         typename AmrMultiGridLevel::indices_t& indices,
                         typename AmrMultiGridLevel::coefficients_t& values,
                         int dir, int shift, const amrex::BoxArray& ba,
                         const AmrIntVect_t& riv,
                         AmrMultiGridLevel* mglevel);
    
#if BL_SPACEDIM == 3
private:
    static constexpr qpattern_t qpattern_ms {
        473536,                             ///< cross pattern
        14798,                              ///< T pattern
        15153152,                           ///< T on head pattern
        918428,                             ///< upper right corner pattern
        7399,                               ///< upper left corner pattern
        7576576,                            ///< mirrored L pattern
        30306304,                           ///< L pattern
        947072,                             ///< left hammer pattern
        236768                              ///< right hammer pattern
    };
    
    static constexpr lpattern_t lpattern_ms {
        54,                                 ///< corner top left pattern
        27,                                 ///< corner top right pattern
        432,                                ///< corner bottom left pattern
        216                                 ///< corner bottom right pattern
    };
#endif
};

#include "AmrLagrangeInterpolater.hpp"

#endif
