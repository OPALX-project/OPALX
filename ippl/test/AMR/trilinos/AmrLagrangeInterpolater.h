#ifndef AMR_LAGRANGE_INTERPOLATER_H
#define AMR_LAGRANGE_INTERPOLATER_H

#include "AmrInterpolater.h"

#if AMREX_SPACEDIM == 3
    #include <bitset>
    #include <iterator>
    #include <array>
#endif

template <class AmrMultiGridLevel>
class AmrLagrangeInterpolater : public AmrInterpolater<AmrMultiGridLevel>
{
public:
    
    typedef typename AmrMultiGridLevel::global_ordinal_t go_t;
    typedef typename AmrMultiGridLevel::lo_t lo_t;
    typedef typename AmrMultiGridLevel::scalar_t scalar_t;
    typedef typename AmrMultiGridLevel::umap_t umap_t;

    enum Order {
        LINEAR = 1,
        QUADRATIC
    };
    
#if AMREX_SPACEDIM == 3
    typedef std::bitset<25> qbits_t; ///< for checking the neighbour cells (quadratic)
    typedef std::bitset<9> lbits_t; ///< for checking the neighbour cells (linear)
    typedef std::array<unsigned int long, 9> qpattern_t;    ///< quadratic pattern
    typedef std::array<unsigned int long, 4> lpattern_t;    ///< linear pattern
#endif
    
public:
    
    AmrLagrangeInterpolater(Order order);
    
    void stencil(const AmrIntVect_t& iv,
                 umap_t& map,
                 const scalar_t& scale,
                 AmrMultiGridLevel* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                umap_t& map,
                const scalar_t& scale,
                lo_t dir, lo_t shift, const amrex::BoxArray& ba,
                const AmrIntVect_t& riv,
                AmrMultiGridLevel* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              umap_t& map,
	      const scalar_t& scale,
              lo_t dir, lo_t shift, const amrex::BoxArray& ba,
              AmrMultiGridLevel* mglevel);
    
private:
    
    /*!
     * First order interpolation on fine cell interface side
     * 
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
    void fineLinear_m(const AmrIntVect_t& iv,
		      umap_t& map,
		      const scalar_t& scale,
		      lo_t dir, lo_t shift, const amrex::BoxArray& ba,
		      AmrMultiGridLevel* mglevel);
    
    /*!
     * Second order interpolation on fine cell interface side
     * 
     * @param iv is the fine ghost cell at the interface (on coarse cell that is not
     * refined)
     * @param map with global matrix indices of fine level cells and
     * values matrix entries of fine level cells (coefficients)
     * @param scale of matrix values
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param ba contains all coarse cells that got refined
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary avlues at physical domain, e.g. Dirichlet, open BC
     */
    void fineQuadratic_m(const AmrIntVect_t& iv,
			 umap_t& map,
			 const scalar_t& scale,
			 lo_t dir, lo_t shift, const amrex::BoxArray& ba,
			 AmrMultiGridLevel* mglevel);
    
    /*!
     * First oder interpolation on coarse cell interface side
     * @param iv is the coarse cell at the interface (center cell of Laplacian)
     * @param map with global matrix indices of coarse level cells and
     * values matrix entries of coarse level cells (coefficients)
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
    void crseLinear_m(const AmrIntVect_t& iv,
                      umap_t& map,
		      const scalar_t& scale,
                      lo_t dir, lo_t shift, const amrex::BoxArray& ba,
                      const AmrIntVect_t& riv,
                      AmrMultiGridLevel* mglevel);
    
    /*!
     * Second order interpolation on coarse cell interface side
     * 
     * @param iv is the coarse cell at the interface (center cell of Laplacian)
     * @param map with global matrix indices of coarse level cells and
     * values matrix entries of coarse level cells (coefficients)
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
    void crseQuadratic_m(const AmrIntVect_t& iv,
                         umap_t& map,
                         const scalar_t& scale,
                         lo_t dir, lo_t shift, const amrex::BoxArray& ba,
                         const AmrIntVect_t& riv,
                         AmrMultiGridLevel* mglevel);
    
#if AMREX_SPACEDIM == 3
private:
    static constexpr qpattern_t qpattern_ms {
        473536,                             ///< cross pattern
        14798,                              ///< T pattern
        236768,                             ///< right hammer pattern
        15153152,                           ///< T on head pattern
        947072,                             ///< left hammer pattern
        7399,                               ///< upper left corner pattern
        29596,                              ///< upper right corner pattern
        7576576,                            ///< mirrored L pattern
        30306304                            ///< L pattern
    };
    
    static constexpr lpattern_t lpattern_ms {
        27,                                 ///< corner top right pattern
        216,                                ///< corner bottom right pattern
        432,                                ///< corner bottom left pattern
        54                                  ///< corner top left pattern
    };
#endif
};

#include "AmrLagrangeInterpolater.hpp"

#endif
