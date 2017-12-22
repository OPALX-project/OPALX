#ifndef AMR_MULTI_GRID_LEVEL
#define AMR_MULTI_GRID_LEVEL

#include <vector>

#include <AMReX_IntVect.H>

#include "Ippl.h"

#include "AmrMultiGridCore.h"

#include <unordered_map>

template <class MatrixType, class VectorType>
class AmrMultiGridLevel {
    
public:
    typedef amrex::MultiFab AmrField_t;
    typedef amrex::Geometry AmrGeometry_t;
    typedef std::unique_ptr<AmrField_t> AmrField_u;
    typedef std::shared_ptr<AmrField_t> AmrField_s;
    typedef amrex::IntVect AmrIntVect_t;
    typedef MatrixType matrix_t;
    typedef VectorType vector_t;
    typedef amrex::BaseFab<int> basefab_t;
    typedef amrex::FabArray<basefab_t> mask_t;
    typedef std::shared_ptr<AmrBoundary<AmrMultiGridLevel<MatrixType,
                                                          VectorType
                                                          >
                                        >
                            > boundary_t;
    
    typedef amr::comm_t comm_t;
    typedef amr::dmap_t dmap_t;
    typedef amr::node_t node_t;
    typedef amr::global_ordinal_t global_ordinal_t;
    typedef amr::scalar_t scalar_t;
    typedef amr::local_ordinal_t lo_t;

    typedef std::vector<lo_t>                  indices_t;
    typedef std::vector<scalar_t>              coefficients_t;
    typedef std::unordered_map<lo_t, scalar_t> umap_t;
    
    // covered   : ghost cells covered by valid cells of this FabArray
    //             (including periodically shifted valid cells)
    // notcovered: ghost cells not covered by valid cells
    //             (including ghost cells outside periodic boundaries) (BNDRY)
    // physbnd   : boundary cells outside the domain (excluding periodic boundaries)
    // interior  : interior cells (i.e., valid cells)
    enum Mask {
        COVERED   = -1,
        INTERIOR  =  0,
        BNDRY     =  1,
        PHYSBNDRY =  2
    };

    // NO   : not a refined cell 
    // YES  : cell got refined
    enum Refined {
        YES = 0,
        NO  = 1
    };
    
public:
    AmrMultiGridLevel(const amrex::BoxArray& _grids,
                      const amrex::DistributionMapping& _dmap,
                      const AmrGeometry_t& _geom,
                      const AmrIntVect_t& rr,
                      const boundary_t* bc,
                      const Teuchos::RCP<comm_t>& comm,
                      const Teuchos::RCP<node_t>& node);
    
    ~AmrMultiGridLevel();
    
    int serialize(const AmrIntVect_t& iv) const;
    
    bool isBoundary(const AmrIntVect_t& iv) const;
    
    bool applyBoundary(const AmrIntVect_t& iv,
                       umap_t& map,
                       const scalar_t& value);
    
    bool applyBoundary(const AmrIntVect_t& iv,
                       const basefab_t& fab,
                       umap_t& map,
                       const scalar_t& value);

    void applyBoundary(const AmrIntVect_t& iv,
                       const lo_t& dir,
                       umap_t& map,
                       const scalar_t& value);
    
    const AmrIntVect_t& refinement() const;
    
private:
    void buildLevelMask_m();
    
    void buildMap_m(const Teuchos::RCP<comm_t>& comm,
                    const Teuchos::RCP<node_t>& node);
    
public:
    const amrex::BoxArray& grids;
    const amrex::DistributionMapping& dmap;
    const AmrGeometry_t& geom;
    
    Teuchos::RCP<dmap_t> map_p;         ///< core map
    
    Teuchos::RCP<matrix_t> Anf_p;       ///< no fine Poisson matrix
    Teuchos::RCP<matrix_t> R_p;         ///< restriction matrix
    Teuchos::RCP<matrix_t> I_p;         ///< interpolation matrix
    Teuchos::RCP<matrix_t> Bcrse_p;     ///< boundary from coarse cells
    Teuchos::RCP<matrix_t> Bfine_p;     ///< boundary from fine cells
    Teuchos::RCP<matrix_t> Awf_p;       ///< composite Poisson matrix
    
    /// gradient matrices in x, y, and z to compute electric field
    Teuchos::RCP<matrix_t> G_p[AMREX_SPACEDIM];
    
    Teuchos::RCP<vector_t> rho_p;       ///< charge density
    Teuchos::RCP<vector_t> phi_p;       ///< potential vector
    Teuchos::RCP<vector_t> residual_p;  ///< residual over all cells
    Teuchos::RCP<vector_t> error_p;     ///< error over all cells
    Teuchos::RCP<matrix_t> UnCovered_p; ///< uncovered cells
    
    std::unique_ptr<mask_t> mask;       ///< interior, phys boundary, interface, covered
    std::unique_ptr<mask_t> refmask;    ///< covered (i.e. refined) or not-covered
    std::unique_ptr<mask_t> crsemask;
    
private:
    int nr_m[AMREX_SPACEDIM];           ///< number of grid points
    
    AmrIntVect_t rr_m;                  ///< refinement
    
    boundary_t bc_mp[AMREX_SPACEDIM];
};


#include "AmrMultiGridLevel.hpp"

#endif
