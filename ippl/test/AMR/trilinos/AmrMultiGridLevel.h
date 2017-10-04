#ifndef AMR_MULTI_GRID_LEVEL
#define AMR_MULTI_GRID_LEVEL

#include <vector>

#include <AMReX_IntVect.H>

#include "Ippl.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Epetra_MpiComm.h>

#include "AmrMultiGridCore.h"

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
    typedef amrex::FabArray<amrex::BaseFab<int> > mask_t;
    
    typedef std::vector<int>        indices_t;
    typedef std::vector<double>     coefficients_t;
    
    typedef Vektor<double, 3> Vector_t;
    
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
    
public:
    AmrMultiGridLevel(const amrex::BoxArray& _grids,
                      const amrex::DistributionMapping& _dmap,
                      const AmrGeometry_t& _geom,
                      const AmrIntVect_t& rr,
                      AmrBoundary<AmrMultiGridLevel<MatrixType, VectorType> >* bc,
                      Epetra_MpiComm& comm);
    
    ~AmrMultiGridLevel();
    
    int serialize(const AmrIntVect_t& iv) const;
    
    bool isBoundary(const AmrIntVect_t& iv) const;
    
    void applyBoundary(const AmrIntVect_t& iv,
                       indices_t& indices,
                       coefficients_t& values,
                       const double& value);
    
    const AmrIntVect_t& refinement() const;
    
private:
    void buildLevelMask_m();
    
    void buildMap_m(const Epetra_MpiComm& comm);
    
public:
    Teuchos::RCP<Epetra_Map> map_p;     ///< core map
    
    Teuchos::RCP<matrix_t> A_p;         ///< Poisson matrix
    Teuchos::RCP<matrix_t> R_p;         ///< restriction matrix
    Teuchos::RCP<matrix_t> I_p;         ///< interpolation matrix
    Teuchos::RCP<matrix_t> Bcrse_p;
    Teuchos::RCP<matrix_t> Bfine_p;     ///< physical boundary vector
    Teuchos::RCP<matrix_t> As_p;        ///< special Poisson matrix
    
    Teuchos::RCP<matrix_t> S_p;         ///< for Gauss-Seidel (smoothing)
    
    Teuchos::RCP<vector_t> rho_p;       ///< charge density
    Teuchos::RCP<vector_t> phi_p;       ///< potential vector
    Teuchos::RCP<vector_t> residual_p;
    Teuchos::RCP<vector_t> error_p;
    
    Teuchos::RCP<matrix_t> UnCovered_p;
    
    std::unique_ptr<mask_t> mask;       ///< interior, phys boundary, interface, covered
    
    const amrex::BoxArray& grids;
    const amrex::DistributionMapping& dmap;
    const AmrGeometry_t& geom;
    
private:
    int nr_m[BL_SPACEDIM];
    
    AmrIntVect_t rr_m;
    
    std::unique_ptr<AmrBoundary<AmrMultiGridLevel<MatrixType, VectorType> > > bc_mp;
    
};


#include "AmrMultiGridLevel.hpp"

#endif
