#ifndef AMR_MULTI_GRID_LEVEL
#define AMR_MULTI_GRID_LEVEL

#include <AMReX_IntVect.H>

#include "Ippl.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Epetra_MpiComm.h>

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
    
    AmrMultiGridLevel(const AmrField_u& _rho,
                      const AmrField_u& _phi,
                      const AmrGeometry_t& _geom,
                      Epetra_MpiComm& comm);
    
    ~AmrMultiGridLevel();
    
    
    Teuchos::RCP<matrix_t> A_p;
    Teuchos::RCP<vector_t> rho_p;
    Teuchos::RCP<vector_t> phi_p;
    Teuchos::RCP<vector_t> r_p;
    Teuchos::RCP<vector_t> err_p;
    
    AmrField_u error;
    AmrField_u phi_saved;
    AmrField_t* phi;
    AmrField_t* rho;
    AmrField_u residual;
    
    const amrex::BoxArray& grids;
    const amrex::DistributionMapping& dmap;
    const AmrGeometry_t& geom;
    
    void buildRHS(/*const AmrField_u rho*/);
    
    void copyBack(AmrField_t& mf, const Teuchos::RCP<vector_t>& mv);

private:
    int serialize_m(const AmrIntVect_t& iv) const;
    
    
    
    void buildLevelMask_m();
    
    void buildMap(const AmrField_u& phi, Epetra_MpiComm& comm);
    
    void buildMatrix_m();
    
private:
    Vector_t nr_m;
    std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > masks_m;
    
    Teuchos::RCP<Epetra_Map> map_mp;
};


#include "AmrMultiGridLevel.hpp"

#endif
