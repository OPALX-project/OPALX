#ifndef AMR_TRILINOS_H
#define AMR_TRILINOS_H

#include <AMReX.H>

#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>


#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <map>
#include <memory>

#include "EllipticDomain.h"

using namespace amrex;

class AmrTrilinos {
    
public:
    typedef std::unique_ptr<MultiFab>   AmrField_t;
    typedef Array<AmrField_t>           container_t;
    
public:
    
    AmrTrilinos();
    
    void solveWithGeometry(const container_t& rho,
               container_t& phi,
               container_t& efield,
               const Array<Geometry>& geom,
               int lbase,
               int lfine);
    
    void solve(const container_t& rho,
               container_t& phi,
               container_t& efield,
               const Array<Geometry>& geom,
               int lbase,
               int lfine,
               double scale);
    
    
private:
    
    void solveLevel_m(AmrField_t& phi,
                          const AmrField_t& rho,
                          const Geometry& geom, int lev);
    
    // single level solve
    void linSysSolve_m(AmrField_t& phi, const Geometry& geom);
    
//     void solveGrid_m(
    
    // create map + right-hand side
    void build_m(const AmrField_t& rho, const AmrField_t& phi, const Geometry& geom, int lev);
    
    void fillMatrix_m(const Geometry& geom, const AmrField_t& phi, int lev);
    
    int serialize_m(const IntVect& iv) const;
    
//     IntVect deserialize_m(int idx) const;
    
    bool isInside_m(const IntVect& iv) const;
    
//     bool isBoundary_m(const IntVect& iv) const;
    
//     void findBoundaryBoxes_m(const BoxArray& ba);
    
    void copyBack_m(AmrField_t& phi,
                    const Teuchos::RCP<Epetra_MultiVector>& sol,
                    const Geometry& geom);
    
    
//     void grid_m(AmrField_t& rhs, const container_t& rho,
//                 const Geometry& geom, int lev);
    
    void interpFromCoarseLevel_m(container_t& phi,
                                 const Array<Geometry>& geom,
                                 int lev);
    
    // efield = -grad phi
    void getGradient(AmrField_t& phi,
                     AmrField_t& efield,
                     const Geometry& geom, int lev);
    
    
    void buildLevelMask_m(const BoxArray& grids,
                          const DistributionMapping& dmap,
                          Geometry& geom,
                          int lev)
    {
        
        masks_m[lev].reset(new FabArray<BaseFab<int> >(grids, dmap, 1, 1));
        
        // covered   : ghost cells covered by valid cells of this FabArray
        //             (including periodically shifted valid cells)
        // notcovered: ghost cells not covered by valid cells
        //             (including ghost cells outside periodic boundaries)
        // physbnd   : boundary cells outside the domain (excluding periodic boundaries)
        // interior  : interior cells (i.e., valid cells)
        int covered = -1;
        int interior = 0;
        int notcovered = 1; // -> boundary
        int physbnd = 2;
        
        masks_m[lev]->BuildMask(geom.Domain(), geom.periodicity(),
                                covered, notcovered, physbnd, interior);
    }
    
    
    void buildGeometry_m(const AmrField_t& rho, const AmrField_t& phi,
                         const Geometry& geom, int lev, EllipticDomain& bp);
    
    void solveLevelGeometry_m(AmrField_t& phi,
                              const AmrField_t& rho,
                              const Geometry& geom, int lev, EllipticDomain& bp);
    
    void fillMatrixGeometry_m(const Geometry& geom, const AmrField_t& phi, int lev, EllipticDomain& bp);
    
    
private:
    
    Array<std::unique_ptr<FabArray<BaseFab<int> > > > masks_m;
    
    std::set<IntVect> bc_m;
    
    Epetra_MpiComm epetra_comm_m;
    
    Teuchos::RCP<Epetra_Map> map_mp;
    
    int nLevel_m;
    
    // info for a level
    IntVect nPoints_m;
    
    // (x, y, z) --> idx
    std::map<int, IntVect> indexMap_m;
    
    
    Teuchos::RCP<Epetra_CrsMatrix> A_mp;
    Teuchos::RCP<Epetra_Vector> rho_mp;
    Teuchos::RCP<Epetra_Vector> x_mp;
    
};

#endif
