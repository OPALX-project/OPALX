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

using namespace amrex;

class AmrTrilinos {
    
public:
    typedef std::unique_ptr<MultiFab>   AmrField_t;
    typedef Array<AmrField_t>           container_t;
    
public:
    
    AmrTrilinos();
    
    void solve(const container_t& rho,
               container_t& phi,
               const Array<Geometry>& geom,
               int lbase,
               int lfine);
    
    
private:
    
    void solveZeroLevel_m(AmrField_t& phi,
                          const AmrField_t& rho,
                          const Geometry& geom);
    
    // single level solve
    void solve_m(AmrField_t& phi);
    
//     void solveGrid_m(
    
    // create map + right-hand side
    void build_m(const AmrField_t& rho);
    
    void fillMatrix_m(const Geometry& geom);
    
    int serialize_m(const IntVect& iv) const;
    
    IntVect deserialize_m(int idx) const;
    
    bool isInside_m(const IntVect& iv) const;
    
    
    void copyBack_m(AmrField_t& phi,
                    const Teuchos::RCP<Epetra_MultiVector>& sol);
    
    
//     void grid_m(AmrField_t& rhs, const container_t& rho,
//                 const Geometry& geom, int lev);
    
    void interpFromCoarseLevel_m(const container_t& phi,
                                 const Array<Geometry>& geom,
                                 int lev);
    
private:
    
    Epetra_MpiComm epetra_comm_m;
    
    Teuchos::RCP<Epetra_Map> map_mp;
    
    int nLevel_m;
    
    // info for a level
    int nPoints_m[3];
    
    // (x, y, z) --> idx
    std::map<int, IntVect> indexMap_m;
    
    
    Teuchos::RCP<Epetra_CrsMatrix> A_mp;
    Teuchos::RCP<Epetra_Vector> rho_mp;
    
};

#endif
