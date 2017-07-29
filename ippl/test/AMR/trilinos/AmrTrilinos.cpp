#include "AmrTrilinos.h"

#include "Ippl.h"

#include <Epetra_LinearProblem.h>

#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>

#include <BelosBlockCGSolMgr.hpp>


#include <AMReX_MultiFabUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>

// #include "EpetraExt_RowMatrixOut.h"
// #include "EpetraExt_MultiVectorOut.h"
// #include "EpetraExt_VectorOut.h"


AmrTrilinos::AmrTrilinos()
    : epetra_comm_m(Ippl::getComm()),
      A_mp(0),
      rho_mp(0),
      map_mp(0)
{ }



void AmrTrilinos::solve(const container_t& rho,
                        container_t& phi,
                        const Array<Geometry>& geom,
                        int lbase,
                        int lfine)
{
    nLevel_m = lfine - lbase + 1;
    
    
    // solve coarsest level
    if ( lbase == 0 )
        solveZeroLevel_m(phi[0], rho[0], geom[0]);
    
    for (int i = 1; i < nLevel_m; ++i) {
        int lev = i + lbase;
        
        interpFromCoarseLevel_m(phi, geom, lev);
        
        for (int j = 0; j < 3; ++j)
            nPoints_m[j] = geom[lev].Domain().length(j);
        
        buildMap_m(rho[lev]);
        
        fillMatrix_m(geom[lev]);
        solve_m(phi[lev]);
    }
    
    IntVect rr(2, 2, 2);
    
    for (int lev = lfine - 1; lev >= lbase; --lev) {
        amrex::average_down(*phi[lev+1], 
                             *phi[lev], 0, 1, rr);
    }
    
//         }
//     }


//         
//         
//         
//         AmrField_t rhs;
//         
//         
//         for (int j = 0; j < 3; ++j)
//             nPoints_m[j] = geom[i].Domain().length(j);
//         
//         
//         if ( i > 0 )
//             break;
//         
//         build_m(rhs);
//         fillMatrix_m(geom[i]);
//         solve_m(phi[i]);
//     }
}


void AmrTrilinos::buildMap_m(const AmrField_t& phi)
{
    build_m(phi);
    
    /*
    int localNumElements = 0;
    std::vector<double> values;
    std::vector<int> globalindices;
    
    indexMap_m.clear();
    
    IntVect rr(2, 2, 2);
    
    for (MFIter mfi(*phi, false); mfi.isValid(); ++mfi) {
        const Box&          bx  = mfi.validbox();
        const FArrayBox&    fab = (*phi)[mfi];
                
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    IntVect iv(i, j, k);
                    
                    int globalidx = serialize_m(iv);
                    globalindices.push_back(globalidx);
                    
                    indexMap_m[globalidx] = iv;
                    
                    values.push_back(fab(iv));
                    
                    ++localNumElements;
                    
//                     std::cout << fab(iv) << std::endl;
                }
            }
        }
    }
    
    const BoxArray& grids = phi->boxArray();
    int N = grids.numPts();
    
    std::cout << "N = " << N << std::endl;*/
}


void AmrTrilinos::solveZeroLevel_m(AmrField_t& phi,
                                   const AmrField_t& rho,
                                   const Geometry& geom)
{
    for (int j = 0; j < 3; ++j)
        nPoints_m[j] = geom.Domain().length(j);
        
    build_m(rho);
    fillMatrix_m(geom);
    solve_m(phi);
}



void AmrTrilinos::solve_m(AmrField_t& phi) {
    /*
     * unknowns
     */
    // each processor fill in its part of the right-hand side
    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp( new Epetra_Vector(*map_mp, false) ); // no initialization to zero
    x->PutScalar(0.0);
    
    
    
    /*
     * solve linear system Ax = b
     */
    Teuchos::RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator> > problem =
        Teuchos::rcp( new Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>(A_mp, x, rho_mp) );
        
    if ( !problem->setProblem() )
        throw std::runtime_error("Belos::LinearProblem failed to set up correctly!");
    
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp( new Teuchos::ParameterList );
    
//     params->set( "Block Size", 4 );
    params->set( "Maximum Iterations", 10000 );
    params->set("Convergence Tolerance", 1.0e-2);
    
    Belos::BlockCGSolMgr<double, Epetra_MultiVector, Epetra_Operator> solver(problem, params);
    
    Belos::ReturnType ret = solver.solve();
    
    // get the solution from the problem
    if ( ret == Belos::Converged ) {
        Teuchos::RCP<Epetra_MultiVector> sol = problem->getLHS();
        
        copyBack_m(phi, sol);
        
        // print stuff
        if ( epetra_comm_m.MyPID() == 0 ) {
            std::cout << "Achieved tolerance: " << solver.achievedTol() << std::endl
                      << "Number of iterations: " << solver.getNumIters() << std::endl;
        }
        
    } else {
        if ( epetra_comm_m.MyPID() == 0 ) {
            std::cout << "Not converged. Achieved tolerance after " << solver.getNumIters() << " iterations is "
                      << solver.achievedTol() << "." << std::endl;
        }
    }
}


void AmrTrilinos::build_m(const AmrField_t& rho)
{
    if ( Ippl::myNode() == 0 ) {
        std::cout << "nx = " << nPoints_m[0] << std::endl
                  << "ny = " << nPoints_m[1] << std::endl
                  << "nz = " << nPoints_m[2] << std::endl;
    }
    
    int localNumElements = 0;
    std::vector<double> values;
    std::vector<int> globalindices;
    
    for (MFIter mfi(*rho, false); mfi.isValid(); ++mfi) {
        const Box&          bx  = mfi.validbox();
        const FArrayBox&    fab = (*rho)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    IntVect iv(i, j, k);
                    
                    int globalidx = serialize_m(iv);
                    globalindices.push_back(globalidx);
                    
                    indexMap_m[globalidx] = iv;
                    
                    values.push_back(fab(iv));
                    
                    ++localNumElements;
                }
            }
        }
    }
    
    // compute map based on localelements
    // create map that specifies which processor gets which data
    const int baseIndex = 0;    // where to start indexing
        
    // numGlobalElements == N
    const BoxArray& grids = rho->boxArray();
    int N = grids.numPts();
    
    if ( Ippl::myNode() == 0 )
        std::cout << "N = " << N << std::endl;
    
    map_mp = Teuchos::rcp( new Epetra_Map(N, localNumElements,
                                         &globalindices[0], baseIndex, epetra_comm_m) );
    
    // each processor fill in its part of the right-hand side
    rho_mp = Teuchos::rcp( new Epetra_Vector(*map_mp, false) );
        
    // fill vector
    int success = rho_mp->ReplaceGlobalValues(localNumElements,
                                              &values[0],
                                              &globalindices[0]);
    
    if ( success == 1 )
        throw std::runtime_error("Error in filling the vector!");
}


void AmrTrilinos::fillMatrix_m(const Geometry& geom) {
    // 3D --> 7 elements per row
    A_mp = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map_mp, 7) );
    
    int indices[7];
    double values[7];
    
    const double* dx = geom.CellSize();
    
    int * myGlobalElements = map_mp->MyGlobalElements();
    for ( int i = 0; i < map_mp->NumMyElements(); ++i) {
        /*
         * GlobalRow	- (In) Row number (in global coordinates) to put elements.
         * NumEntries	- (In) Number of entries.
         * Values	- (In) Values to enter.
         * Indices	- (In) Global column indices corresponding to values.
         */
        int globalRow = myGlobalElements[i];
        int numEntries = 0;
        
        IntVect iv = indexMap_m[globalRow]/*deserialize_m(globalRow)*/;
        
        // check left neighbor in x-direction
        IntVect xl(iv[0] - 1, iv[1], iv[2]);
        int gidx = serialize_m(xl);
        if ( isInside_m(xl) ) {
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
            ++numEntries;
        }
        
        // check right neighbor in x-direction
        IntVect xr(iv[0] + 1, iv[1], iv[2]);
        gidx = serialize_m(xr);
        if ( isInside_m(xr) ) {
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
            ++numEntries;
        }
        
        // check lower neighbor in y-direction
        IntVect yl(iv[0], iv[1] - 1, iv[2]);
        gidx = serialize_m(yl);
        if ( isInside_m(yl) ) {
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
            ++numEntries;
        }
        
        // check upper neighbor in y-direction
        IntVect yr(iv[0], iv[1] + 1, iv[2]);
        gidx = serialize_m(yr);
        if ( isInside_m(yr) ) {
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
            ++numEntries;
        }
        
        // check front neighbor in z-direction
        IntVect zl(iv[0], iv[1], iv[2] - 1);
        gidx = serialize_m(zl);
        if ( isInside_m(zl) ) {
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
            ++numEntries;
        }
        
        // check back neighbor in z-direction
        IntVect zr(iv[0], iv[1], iv[2] + 1);
        gidx = serialize_m(zr);
        if ( isInside_m(zr) ) {
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
            ++numEntries;
        }
        
        // check center
        gidx = globalRow;
        if ( isInside_m(iv) ) {
            indices[numEntries] = globalRow;
            values[numEntries]  = 2.0 / ( dx[0] * dx[0] ) +
                               2.0 / ( dx[1] * dx[1] ) +
                               2.0 / ( dx[2] * dx[2] );
            
            ++numEntries;
        }
        
        int error = A_mp->InsertGlobalValues(globalRow, numEntries, &values[0], &indices[0]);
        
        if ( error != 0 )
            throw std::runtime_error("Error in filling the matrix!");
    }
    
    A_mp->FillComplete(true);
    
//     A_mp->OptimizeStorage();
    
    /*
     * some printing
     */
    
//     EpetraExt::RowMatrixToMatlabFile("matrix.txt", *A_mp);
    
    if ( Ippl::myNode() == 0 ) {
        std::cout << "Global info" << std::endl
                  << "Number of rows:      " << A_mp->NumGlobalRows() << std::endl
                  << "Number of cols:      " << A_mp->NumGlobalCols() << std::endl
                  << "Number of diagonals: " << A_mp->NumGlobalDiagonals() << std::endl
                  << "Number of non-zeros: " << A_mp->NumGlobalNonzeros() << std::endl
                  << std::endl;
    }
    
    Ippl::Comm->barrier();
    
    for (int i = 0; i < Ippl::getNodes(); ++i) {
        
        if ( i == Ippl::myNode() ) {
            std::cout << "Rank:                "
                      << i << std::endl
                      << "Number of rows:      "
                      << A_mp->NumMyRows() << std::endl
                      << "Number of cols:      "
                      << A_mp->NumMyCols() << std::endl
                      << "Number of diagonals: "
                      << A_mp->NumMyDiagonals() << std::endl
                      << "Number of non-zeros: "
                      << A_mp->NumMyNonzeros() << std::endl
                      << std::endl;
        }
        Ippl::Comm->barrier();
    }
}



inline int AmrTrilinos::serialize_m(const IntVect& iv) const {
    return iv[0] + nPoints_m[0] * iv[1] + nPoints_m[0] * nPoints_m[1] * iv[2];
}


inline IntVect AmrTrilinos::deserialize_m(int idx) const {
    int i = idx % nPoints_m[0];
    int j = ( idx / nPoints_m[0] ) % nPoints_m[1];
    int k = ( idx / ( nPoints_m[0] * nPoints_m[1] ) ) % nPoints_m[2];
    
    return IntVect(i, j, k);
}


inline bool AmrTrilinos::isInside_m(const IntVect& iv) const {
    
    return ( iv[0] > -1 && iv[0] < nPoints_m[0] &&
             iv[1] > -1 && iv[1] < nPoints_m[1] &&
             iv[2] > -1 && iv[2] < nPoints_m[2]);
}


void AmrTrilinos::copyBack_m(AmrField_t& phi,
                             const Teuchos::RCP<Epetra_MultiVector>& sol)
{
    int idx = 0;
    for (MFIter mfi(*phi, false); mfi.isValid(); ++mfi) {
        const Box&          bx  = mfi.validbox();
        FArrayBox&          fab = (*phi)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    IntVect iv(i, j, k);
                    
//                     int globalidx = serialize_m(iv);
                    
//                     std::cout << iv << " " << globalidx << std::endl; std::cin.get();
                    
                    fab(iv) = (*sol)[0][idx++/*globalidx*/];
                }
            }
        }
    }
}


void AmrTrilinos::interpFromCoarseLevel_m(const container_t& phi,
                                          const Array<Geometry>& geom,
                                          int lev)
{
    PhysBCFunct cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
    PCInterp mapper;
    
    IntVect rr(2, 2, 2);
    
    const BoxArray& ba = phi[lev]->boxArray();
    const DistributionMapping& dm = phi[lev]->DistributionMap();
    AmrField_t fphi = AmrField_t(new MultiFab(ba, dm, 1, 0));
    fphi->setVal(0.0);
    
    amrex::InterpFromCoarseLevel(*fphi, 0.0,
                                 *phi[lev-1],
                                 0, 0, 1,
                                 geom[lev-1], geom[lev],
                                 cphysbc, fphysbc,
                                 rr, &mapper, bcs);
    
    MultiFab::Copy(*phi[lev], *fphi, 0, 0, 1, 0);
}

// void AmrTrilinos::grid_m(AmrField_t& rhs,
//                          const container_t& rho,
//                          const Geometry& geom,
//                          int lev)
// {
//     if ( lev == 0 ) {
//         rhs.reset( new MultiFab(rho[lev]->boxArray(),
//                                 rho[lev]->DistributionMap(),
//                                 rho[lev]->nComp(),
//                                 rho[lev]->nGrow()));
//         MultiFab::Copy(*rhs, *rho[lev], 0, 0,
//                        rho[lev]->nComp(),
//                        rho[lev]->nGrow());
//         
//         return;
//     }
//     
//     Box bx = geom.Domain();
//     
//     BoxArray
//     
//     rhs.reset( new MultiFab(
//     
//     std::cout << bx << std::endl;
//     
// }
