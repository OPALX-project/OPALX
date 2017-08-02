#include "AmrTrilinos.h"

#include "Ippl.h"

#include <Epetra_LinearProblem.h>

#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>

#include <BelosBlockCGSolMgr.hpp>


#include <AMReX_MultiFabUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>

#include <AMReX_MacBndry.H>

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
                        int lfine,
                        double scale)
{
    nLevel_m = lfine - lbase + 1;
    
    for (int ntimes = 0; ntimes < 1; ++ntimes) {
    
    // solve coarsest level
    if ( lbase == 0 ) {
        
//         const BoxArray ba = findBoundaryBoxes_m(rho[0]->boxArray());
        
//         AmrField_t rho_tmp = AmrField_t(new MultiFab(ba,
//                                                  rho[0]->DistributionMap(),
//                                                  1, 0));
//         rho_tmp->setVal(0.0, 0);
//         /*
//          * MultiFab::Copy(MultiFab&       dst,
//          *                const MultiFab& src,
//          *                int             srccomp,
//          *                int             dstcomp,
//          *                int             numcomp,
//          *                int             nghost);
//         */
//         MultiFab::Copy(*tmp, *rho[0], 0, 0, 1, 0);
//         rho_tmp->copy(*rho[0], 0, 0, 1);
//         
//         
//         
//         AmrField_t phi_tmp = AmrField_t(new MultiFab(ba,
//                                                      phi[0]->DistributionMap(),
//                                                      1, 0));
//         
//         
//         nPoints_m = getDimensions_m(ba);
//         
//         phi_tmp->setVal(0.0, 0);
//         
//         solveZeroLevel_m(phi_tmp, rho_tmp, geom[0]);
        
        phi[0]->setVal(0.0, 1);
        
        solveZeroLevel_m(phi[0], rho[0], geom[0], scale);
        
        
        /*
         * void copy (const FabArray<FAB>& src,
         *            int                  src_comp,
         *            int                  dest_comp,
         *            int                  num_comp,
         *            int                  src_nghost,
         *            int                  dst_nghost,
         * ...
         */
//         phi[0]->copy(*phi_tmp, 0, 0, 1, 0, 0);
        
    }
    
    for (int i = 1; i < nLevel_m; ++i) {
        int lev = i + lbase;
        
        phi[lev]->setVal(0.0, 1);
        
        interpFromCoarseLevel_m(phi, geom, lev-1);
        
        std::cout << "after: " << phi[lev]->min(0) << " "
                  << phi[lev]->max(0) << " " << phi[lev]->sum(0) << std::endl;
        
        const BoxArray ba = findBoundaryBoxes_m(phi[lev]->boxArray(), phi[lev]);
        
        for (int j = 0; j < 3; ++j)
            nPoints_m[j] = geom[i].Domain().length(j);
        
        build_m(rho[lev], phi[lev], geom[lev], lev);
        
        fillMatrix_m(geom[lev], phi[lev], lev);
        solve_m(phi[lev]);
        
//         phi[lev]->copy(*phi_tmp, 0, 0, 1);
    }
    
//     IntVect rr(2, 2, 2);
//     
//     for (int lev = lfine - 1; lev >= lbase; --lev) {
//         amrex::average_down(*phi[lev+1], 
//                              *phi[lev], 0, 1, rr);
//     }
    
    }
}


void AmrTrilinos::solveZeroLevel_m(AmrField_t& phi,
                                   const AmrField_t& rho,
                                   const Geometry& geom, double scale)
{
    for (int j = 0; j < 3; ++j)
        nPoints_m[j] = geom.Domain().length(j);
        
    build_m(rho, phi, geom, 0);
    fillMatrix_m(geom, phi, 0);
    solve_m(phi);
}



void AmrTrilinos::solve_m(AmrField_t& phi) {
    /*
     * unknowns
     */
    
    
    
    /*
     * solve linear system Ax = b
     */
    Teuchos::RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator> > problem =
        Teuchos::rcp( new Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>(A_mp, x_mp, rho_mp) );
        
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


void AmrTrilinos::build_m(const AmrField_t& rho, const AmrField_t& phi, const Geometry& geom, int lev)
{
    if ( Ippl::myNode() == 0 ) {
        std::cout << "nx = " << nPoints_m[0] << std::endl
                  << "ny = " << nPoints_m[1] << std::endl
                  << "nz = " << nPoints_m[2] << std::endl;
    }
    
    indexMap_m.clear();

    int localNumElements = 0;
    std::vector<double> values, xvalues;
    std::vector<int> globalindices;
    
    const double* dx = geom.CellSize();
    
    for (MFIter mfi(*rho, false); mfi.isValid(); ++mfi) {
        const Box&          bx  = mfi.validbox();
        const FArrayBox&    fab = (*rho)[mfi];
        const FArrayBox&    pot = (*phi)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    
                    
                    IntVect iv(i, j, k);
                    
                    double val = 0.0;
                    
                    IntVect xl(iv[0] - 1, iv[1], iv[2]);
                    if ( isBoundary_m(xl) ) {
                        std::cout << "xl = " << xl << " " << pot(xl) <<std::endl;
                        val += pot(xl) / ( dx[0] * dx[0] );
                    }
                        
                    IntVect xr(iv[0] + 1, iv[1], iv[2]);
                    if ( isBoundary_m(xr) ) {
                        std::cout << "xr = " << xr << " " << pot(xr) <<std::endl;
                        val += pot(xr) / ( dx[0] * dx[0] );
                    }
                    
                    IntVect yl(iv[0], iv[1] - 1, iv[2]);
                    if ( isBoundary_m(yl) ) {
                        std::cout << "yl = " << yl << " " << pot(yl) <<std::endl;
                        val += pot(yl) / ( dx[1] * dx[1] );
                    }
                    
                    IntVect yr(iv[0], iv[1] + 1, iv[2]);
                    if ( isBoundary_m(yr) ) {
                        std::cout << "yr = " << yr << " " << pot(yr) <<std::endl;
                        val += pot(yr) / ( dx[1] * dx[1] );
                    }
                    
                    IntVect zl(iv[0], iv[1], iv[2] - 1);
                    if ( isBoundary_m(zl) ) {
                        std::cout << "zl = " << zl << " " << pot(zl) << std::endl;
                        val += pot(zl) / ( dx[2] * dx[2] );
                    }
                    
                    IntVect zr(iv[0], iv[1], iv[2] + 1);
                    if ( isBoundary_m(zr) ) {
                        std::cout << "zr = " << zr << " " << pot(zr) << std::endl;
                        val += pot(zr) / ( dx[2] * dx[2] );
                    }
                    int globalidx = serialize_m(iv/*shift*/);
                    globalindices.push_back(globalidx);
                    
                    indexMap_m[globalidx] = iv;
                    
                    values.push_back(fab(iv) + val);
                    
                    xvalues.push_back(pot(iv));
                    
                    ++localNumElements;
                }
            }
        }
    }
    
    std::cout << indexMap_m.begin()->first << " " << indexMap_m.end()->first << std::endl;
    
    // compute map based on localelements
    // create map that specifies which processor gets which data
    const int baseIndex = 0;    // where to start indexing
        
    // numGlobalElements == N
    const BoxArray& grids = rho->boxArray();
    int N = grids.numPts();
    
    if ( Ippl::myNode() == 0 )
        std::cout << "N = " << N << std::endl;
    
    std::cout << localNumElements << " " << globalindices.size() << std::endl;
    
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
    
    
    // each processor fill in its part of the right-hand side
    x_mp = Teuchos::rcp( new Epetra_Vector(*map_mp, false) ); // no initialization to zero
    success = x_mp->ReplaceGlobalValues(localNumElements,
                                              &xvalues[0],
                                              &globalindices[0]);
    
    if ( success == 1 )
        throw std::runtime_error("Error in filling the vector!");
    
//     std::cout << "Vector done." << std::endl;
}


void AmrTrilinos::fillMatrix_m(const Geometry& geom, const AmrField_t& phi, int lev) {
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
        if ( isInside_m(xl, lev) ) {
//             IntVect shift(xl[0]+1, xl[1]+1, xl[2]+1);
            int gidx = serialize_m(xl/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
            
//             setBoundaryValue_m(xl, values[numEntries]);
            
            ++numEntries;
        } //else if ( isBoundary_m(xl) ) {
            
// //             double val = rho_mp->Values()[i] + bc_m[xl] / ( dx[0] * dx[0] );
//             int gidx = serialize_m(iv/*shift*/);
//             rho_mp->ReplaceGlobalValues(rho_mp->MyLength(),
//                                               &val,
//                                               &gidx);
//             rho_mp->Values()[i] += bc_m[xl] / ( dx[0] * dx[0] );
//         }
        
        // check right neighbor in x-direction
        IntVect xr(iv[0] + 1, iv[1], iv[2]);
        if ( isInside_m(xr, lev) ) {
//             IntVect shift(xr[0]+1, xr[1]+1, xr[2]+1);
            int gidx = serialize_m(xr/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
            
//             setBoundaryValue_m(xr, values[numEntries]);
            
            ++numEntries;
        } else if ( isBoundary_m(xr) ) {
//             rho_mp->Values()[i] += bc_m[xr] / ( dx[0] * dx[0] );
        }
        
        // check lower neighbor in y-direction
        IntVect yl(iv[0], iv[1] - 1, iv[2]);
        if ( isInside_m(yl, lev) ) {
//             IntVect shift(yl[0]+1, yl[1]+1, yl[2]+1);
            int gidx = serialize_m(yl/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
            
//             setBoundaryValue_m(yl, values[numEntries]);
            
            ++numEntries;
        } else if ( isBoundary_m(yl) ) {
//             rho_mp->Values()[i] += bc_m[yl] / ( dx[1] * dx[1] );
        }
        
        // check upper neighbor in y-direction
        IntVect yr(iv[0], iv[1] + 1, iv[2]);
        if ( isInside_m(yr, lev) ) {
//             IntVect shift(yr[0]+1, yr[1]+1, yr[2]+1);
            int gidx = serialize_m(yr/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
            
//             setBoundaryValue_m(yr, values[numEntries]);
            
            ++numEntries;
        } else if ( isBoundary_m(yr) ) {
//             rho_mp->Values()[i] += bc_m[yr] / ( dx[1] * dx[1] );
        }
        
        // check front neighbor in z-direction
        IntVect zl(iv[0], iv[1], iv[2] - 1);
        if ( isInside_m(zl, lev) ) {
//             IntVect shift(zl[0]+1, zl[1]+1, zl[2]+1);
            int gidx = serialize_m(zl/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
            
//             setBoundaryValue_m(zl, values[numEntries]);
            
            ++numEntries;
        } else if ( isBoundary_m(zl) ) {
//             rho_mp->Values()[i] += bc_m[zl] / ( dx[2] * dx[2] );
        }
        
        // check back neighbor in z-direction
        IntVect zr(iv[0], iv[1], iv[2] + 1);
        if ( isInside_m(zr, lev) ) {
//             IntVect shift(zr[0]+1, zr[1]+1, zr[2]+1);
            int gidx = serialize_m(zr/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
            
//             setBoundaryValue_m(zr, values[numEntries]);
            
            ++numEntries;
        } else if ( isBoundary_m(zr) ) {
//             rho_mp->Values()[i] += bc_m[zr] / ( dx[2] * dx[2] );
        }
        
        // check center
        if ( isInside_m(iv, lev) ) {
            indices[numEntries] = globalRow;
            values[numEntries]  = 2.0 / ( dx[0] * dx[0] ) +
                               2.0 / ( dx[1] * dx[1] ) +
                               2.0 / ( dx[2] * dx[2] );
            
//             setBoundaryValue_m(iv, values[numEntries]);                   
            
            ++numEntries;
        }
        
//         std::cout << iv << " " << numEntries << " " << globalRow << std::endl; std::cin.get();
        
        int error = A_mp->InsertGlobalValues(globalRow, numEntries, &values[0], &indices[0]);
        
        if ( error != 0 )
            throw std::runtime_error("Error in filling the matrix!");
    }
    
    A_mp->FillComplete(true);
    
//     std::cout << "Matrix done." << std::endl;
    
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


inline bool AmrTrilinos::isInside_m(const IntVect& iv, int lev) const {
    
    return ( iv[0] > -1/*-2*/ && iv[0] < nPoints_m[0]/*-1*/ &&
             iv[1] > -1/*-2*/ && iv[1] < nPoints_m[1]/*-1*/ &&
             iv[2] > -1/*-2*/ && iv[2] < nPoints_m[2]/*-1*/);
}


inline bool AmrTrilinos::isBoundary_m(const IntVect& iv) const {
    return ( bc_m.find(iv) != bc_m.end() );
}


BoxArray AmrTrilinos::findBoundaryBoxes_m(const BoxArray& ba, AmrField_t& phi)
{
    // ghosts
    BoxList boundary_cells = amrex::GetBndryCells(ba, 1 /*ngrow, nghosts*/);
    
    boundary_m.define(boundary_cells);

    if (boundary_m.size() == 0)
        throw std::runtime_error("No boundary!");    
    
    BoxArray cba = ba;
    
    cba.resize( ba.size() + boundary_m.size() );
    
    for (uint i = 0; i < boundary_m.size(); ++i) {
        cba.set(i + ba.size(), boundary_m[i]);
    }
    
    cba.removeOverlap(true);
    
    bc_m.clear();
    
//     phi->FillBoundary(false);
//     
//     for (MFIter mfi(boundary_m, phi->DistributionMap(), false); mfi.isValid(); ++mfi) {
//         const Box&          bx  = mfi.fabbox();
//         const FArrayBox&    fab = (*phi)[mfi];
// 
//         const int* lo = bx.loVect();
//         const int* hi = bx.hiVect();
//         
//         for (int i = lo[0]; i <= hi[0]; ++i) {
//             for (int j = lo[1]; j <= hi[1]; ++j) {
//                 for (int k = lo[2]; k <= hi[2]; ++k) {                    
//                     IntVect iv(i, j, k);
//                     
//                     std::cout << iv << " " << fab(iv) << std::endl; std::cin.get();
//                 }
//             }
//         }
//     }
    
    for (uint b = 0; b < boundary_m.size(); ++b) {
        const Box& bx = boundary_m[b];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    bc_m.insert(IntVect(i, j, k));
                    
                    
                    
// //                     std::cout << IntVect(i, j, k) << std::endl; std::cin.get();
                }
            }
        }
    }
// //     
//     std::cerr << "done." << std::endl; std::cin.get();
    
    return cba;
    
//     DistributionMapping ndm{nba};
    
//     MultiFab mf(nba, ndm, 1, 0);
}


void AmrTrilinos::applyBoundary_m(const AmrField_t& phi) {
/*    
    
    bc_m.clear();
    
    for (MFIter mfi(boundary_m, phi->DistributionMap()); mfi.isValid(); ++mfi) {
        const Box&          bx  = mfi.validbox();
        const FArrayBox&    fab = (*phi)[mfi];

        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {                    
                    IntVect iv(i, j, k);
                    
                    IntVect shift(i+1, j+1, k+1);
                    
                    int globalidx = serialize_m(shift);
                    
                    std::cout << iv << " " << globalidx << std::endl; //std::cin.get();
                    
                    bc_m[globalidx] = fab(iv);
                }
            }
        }
    }*/
}


// void AmrTrilinos::setBoundaryValue_m(const IntVect& iv, double& value) {
//     IntVect shift(iv[0]+1, iv[1]+1, iv[2]+1);
//     int gidx = serialize_m(shift);
//     
//     value = (bc_m.find(gidx) != bc_m.end() ) ? bc_m[gidx] * value : value; 
// }


IntVect AmrTrilinos::getDimensions_m(const BoxArray& ba) {
    Box bx = ba.minimalBox();
    std::cout << "minbox = " << bx << std::endl;
    return IntVect(bx.length(0), bx.length(1), bx.length(2));
}

void AmrTrilinos::copyBack_m(AmrField_t& phi,
                             const Teuchos::RCP<Epetra_MultiVector>& sol)
{
    int localidx = 0;
    for (MFIter mfi(*phi, false); mfi.isValid(); ++mfi) {
        const Box&          bx  = mfi.validbox();
        FArrayBox&          fab = (*phi)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    IntVect iv(i, j, k);
                    
                    fab(iv) = (*sol)[0][localidx++];
                }
            }
        }
    }
}

class NoOpPhysBC
    : public amrex::PhysBCFunctBase
{
public:
    NoOpPhysBC () {}
    virtual ~NoOpPhysBC () {}
    virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override { }
    using amrex::PhysBCFunctBase::FillBoundary;
};


void AmrTrilinos::interpFromCoarseLevel_m(container_t& phi,
                                          const Array<Geometry>& geom,
                                          int lev)
{
    NoOpPhysBC cphysbc, fphysbc;
// //     PhysBCFunct cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
    /*PCInterp*/CellConservativeLinear mapper;
    
    IntVect rr(2, 2, 2);
    
    phi[lev+1]->setVal(0.0, 1);
    
    amrex::InterpFromCoarseLevel(*phi[lev+1], 0.0,
                                 *phi[lev],
                                 0, 0, 1,
                                 geom[lev], geom[lev+1],
                                 cphysbc, fphysbc,
                                 rr, &mapper, bcs);
    
    for (uint lev = 0; lev < phi.size(); ++lev) {
        const Geometry& gm = geom[lev];
        phi[lev]->FillBoundary(gm.periodicity());
    }
    
    std::cout << "BEFORE SETTING THE BC VALUES" << std::endl;
    for (MFIter mfi(*phi[lev+1], false); mfi.isValid(); ++mfi) {
        const Box&          bx  = mfi.validbox();
        const FArrayBox&    fab = (*phi[lev+1])[mfi];

        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]-1; i <= hi[0]+1; ++i) {
            for (int j = lo[1]-1; j <= hi[1]+1; ++j) {
                for (int k = lo[2]-1; k <= hi[2]+1; ++k) {                    
                    IntVect iv(i, j, k);
                    
                    if ( iv[0] == 0 && iv[1] == 8 && iv[2] == 32 ) {
                    std::cout << iv << " " << fab(iv) << std::endl;  std::cin.get();
                    }
                }
            }
        }
    }
//     
//     std::cout << "AFTER SETTING THE BC VALUES" << std::endl;
//     for (MFIter mfi(*phi[lev+1], false); mfi.isValid(); ++mfi) {
//         const Box&          bx  = mfi.validbox();
//         const FArrayBox&    fab = (*phi[lev+1])[mfi];
// 
//         const int* lo = bx.loVect();
//         const int* hi = bx.hiVect();
//         
//         for (int i = lo[0]-1; i <= hi[0]+1; ++i) {
//             for (int j = lo[1]-1; j <= hi[1]+1; ++j) {
//                 for (int k = lo[2]-1; k <= hi[2]+1; ++k) {                    
//                     IntVect iv(i, j, k);
//                     
//                     std::cout << iv << " " << fab(iv) << std::endl; std::cin.get();
//                 }
//             }
//         }
//     }
}
