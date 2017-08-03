#include "AmrTrilinos.h"

#include "Ippl.h"

#include <Epetra_LinearProblem.h>

#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>

#include <BelosBlockCGSolMgr.hpp>


#include <AMReX_MultiFabUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>

#include <AMReX_BCUtil.H>


// #include "EpetraExt_RowMatrixOut.h"
// #include "EpetraExt_MultiVectorOut.h"
// #include "EpetraExt_VectorOut.h"


class NoOpPhysBC
    : public amrex::PhysBCFunctBase
{
public:
    NoOpPhysBC () {}
    virtual ~NoOpPhysBC () {}
    virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override { }
    using amrex::PhysBCFunctBase::FillBoundary;
};

AmrTrilinos::AmrTrilinos()
    : epetra_comm_m(Ippl::getComm()),
      A_mp(0),
      rho_mp(0),
      map_mp(0)
{ }



void AmrTrilinos::solve(const container_t& rho,
                        container_t& phi,
                        container_t& efield,
                        const Array<Geometry>& geom,
                        int lbase,
                        int lfine,
                        double scale)
{
    nLevel_m = lfine - lbase + 1;
    
    // solve coarsest level
    if ( lbase == 0 ) {
        
        phi[0]->setVal(0.0, 1);
        
        solveZeroLevel_m(phi[0], rho[0], geom[0], scale);
        
//         Array<BCRec> bc( phi[0]->nComp());
//         for (int n = 0; n < phi[0]->nComp(); ++n) {
//             for (int idim = 0; idim < 3; ++idim) {
// //                 if (Geometry::isPeriodic(idim)) {
//                     bc[n].setLo(idim, BCType::int_dir); // interior
//                     bc[n].setHi(idim, BCType::int_dir);
// //                 }
// //                 else {
// //                     bc[n].setLo(idim, BCType::foextrap); //first-order extrapolation
// //                     bc[n].setHi(idim, BCType::foextrap);
// //                 }
//             }
//         }
// 
//         amrex::FillDomainBoundary(*phi[0], geom[0], bc);
        
        if ( phi[0]->contains_nan() )
            throw std::runtime_error("Nan");
        
        
    }
    
    for (int i = 1; i < nLevel_m; ++i) {
        int lev = i + lbase;
        
        phi[lev]->setVal(0.0, 1);
        
        interpFromCoarseLevel_m(phi, geom, lev-1);

        
        if ( phi[lev]->contains_nan() )
            throw std::runtime_error("Nan");
        
        findBoundaryBoxes_m(phi[lev]->boxArray());
        
        for (int j = 0; j < 3; ++j)
            nPoints_m[j] = geom[i].Domain().length(j);
        
        build_m(rho[lev], phi[lev], geom[lev], lev);
        
        fillMatrix_m(geom[lev], phi[lev], lev);
        solve_m(phi[lev], geom[lev]);
        
//         Array<BCRec> bc( phi[lev]->nComp());
//         for (int n = 0; n < phi[lev]->nComp(); ++n) {
//             for (int idim = 0; idim < 3; ++idim) {
// //                 if (Geometry::isPeriodic(idim)) {
//                     bc[n].setLo(idim, BCType::int_dir); // interior
//                     bc[n].setHi(idim, BCType::int_dir);
// //                 }
// //                 else {
// //                     bc[n].setLo(idim, BCType::foextrap); //first-order extrapolation
// //                     bc[n].setHi(idim, BCType::foextrap);
// //                 }
//             }
//         }
//         
//         amrex::FillDomainBoundary(*phi[lev], geom[lev], bc);
        
    }
    
    IntVect rr(2, 2, 2);
    
    for (int lev = lfine - 1; lev >= lbase; --lev) {
        amrex::average_down(*phi[lev+1], 
                             *phi[lev], geom[lev+1], geom[lev], 0, 1, rr);
    }
    
    getGradient(phi, efield, geom);
    
    
    
    for (int lev = lfine - 1; lev >= lbase; --lev) {
        amrex::average_down(*efield[lev+1],
                            *efield[lev],
                            geom[lev+1],
                            geom[lev],
                            0, 3, rr);
    }
    
    
    
//     std::unique_ptr<FluxRegister> flux_reg;
//     for (int lev = lfine - 1; lev >= lbase; --lev) {
//         flux_reg.reset(new FluxRegister(efield[lev+1]->boxArray(),
//                                         efield[lev+1]->DistributionMap(),
//                                         rr,
//                                         lev+1, 3));
//         
//         // update lev based on coarse-fine flux mismatch
//         flux_reg->Reflux(*efield[lev], 1.0, 0, 0, efield[lev]->nComp(), geom[lev]);
//     }
}


void AmrTrilinos::solveZeroLevel_m(AmrField_t& phi,
                                   const AmrField_t& rho,
                                   const Geometry& geom, double scale)
{
    for (int j = 0; j < 3; ++j)
        nPoints_m[j] = geom.Domain().length(j);
    
    build_m(rho, phi, geom, 0);
    fillMatrix_m(geom, phi, 0);
    solve_m(phi, geom);
}



void AmrTrilinos::solve_m(AmrField_t& phi, const Geometry& geom) {
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
    params->set("Convergence Tolerance", 1.0e-12);
    
    Belos::BlockCGSolMgr<double, Epetra_MultiVector, Epetra_Operator> solver(problem, params);
    
    Belos::ReturnType ret = solver.solve();
    
    // get the solution from the problem
    if ( ret == Belos::Converged ) {
        Teuchos::RCP<Epetra_MultiVector> sol = problem->getLHS();
        
        copyBack_m(phi, sol, geom);
        
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
                        val += pot(xl) / ( dx[0] * dx[0] );
                    }
                        
                    IntVect xr(iv[0] + 1, iv[1], iv[2]);
                    if ( isBoundary_m(xr) ) {
                        val += pot(xr) / ( dx[0] * dx[0] );
                    }
                    
                    IntVect yl(iv[0], iv[1] - 1, iv[2]);
                    if ( isBoundary_m(yl) ) {
                        val += pot(yl) / ( dx[1] * dx[1] );
                    }
                    
                    IntVect yr(iv[0], iv[1] + 1, iv[2]);
                    if ( isBoundary_m(yr) ) {
                        val += pot(yr) / ( dx[1] * dx[1] );
                    }
                    
                    IntVect zl(iv[0], iv[1], iv[2] - 1);
                    if ( isBoundary_m(zl) ) {
                        val += pot(zl) / ( dx[2] * dx[2] );
                    }
                    
                    IntVect zr(iv[0], iv[1], iv[2] + 1);
                    if ( isBoundary_m(zr) ) {
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
        if ( isInside_m(xl) ) {
            int gidx = serialize_m(xl/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
            ++numEntries;
        } 
        
        // check right neighbor in x-direction
        IntVect xr(iv[0] + 1, iv[1], iv[2]);
        if ( isInside_m(xr) ) {
            int gidx = serialize_m(xr/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
            ++numEntries;
        }

        // check lower neighbor in y-direction
        IntVect yl(iv[0], iv[1] - 1, iv[2]);
        if ( isInside_m(yl) ) {
            int gidx = serialize_m(yl/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
            ++numEntries;
        }

        // check upper neighbor in y-direction
        IntVect yr(iv[0], iv[1] + 1, iv[2]);
        if ( isInside_m(yr) ) {
            int gidx = serialize_m(yr/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
            ++numEntries;
        }
        
        // check front neighbor in z-direction
        IntVect zl(iv[0], iv[1], iv[2] - 1);
        if ( isInside_m(zl) ) {
            int gidx = serialize_m(zl/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
            ++numEntries;
        }

        // check back neighbor in z-direction
        IntVect zr(iv[0], iv[1], iv[2] + 1);
        if ( isInside_m(zr) ) {
            int gidx = serialize_m(zr/*shift*/);
            indices[numEntries] = gidx;
            values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
            ++numEntries;
        }
        
        // check center
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


inline bool AmrTrilinos::isBoundary_m(const IntVect& iv) const {
    return ( bc_m.find(iv) != bc_m.end() );
}


void AmrTrilinos::findBoundaryBoxes_m(const BoxArray& ba)
{
    // ghosts
    BoxList boundary_cells = amrex::GetBndryCells(ba, 1 /*ngrow, nghosts*/);
    
    BoxArray boundary(boundary_cells);

    if (boundary.size() == 0)
        throw std::runtime_error("No boundary!");    
    
    bc_m.clear();

    for (uint b = 0; b < boundary.size(); ++b) {
        const Box& bx = boundary[b];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    bc_m.insert(IntVect(i, j, k));
                }
            }
        }
    }
}


void AmrTrilinos::copyBack_m(AmrField_t& phi,
                             const Teuchos::RCP<Epetra_MultiVector>& sol,
                             const Geometry& geom)
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
    
    phi->FillBoundary(geom.periodicity());
    
}


void AmrTrilinos::interpFromCoarseLevel_m(container_t& phi,
                                          const Array<Geometry>& geom,
                                          int lev)
{
    NoOpPhysBC cphysbc, fphysbc;
//     PhysBCFunct cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
    CellConservativeProtected mapper;
    
    IntVect rr(2, 2, 2);
    
    phi[lev+1]->setVal(0.0, 1);
    
    phi[lev+1]->FillBoundary(geom[lev+1].periodicity());
    
    amrex::InterpFromCoarseLevel(*phi[lev+1], 0.0,
                                 *phi[lev],
                                 0, 0, 1,
                                 geom[lev], geom[lev+1],
                                 cphysbc, fphysbc,
                                 rr, &mapper, bcs);
    
    phi[lev+1]->FillBoundary(geom[lev+1].periodicity());
}


void AmrTrilinos::getGradient(const container_t& phi,
                              container_t& efield,
                              const Array<Geometry> geom)
{
    for (uint lev = 0; lev < phi.size(); ++lev) {
        const double* dx = geom[lev].CellSize();
        
        for (MFIter mfi(*phi[lev], false); mfi.isValid(); ++mfi) {
            const Box&          bx   = mfi.validbox();
            FArrayBox&          pfab = (*phi[lev])[mfi];
            FArrayBox&          efab = (*efield[lev])[mfi];
            
            const int* lo = bx.loVect();
            const int* hi = bx.hiVect();
            
            
            for (int i = lo[0]; i <= hi[0]; ++i) {
                for (int j = lo[1]; j <= hi[1]; ++j) {
                    for (int k = lo[2]; k <= hi[2]; ++k) {
                        
                        IntVect iv(i, j, k);
                        
                        // x direction
                        IntVect left(i-1, j, k);
                        IntVect right(i+1, j, k);
                        efab(iv, 0) = 0.5 * (pfab(left) - pfab(right)) / dx[0];
                        
                        // y direction
                        IntVect up(i, j+1, k);
                        IntVect down(i, j-1, k);
                        efab(iv, 1) = 0.5 * (pfab(down) - pfab(up)) / dx[1];
                        
                        // z direction
                        IntVect back(i, j, k+1);
                        IntVect front(i, j, k-1);
                        efab(iv, 2) = 0.5 * (pfab(front) - pfab(back)) / dx[2];
                    }
                }
            }
        }        
//         efield[lev]->FillBoundary(0, 3, geom[lev].periodicity());
    }
}