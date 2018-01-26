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
    
    
    
    /*
     * find boundary cells --> mask
     */
    Array<BoxArray> grids(nLevel_m);
    Array<DistributionMapping> dmap(nLevel_m);
    Array<Geometry> gm(nLevel_m);
    masks_m.resize(nLevel_m);
    
    for (int lev = 0; lev < nLevel_m; ++lev) {
        grids[lev]   = rho[lev/* + lbase*/]->boxArray();
        dmap[lev]    = rho[lev/* + lbase*/]->DistributionMap();
        gm[lev]      = geom[lev/* + lbase*/];
        
        efield[lev]->setVal(0.0, 1);
        phi[lev]->setVal(0.0, 1);
        phi[lev]->setBndry(0.0);
        
        phi[lev]->FillBoundary(/*gm[lev].periodicity()*/);
        
        this->buildLevelMask_m(grids[lev], dmap[lev], gm[lev], lev);
    }
    
    
//     IntVect ratio(D_DECL(2, 2, 2));
//     fixRHSForSolve(rho, masks_m, geom, ratio);
    
    IntVect rr(2, 2, 2);
    
    bndry_m.reset(new MacBndry(grids[0], dmap[0], 1 /*ncomp*/, gm[0]));
    
    // initializes all ghosts to zero
    initGhosts_m();
    
    this->setBoundaryValue_m(phi[0]);
    
    
    for (OrientationIter oitr; oitr; ++oitr)
    {
        const Orientation face  = oitr();
        const FabSet& fs = bndry_m->bndryValues(oitr());
//         fs.plusTo(*phi[0], 1, 0, 0, 1, geom[0].periodicity());
        for (MFIter umfi(*(phi[0])); umfi.isValid(); ++umfi)
        {
            FArrayBox& dest = (*(phi[0]))[umfi];
            dest.copy(fs[umfi],fs[umfi].box());
        }
    }
    
    for (int i = 0; i < nLevel_m; ++i) {
//         fillDomainBoundary_m(*phi[i], geom[i]);
        phi[i]->FillBoundary();
    }
    
    for (int i = 0; i < nLevel_m; ++i) {
        int lev = i + lbase;
        
        
        if ( lev > 0 ) {
            // crse-fine interface
            interpFromCoarseLevel_m(phi, geom, lev-1);
//             phi[lev]->setBndry(0.0);
            phi[lev]->FillBoundary();
            
            
            
            bndry_m.reset(new MacBndry(grids[lev], dmap[lev], 1 /*ncomp*/, gm[lev]));
            initGhosts_m();
            this->setBoundaryValue_m(phi[lev], phi[lev-1], rr);
            
            for (OrientationIter oitr; oitr; ++oitr)
            {
                const Orientation face  = oitr();
                const FabSet& fs = bndry_m->bndryValues(oitr());
//                 fs.plusTo(*phi[lev], 1, 0, 0, 1, geom[lev].periodicity());
                for (MFIter umfi(*(phi[lev])); umfi.isValid(); ++umfi)
                {
                    FArrayBox& dest = (*(phi[lev]))[umfi];
                    dest.copy(fs[umfi],fs[umfi].box());
                }
            }
            
//             print_m(phi[lev]);
        }
        
        if ( phi[lev]->contains_nan() )
            throw std::runtime_error("Nan");
        
        solveLevel_m(phi[lev], rho[lev], geom[lev], lev);
        
        getGradient_1(phi[lev], efield[lev], geom[lev],  lev);
    }
    
    
//     for (int i = 0; i < nLevel_m; ++i) {
//         getGradient(phi[i], efield[i], geom[i],  i);
//     }
    
    // average down and compute electric field
//     for (int lev = lfine - 1; lev >= lbase; --lev) {
// //         amrex::average_down(*phi[lev+1], 
// //                              *phi[lev], geom[lev+1], geom[lev], 0, 1, rr);
//         amrex::average_down(*efield[lev+1],
//                             *efield[lev],
//                             geom[lev+1],
//                             geom[lev],
//                             0, 3, rr);
//     }
    
    
//     for (int lev = lfine - 1; lev >= lbase; --lev) {
//         amrex::average_down(*phi[lev+1], 
//                              *phi[lev], geom[lev+1], geom[lev], 0, 1, rr);
//         
//     }
}


void AmrTrilinos::solveLevel_m(AmrField_t& phi,
                               const AmrField_t& rho,
                               const Geometry& geom, int lev)
{
    for (int j = 0; j < 3; ++j)
        nPoints_m[j] = geom.Domain().length(j);
    
    build_m(rho, phi, geom, lev);
    fillMatrix_m(geom, phi, lev);
    linSysSolve_m(phi, geom);
}



void AmrTrilinos::linSysSolve_m(AmrField_t& phi, const Geometry& geom) {
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
    params->set("Convergence Tolerance", 1.0e-8);
    params->set("Block Size", 32);
    
    
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
#ifdef DEBUG
    if ( Ippl::myNode() == 0 ) {
        std::cout << "nx = " << nPoints_m[0] << std::endl
                  << "ny = " << nPoints_m[1] << std::endl
                  << "nz = " << nPoints_m[2] << std::endl;
    }
#endif
    
    indexMap_m.clear();

    int localNumElements = 0;
    std::vector<double> values, xvalues;
    std::vector<int> globalindices;
    
    const double* dx = geom.CellSize();
    
    
    
    for (MFIter mfi(*rho, false); mfi.isValid(); ++mfi) {
        const Box&          bx  = mfi.validbox();
        const FArrayBox&    fab = (*rho)[mfi];
        const FArrayBox&    pot = (*phi)[mfi];
        const BaseFab<int>& mfab = (*masks_m[lev])[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    
                    
                    IntVect iv(i, j, k);
                    
                    double val = 0.0;
                    
                    /* The boundary values (ghost cells) of level > 0 need
                     * to be the refined values of the next coarser
                     * level. --> add sum to right-hand side (Dirichlet)
                     */
                    IntVect xl(i-1, j, k);
                    if ( mfab(xl) > 0/*isBoundary_m(xl)*/ ) {
                        val += pot(xl) / ( dx[0] * dx[0] );
                    }
                        
                    IntVect xr(i+1, j, k);
                    if ( mfab(xr) > 0/*isBoundary_m(xr)*/ ) {
                        val += pot(xr) / ( dx[0] * dx[0] );
                    }
                    
                    IntVect yl(i, j-1, k);
                    if ( mfab(yl) > 0/*isBoundary_m(yl)*/ ) {
                        val += pot(yl) / ( dx[1] * dx[1] );
                    }
                    
                    IntVect yr(i, j+1, k);
                    if ( mfab(yr) > 0/*isBoundary_m(yr)*/ ) {
                        val += pot(yr) / ( dx[1] * dx[1] );
                    }
                    
                    IntVect zl(i, j, k-1);
                    if ( mfab(zl) > 0/*isBoundary_m(zl)*/ ) {
                        val += pot(zl) / ( dx[2] * dx[2] );
                    }
                    
                    IntVect zr(i, j, k+1);
                    if ( mfab(zr) > 0/*isBoundary_m(zr)*/ ) {
                        val += pot(zr) / ( dx[2] * dx[2] );
                    }
                    int globalidx = serialize_m(iv/*shift*/);
                    globalindices.push_back(globalidx);
                    
                    
                    // Invect (i, j, k) index --> 1D index
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
    
#ifdef DEBUG
    if ( Ippl::myNode() == 0 )
        std::cout << "N = " << N << std::endl;
#endif
    
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
    
    int indices[7] = {0, 0, 0, 0, 0, 0, 0};
    double values[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    const double* dx = geom.CellSize();
    
    for (MFIter mfi(*masks_m[lev], false); mfi.isValid(); ++mfi) {
        const Box&          bx  = mfi.validbox();
        const BaseFab<int>& mfab = (*masks_m[lev])[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    
                    double subx = 0.0;
                    double suby = 0.0;
                    double subz = 0.0;
                    
                    int numEntries = 0;
                    IntVect iv(i, j, k);
                    int globalRow = serialize_m(iv);
                    
                    // check left neighbor in x-direction
                    IntVect xl(i-1, j, k);
                    if ( mfab(xl) < 1 ) {
                        int gidx = serialize_m(xl/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
                        ++numEntries;
                    } else if ( mfab(xl) == 2 )
                        subx = 1.0;
                    
                    
                    // check right neighbor in x-direction
                    IntVect xr(i+1, j, k);
                    if ( mfab(xr) < 1 ) {
                        int gidx = serialize_m(xr/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
                        ++numEntries;
                    } else if ( mfab(xr) == 2 )
                        subx += 1.0;
                    
                    // check lower neighbor in y-direction
                    IntVect yl(i, j-1, k);
                    if ( mfab(yl) < 1 ) {
                        int gidx = serialize_m(yl/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
                        ++numEntries;
                    } else if ( mfab(yl) == 2 )
                        suby += 1.0;
                    
                    // check upper neighbor in y-direction
                    IntVect yr(i, j+1, k);
                    if ( mfab(yr) < 1 ) {
                        int gidx = serialize_m(yr/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
                        ++numEntries;
                    } else if ( mfab(yr) == 2 )
                        suby += 1.0;
                    
                    // check front neighbor in z-direction
                    IntVect zl(i, j, k-1);
                    if ( mfab(zl) < 1 ) {
                        int gidx = serialize_m(zl/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
                        ++numEntries;
                    } else if ( mfab(zl) == 2 )
                        subz += 1.0;
                    
                    // check back neighbor in z-direction
                    IntVect zr(i, j, k+1);
                    if ( mfab(zr) < 1 ) {
                        int gidx = serialize_m(zr/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
                        ++numEntries;
                    } else if ( mfab(zr) == 2 )
                        subz += 1.0;
                    
                    // check center
                    if ( mfab(iv) == 0 ) {
                        indices[numEntries] = globalRow;
                        values[numEntries]  = (2.0 + subx) / ( dx[0] * dx[0] ) +
                                        (2.0 + suby) / ( dx[1] * dx[1] ) +
                                        (2.0 + subz) / ( dx[2] * dx[2] );
                        ++numEntries;
                    }
                    
                    
                    int error = A_mp->InsertGlobalValues(globalRow, numEntries, &values[0], &indices[0]);
                    
                    if ( error != 0 )
                        throw std::runtime_error("Error in filling the matrix!");
                }
            }
        }
    }
    
    A_mp->FillComplete(true);
    
    /*
     * some printing
     */
    
//     EpetraExt::RowMatrixToMatlabFile("matrix.txt", *A_mp);
    
#ifdef DEBUG
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
#endif
}


// void AmrTrilinos::fillMatrix_m(const Geometry& geom, const AmrField_t& phi, int lev) {
//     // 3D --> 7 elements per row
//     A_mp = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map_mp, 7) );
//     
//     int indices[7] = {0, 0, 0, 0, 0, 0, 0};
//     double values[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//     
//     const double* dx = geom.CellSize();
//     
//     int * myGlobalElements = map_mp->MyGlobalElements();
//     for ( int i = 0; i < map_mp->NumMyElements(); ++i) {
//         /*
//          * GlobalRow	- (In) Row number (in global coordinates) to put elements.
//          * NumEntries	- (In) Number of entries.
//          * Values	- (In) Values to enter.
//          * Indices	- (In) Global column indices corresponding to values.
//          */
//         int globalRow = myGlobalElements[i];
//         int numEntries = 0;
//         
//         IntVect iv = indexMap_m[globalRow]/*deserialize_m(globalRow)*/;
//         
//         // check left neighbor in x-direction
//         IntVect xl(iv[0] - 1, iv[1], iv[2]);
//         if ( isInside_m(xl) ) {
//             int gidx = serialize_m(xl/*shift*/);
//             indices[numEntries] = gidx;
//             values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
//             ++numEntries;
//         } 
//         
//         // check right neighbor in x-direction
//         IntVect xr(iv[0] + 1, iv[1], iv[2]);
//         if ( isInside_m(xr) ) {
//             int gidx = serialize_m(xr/*shift*/);
//             indices[numEntries] = gidx;
//             values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
//             ++numEntries;
//         }
// 
//         // check lower neighbor in y-direction
//         IntVect yl(iv[0], iv[1] - 1, iv[2]);
//         if ( isInside_m(yl) ) {
//             int gidx = serialize_m(yl/*shift*/);
//             indices[numEntries] = gidx;
//             values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
//             ++numEntries;
//         }
// 
//         // check upper neighbor in y-direction
//         IntVect yr(iv[0], iv[1] + 1, iv[2]);
//         if ( isInside_m(yr) ) {
//             int gidx = serialize_m(yr/*shift*/);
//             indices[numEntries] = gidx;
//             values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
//             ++numEntries;
//         }
//         
//         // check front neighbor in z-direction
//         IntVect zl(iv[0], iv[1], iv[2] - 1);
//         if ( isInside_m(zl) ) {
//             int gidx = serialize_m(zl/*shift*/);
//             indices[numEntries] = gidx;
//             values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
//             ++numEntries;
//         }
// 
//         // check back neighbor in z-direction
//         IntVect zr(iv[0], iv[1], iv[2] + 1);
//         if ( isInside_m(zr) ) {
//             int gidx = serialize_m(zr/*shift*/);
//             indices[numEntries] = gidx;
//             values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
//             ++numEntries;
//         }
//         
//         // check center
//         if ( isInside_m(iv) ) {
//             indices[numEntries] = globalRow;
//             values[numEntries]  = 2.0 / ( dx[0] * dx[0] ) +
//                                2.0 / ( dx[1] * dx[1] ) +
//                                2.0 / ( dx[2] * dx[2] );
//             ++numEntries;
//         }
//         
//         
//         int error = A_mp->InsertGlobalValues(globalRow, numEntries, &values[0], &indices[0]);
//         
//         if ( error != 0 )
//             throw std::runtime_error("Error in filling the matrix!");
//     }
//     
//     A_mp->FillComplete(true);
//     
//     /*
//      * some printing
//      */
//     
// //     EpetraExt::RowMatrixToMatlabFile("matrix.txt", *A_mp);
//     
// #ifdef DEBUG
//     if ( Ippl::myNode() == 0 ) {
//         std::cout << "Global info" << std::endl
//                   << "Number of rows:      " << A_mp->NumGlobalRows() << std::endl
//                   << "Number of cols:      " << A_mp->NumGlobalCols() << std::endl
//                   << "Number of diagonals: " << A_mp->NumGlobalDiagonals() << std::endl
//                   << "Number of non-zeros: " << A_mp->NumGlobalNonzeros() << std::endl
//                   << std::endl;
//     }
//     
//     Ippl::Comm->barrier();
//     
//     for (int i = 0; i < Ippl::getNodes(); ++i) {
//         
//         if ( i == Ippl::myNode() ) {
//             std::cout << "Rank:                "
//                       << i << std::endl
//                       << "Number of rows:      "
//                       << A_mp->NumMyRows() << std::endl
//                       << "Number of cols:      "
//                       << A_mp->NumMyCols() << std::endl
//                       << "Number of diagonals: "
//                       << A_mp->NumMyDiagonals() << std::endl
//                       << "Number of non-zeros: "
//                       << A_mp->NumMyNonzeros() << std::endl
//                       << std::endl;
//         }
//         Ippl::Comm->barrier();
//     }
// #endif
// }



inline int AmrTrilinos::serialize_m(const IntVect& iv) const {
    return iv[0] + nPoints_m[0] * iv[1] + nPoints_m[0] * nPoints_m[1] * iv[2];
}


// inline IntVect AmrTrilinos::deserialize_m(int idx) const {
//     int i = idx % nPoints_m[0];
//     int j = ( idx / nPoints_m[0] ) % nPoints_m[1];
//     int k = ( idx / ( nPoints_m[0] * nPoints_m[1] ) ) % nPoints_m[2];
//     
//     return IntVect(i, j, k);
// }


inline bool AmrTrilinos::isInside_m(const IntVect& iv) const {
    
    return ( iv[0] > -1 && iv[0] < nPoints_m[0] &&
             iv[1] > -1 && iv[1] < nPoints_m[1] &&
             iv[2] > -1 && iv[2] < nPoints_m[2]);
}


// inline bool AmrTrilinos::isBoundary_m(const IntVect& iv) const {
//     return ( bc_m.find(iv) != bc_m.end() );
// }
// 
// 
// void AmrTrilinos::findBoundaryBoxes_m(const BoxArray& ba)
// {
//     // ghosts
//     BoxList boundary_cells = amrex::GetBndryCells(ba, 1 /*ngrow, nghosts*/);
//     
//     BoxArray boundary(boundary_cells);
// 
//     if (boundary.size() == 0)
//         throw std::runtime_error("No boundary!");    
//     
//     bc_m.clear();
// 
//     for (uint b = 0; b < boundary.size(); ++b) {
//         const Box& bx = boundary[b];
//         
//         const int* lo = bx.loVect();
//         const int* hi = bx.hiVect();
//         
//         for (int i = lo[0]; i <= hi[0]; ++i) {
//             for (int j = lo[1]; j <= hi[1]; ++j) {
//                 for (int k = lo[2]; k <= hi[2]; ++k) {
//                     bc_m.insert(IntVect(i, j, k));
//                 }
//             }
//         }
//     }
// }


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
    phi->FillBoundary();
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
    /*PCInterp*/CellConservativeProtected mapper;
    
    IntVect rr(2, 2, 2);
    
    amrex::InterpFromCoarseLevel(*phi[lev+1], 0.0,
                                 *phi[lev],
                                 0, 0, 1,
                                 geom[lev], geom[lev+1],
                                 cphysbc, fphysbc,
                                 rr, &mapper, bcs);
}


void AmrTrilinos::getGradient(AmrField_t& phi,
                              AmrField_t& efield,
                              const Geometry& geom, int lev)
{
    const double* dx = geom.CellSize();
    
    for (MFIter mfi(*masks_m[lev], false); mfi.isValid(); ++mfi) {
        const Box&          bx   = mfi.validbox();
        const FArrayBox&    pfab = (*phi)[mfi];
        FArrayBox&          efab = (*efield)[mfi];
        const BaseFab<int>& mfab = (*masks_m[lev])[mfi];
        
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
}



void AmrTrilinos::getGradient_1(AmrField_t& phi,
                               AmrField_t& efield,
                               const Geometry& geom, int lev)
{
    container_t grad_phi_edge(3);
    
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        BoxArray ba = phi->boxArray();
        const DistributionMapping& dmap = phi->DistributionMap();
        grad_phi_edge[n].reset(new MultiFab(ba.surroundingNodes(n), dmap, 1, 1));
        grad_phi_edge[n]->setVal(0.0, 1);
    }
    
    const double* dx = geom.CellSize();
    
    phi->FillBoundary();
    
    for (MFIter mfi(*masks_m[lev], false); mfi.isValid(); ++mfi) {
        const Box&          bx   = mfi.validbox();
        const FArrayBox&    pfab = (*phi)[mfi];
        FArrayBox&          xface = (*grad_phi_edge[0])[mfi];
        FArrayBox&          yface = (*grad_phi_edge[1])[mfi];
        FArrayBox&          zface = (*grad_phi_edge[2])[mfi];
        const BaseFab<int>& mfab = (*masks_m[lev])[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
            
        
        for (int i = lo[0]; i <= hi[0]+1; ++i) {
            for (int j = lo[1]; j <= hi[1]+1; ++j) {
                for (int k = lo[2]; k <= hi[2]+1; ++k) {
                    
                    IntVect iv(i, j, k);
                    
                    
                    
                    // x direction
                    IntVect left(i-1, j, k);
                    if ( mfab(left) == 2)
                        xface(iv) = -2.0 * pfab(iv) / dx[0];
                    else if ( mfab(iv) == 2 )
                        xface(iv) = 2.0 * pfab(left) / dx[0];
                    else
                        xface(iv) = (pfab(left) - pfab(iv)) / dx[0];
                    
                    // y direction
                    IntVect down(i, j-1, k);
                    
                    if ( mfab(down) == 2)
                        yface(iv) = -2.0 * pfab(iv) / dx[1];
                    else if ( mfab(iv) == 2 )
                        yface(iv) = 2.0 * pfab(down) / dx[1];
                    else
                        yface(iv) = (pfab(down) - pfab(iv)) / dx[1];
                    
                    // z direction
                    IntVect front(i, j, k-1);
                    
                    if ( mfab(front) == 2)
                        zface(iv) = -2.0 * pfab(iv) / dx[2];
                    else if ( mfab(iv) == 2 )
                        zface(iv) = 2.0 * pfab(front) / dx[2];
                    else
                        zface(iv) = (pfab(front) - pfab(iv)) / dx[2];
                }
            }
        }
    }
    
    
    amrex::average_face_to_cellcenter(*(efield.get()),
                                      amrex::GetArrOfConstPtrs(grad_phi_edge),
                                      geom);
        
    efield->FillBoundary(0, AMREX_SPACEDIM, geom.periodicity());
}


// void AmrTrilinos::solveWithGeometry(const container_t& rho,
//                                     container_t& phi,
//                                     container_t& efield,
//                                     const Array<Geometry>& geom,
//                                     int lbase,
//                                     int lfine)
// {
//     
//     double semimajor = 0.5;
//     double semiminor = 0.5;
//     Vector_t nr = Vector_t(geom[0].Domain().length(0),
//                            geom[0].Domain().length(1),
//                            geom[0].Domain().length(2));
//     
//     Vector_t hr = Vector_t(geom[0].CellSize(0),
//                            geom[0].CellSize(1),
//                            geom[0].CellSize(2));
//     
// //     std::cout << " nr = " << nr << " hr = " << hr << std::endl;
//     
//     std::string interpl = "CONSTANT";
//     EllipticDomain bp = EllipticDomain(semimajor, semiminor, nr, hr, interpl);
//     
//     bp.compute(hr);
//     
//     nLevel_m = lfine - lbase + 1;
//     
//     /*
//      * find boundary cells --> mask
//      */
//     Array<BoxArray> grids(nLevel_m);
//     Array<DistributionMapping> dmap(nLevel_m);
//     Array<Geometry> gm(nLevel_m);
//     masks_m.resize(nLevel_m);
//     
//     for (int lev = 0; lev < nLevel_m; ++lev) {
//         grids[lev]   = rho[lev/* + lbase*/]->boxArray();
//         dmap[lev]    = rho[lev/* + lbase*/]->DistributionMap();
//         gm[lev]      = geom[lev/* + lbase*/];
//         
//         efield[lev]->setVal(0.0, 1);
//         phi[lev]->setVal(0.0, 1);
//         
//         this->buildLevelMask_m(grids[lev], dmap[lev], gm[lev], lev);
//     }
//     
//     
//     for (int i = 0; i < nLevel_m; ++i) {
//         int lev = i + lbase;
//         
//         
//         if ( lev > 0 ) {
//             interpFromCoarseLevel_m(phi, geom, lev-1);
//         }
//         
//         if ( phi[lev]->contains_nan() )
//             throw std::runtime_error("Nan");
//         
//         solveLevelGeometry_m(phi[lev], rho[lev], geom[lev], lev, bp);
//     }
//     
//     
//     // average down and compute electric field
//     IntVect rr(2, 2, 2);
//     
//     
//     for (int i = 0; i < nLevel_m; ++i)
//         getGradient(phi[i], efield[i], geom[i],  i);
//     
//     
//     for (int lev = lfine - 1; lev >= lbase; --lev) {
//         amrex::average_down(*phi[lev+1], 
//                              *phi[lev], geom[lev+1], geom[lev], 0, 1, rr);
//         amrex::average_down(*efield[lev+1],
//                             *efield[lev],
//                             geom[lev+1],
//                             geom[lev],
//                             0, 3, rr);
//     }
// }
// 
// 
// void AmrTrilinos::buildGeometry_m(const AmrField_t& rho, const AmrField_t& phi,
//                           const Geometry& geom, int lev, EllipticDomain& bp)
// {
//     if ( Ippl::myNode() == 0 ) {
//         std::cout << "nx = " << nPoints_m[0] << std::endl
//                   << "ny = " << nPoints_m[1] << std::endl
//                   << "nz = " << nPoints_m[2] << std::endl;
//     }
//     
//     indexMap_m.clear();
// 
//     int localNumElements = 0;
//     std::vector<double> values, xvalues;
//     std::vector<int> globalindices;
//     
//     const double* dx = geom.CellSize();
//     
//     for (MFIter mfi(*rho, false); mfi.isValid(); ++mfi) {
//         const Box&          bx  = mfi.validbox();
//         const FArrayBox&    fab = (*rho)[mfi];
//         const FArrayBox&    pot = (*phi)[mfi];
// //         const BaseFab<int>& mfab = (*masks_m[lev])[mfi];
//             
//         const int* lo = bx.loVect();
//         const int* hi = bx.hiVect();
//         
//         for (int i = lo[0]; i <= hi[0]; ++i) {
//             for (int j = lo[1]; j <= hi[1]; ++j) {
//                 for (int k = lo[2]; k <= hi[2]; ++k) {
//                     
//                     
//                     
//                     IntVect iv(i, j, k);
//                     
//                     double val = 0.0;
//                     
// //                     /* The boundary values (ghost cells) of level > 0 need
// //                      * to be the refined values of the next coarser
// //                      * level. --> add sum to right-hand side (Dirichlet)
// //                      */
// //                     IntVect xl(iv[0] - 1, iv[1], iv[2]);
// //                     if ( mfab(xl) == 1/*isBoundary_m(xl)*/ ) {
// //                         val += pot(xl) / ( dx[0] * dx[0] );
// //                     }
// //                         
// //                     IntVect xr(iv[0] + 1, iv[1], iv[2]);
// //                     if ( mfab(xr) == 1/*isBoundary_m(xr)*/ ) {
// //                         val += pot(xr) / ( dx[0] * dx[0] );
// //                     }
// //                     
// //                     IntVect yl(iv[0], iv[1] - 1, iv[2]);
// //                     if ( mfab(yl) == 1/*isBoundary_m(yl)*/ ) {
// //                         val += pot(yl) / ( dx[1] * dx[1] );
// //                     }
// //                     
// //                     IntVect yr(iv[0], iv[1] + 1, iv[2]);
// //                     if ( mfab(yr) == 1/*isBoundary_m(yr)*/ ) {
// //                         val += pot(yr) / ( dx[1] * dx[1] );
// //                     }
// //                     
// //                     IntVect zl(iv[0], iv[1], iv[2] - 1);
// //                     if ( mfab(zl) == 1/*isBoundary_m(zl)*/ ) {
// //                         val += pot(zl) / ( dx[2] * dx[2] );
// //                     }
// //                     
// //                     IntVect zr(iv[0], iv[1], iv[2] + 1);
// //                     if ( mfab(zr) == 1/*isBoundary_m(zr)*/ ) {
// //                         val += pot(zr) / ( dx[2] * dx[2] );
// //                     }
//                     int globalidx = serialize_m(iv/*shift*/);
//                     globalindices.push_back(globalidx);
//                     
//                     
//                     // Invect (i, j, k) index --> 1D index
//                     indexMap_m[globalidx] = iv;
//                     
//                     if (bp.isInside(iv[0], iv[1], iv[2]))
//                         values.push_back(fab(iv) + val);
//                     else
//                         values.push_back(0.0);
//                     
//                     xvalues.push_back(pot(iv));
//                     
//                     ++localNumElements;
//                 }
//             }
//         }
//     }
//     
//     // compute map based on localelements
//     // create map that specifies which processor gets which data
//     const int baseIndex = 0;    // where to start indexing
//         
//     // numGlobalElements == N
//     const BoxArray& grids = rho->boxArray();
//     int N = grids.numPts();
//     
//     if ( Ippl::myNode() == 0 )
//         std::cout << "N = " << N << std::endl;
//     
//     std::cout << localNumElements << " " << globalindices.size() << std::endl;
//     
//     map_mp = Teuchos::rcp( new Epetra_Map(N, localNumElements,
//                                          &globalindices[0], baseIndex, epetra_comm_m) );
//     
//     // each processor fill in its part of the right-hand side
//     rho_mp = Teuchos::rcp( new Epetra_Vector(*map_mp, false) );
//     
//     // fill vector
//     int success = rho_mp->ReplaceGlobalValues(localNumElements,
//                                               &values[0],
//                                               &globalindices[0]);
//     
//     if ( success == 1 )
//         throw std::runtime_error("Error in filling the vector!");
//     
//     
//     // each processor fill in its part of the right-hand side
//     x_mp = Teuchos::rcp( new Epetra_Vector(*map_mp, false) ); // no initialization to zero
//     success = x_mp->ReplaceGlobalValues(localNumElements,
//                                               &xvalues[0],
//                                               &globalindices[0]);
//     
//     if ( success == 1 )
//         throw std::runtime_error("Error in filling the vector!");
//     
// //     std::cout << "Vector done." << std::endl;
// }
// 
// 
// void AmrTrilinos::solveLevelGeometry_m(AmrField_t& phi,
//                                        const AmrField_t& rho,
//                                        const Geometry& geom, int lev,
//                                        EllipticDomain& bp)
// {
//     for (int j = 0; j < 3; ++j)
//         nPoints_m[j] = geom.Domain().length(j);
//     
//     buildGeometry_m(rho, phi, geom, lev, bp);
//     fillMatrixGeometry_m(geom, phi, lev, bp);
//     std::cout << "Start solving." << std::endl;
//     linSysSolve_m(phi, geom);
// }
// 
// 
// void AmrTrilinos::fillMatrixGeometry_m(const Geometry& geom, const AmrField_t& phi, int lev,
//                                        EllipticDomain& bp)
// {
//     // 3D --> 7 elements per row
//     A_mp = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map_mp, 7) );
//     
//     int indices[7];
//     double values[7];
//     
//     const double* dx = geom.CellSize();
//     
//     int * myGlobalElements = map_mp->MyGlobalElements();
//     for ( int i = 0; i < map_mp->NumMyElements(); ++i) {
//         /*
//          * GlobalRow	- (In) Row number (in global coordinates) to put elements.
//          * NumEntries	- (In) Number of entries.
//          * Values	- (In) Values to enter.
//          * Indices	- (In) Global column indices corresponding to values.
//          */
//         int globalRow = myGlobalElements[i];
//         int numEntries = 0;
//         
//         IntVect iv = indexMap_m[globalRow]/*deserialize_m(globalRow)*/;
//         
//         
//         double WV, EV, SV, NV, FV, BV, CV, scaleFactor=1.0;
//         int W, E, S, N, F, B;
//         
//         bp.getBoundaryStencil(iv[0], iv[1], iv[2]/*myGlobalElements[i]*/, WV, EV, SV, NV, FV, BV, CV, scaleFactor);
//         
//         rho_mp->Values()[i] *= scaleFactor;
// 
//         bp.getNeighbours(iv[0], iv[1], iv[2]/*myGlobalElements[i]*/, W, E, S, N, F, B);
//         
//         if(E != -1) {
//             indices[numEntries] = E;
//             values[numEntries++] = EV;
//         }
//         if(W != -1) {
//             indices[numEntries] = W;
//             values[numEntries++] = WV;
//         }
//         if(S != -1) {
//             indices[numEntries] = S;
//             values[numEntries++] = SV;
//         }
//         if(N != -1) {
//             indices[numEntries] = N;
//             values[numEntries++] = NV;
//         }
//         if(F != -1) {
//             indices[numEntries] = F;
//             values[numEntries++] = FV;
//         }
//         if(B != -1) {
//             indices[numEntries] = B;
//             values[numEntries++] = BV;
//         }
//         
//         indices[numEntries] = globalRow;
//         values[numEntries++] = CV;
//         
//         int error = A_mp->InsertGlobalValues(globalRow, numEntries, &values[0], &indices[0]);
//         
//         if ( error != 0 )
//             throw std::runtime_error("Error in filling the matrix!");
//     }
//     
//     A_mp->FillComplete(true);
//     
//     /*
//      * some printing
//      */
//     
// //     EpetraExt::RowMatrixToMatlabFile("matrix.txt", *A_mp);
//     
//     if ( Ippl::myNode() == 0 ) {
//         std::cout << "Global info" << std::endl
//                   << "Number of rows:      " << A_mp->NumGlobalRows() << std::endl
//                   << "Number of cols:      " << A_mp->NumGlobalCols() << std::endl
//                   << "Number of diagonals: " << A_mp->NumGlobalDiagonals() << std::endl
//                   << "Number of non-zeros: " << A_mp->NumGlobalNonzeros() << std::endl
//                   << std::endl;
//     }
//     
//     Ippl::Comm->barrier();
//     
//     for (int i = 0; i < Ippl::getNodes(); ++i) {
//         
//         if ( i == Ippl::myNode() ) {
//             std::cout << "Rank:                "
//                       << i << std::endl
//                       << "Number of rows:      "
//                       << A_mp->NumMyRows() << std::endl
//                       << "Number of cols:      "
//                       << A_mp->NumMyCols() << std::endl
//                       << "Number of diagonals: "
//                       << A_mp->NumMyDiagonals() << std::endl
//                       << "Number of non-zeros: "
//                       << A_mp->NumMyNonzeros() << std::endl
//                       << std::endl;
//         }
//         Ippl::Comm->barrier();
//     }
// }
