#include "Ippl.h"
#include <AMReX.H>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#include <Amesos.h> // direct Ax = b solver

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>

#include <BelosBlockCGSolMgr.hpp>

#include <memory>

using namespace amrex;

struct TestParams {
    int nx;
    int ny;
    int nz;
    int max_grid_size;
    int nlevs;
};

Vector<std::unique_ptr<MultiFab> > initData(TestParams& parms,
                                           Vector<BoxArray>& ba,
                                           Vector<Geometry>& geom)
{
    
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, -1.0);
        real_box.setHi(n,  1.0);
    }
    
    
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(parms.nx - 1, parms.ny - 1, parms.nz-1); 
    const Box domain(domain_lo, domain_hi);
    
    
    Vector<int> rr(parms.nlevs-1);
    for (int lev = 1; lev < parms.nlevs; lev++)
        rr[lev-1] = 2;
    
    int coord = 0;
    
    int is_per[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) 
        is_per[i] = 0; 

    geom.resize(parms.nlevs);
    geom[0].define(domain, &real_box, coord, is_per);
    for (int lev = 1; lev < parms.nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, coord, is_per);
    }

    ba.resize(parms.nlevs);
    ba[0].define(domain);
    
    // Now we make the refined level be the center eighth of the domain
    if (parms.nlevs > 1) {
        int n_fine = parms.nx*rr[0];
        IntVect refined_lo(3*n_fine/8,3*n_fine/8,3*n_fine/8); 
        IntVect refined_hi(5*n_fine/8-1,5*n_fine/8-1,5*n_fine/8-1);

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < parms.nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
    }
    
    Vector<DistributionMapping> dmap(parms.nlevs);
    
    Vector<std::unique_ptr<MultiFab> > theData(parms.nlevs);
    
    for (int lev = 0; lev < parms.nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
        
        theData[lev].reset(new MultiFab(ba[lev], dmap[lev], 1, 0));

        theData[lev]->setVal(lev + 1.0);
    }
    
    return theData;
}

inline int serialize(int nx, int ny, int nz, int i, int j, int k) {
    return i + nx * j + nx * ny * k;
}

inline IntVect deserialize(int idx, int nx, int ny, int nz) {
    int i = idx % nx;
    int j = ( idx / nx ) % ny;
    int k = ( idx / ( nx * ny ) ) % nz;
    
    return IntVect(i, j, k);
}

// inline int localindex(int nx, int ny, int nz, int i, int j, int k, Epetra_Map& map) {
//     return serialize(nx, ny, nz, i, j, k) % map.NumMyElements();
// }

inline bool isBoundaryCell(IntVect iv, int nx, int ny, int nz) {
    return ((iv[0] % (nx - 1) == 0) ||
            (iv[1] % (ny - 1) == 0) ||
            (iv[2] % (nz - 1) == 0));
}

inline bool isInside(IntVect iv, int nx, int ny, int nz) {
    return ( iv[0] > -1 && iv[0] < nx &&
             iv[1] > -1 && iv[1] < ny &&
             iv[2] > -1 && iv[2] < nz);
}

void amrex3trilinos(const Vector<std::unique_ptr<MultiFab> >& theData,
                    const Vector<BoxArray>& grids,
                    const Vector<Geometry>& geom, int lev)
{
    
    Epetra_MpiComm epetra_comm(Ippl::getComm());
    
    
    Box minbox = grids[lev].minimalBox();
    
    // for boundary conditions
    std::map<int, IntVect> gidx2lidx;
    
    int nx = minbox.length(0);
    int ny = minbox.length(1);
    int nz = minbox.length(2);
    
    if ( Ippl::myNode() == 0 ) {
        std::cout << "nx = " << nx << std::endl
                  << "ny = " << ny << std::endl
                  << "nz = " << nz << std::endl;
    }
    
    const double* dx = geom[lev].CellSize();
    
    
    int localNumElements = 0;
    std::vector<double> values;
    std::vector<int> globalindices;
    
    for (MFIter mfi(*theData[lev], false); mfi.isValid(); ++mfi) {
        const Box&          bx  = mfi.validbox();
        const FArrayBox&    fab = (*theData[lev])[mfi];
        
//         std::cout << mfi.LocalIndex() << std::endl;
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    IntVect iv(i, j, k);
                    
                    int globalidx = serialize(nx, ny, nz, i, j, k);
                    globalindices.push_back(globalidx);
                    
                    gidx2lidx[globalidx] = iv;
                    
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
    int N = grids[lev].numPts();
    
    if ( Ippl::myNode() == 0 )
        std::cout << "N = " << N << std::endl;
        
    Epetra_Map map(N, localNumElements, &globalindices[0], baseIndex, epetra_comm);
    
    // each processor fill in its part of the right-hand side
    Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp( new Epetra_Vector(map, false) );
    
    
    std::cout << Ippl::myNode() << " " << map.NumMyElements() << std::endl;
    
    // fill vector
    int success = b->ReplaceGlobalValues(localNumElements,
                                         &values[0],
                                         &globalindices[0]);
    
    if ( success == 1 )
        std::cout << "Error in filling the vector!" << std::endl;
    
    if ( success == 0 && Ippl::myNode() == 0)
        std::cout << "Successfully filled vector!" << std::endl;
    
    
    /*
     * Fill matrix
     */
    // 3D --> 7 elements per row
    Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, map, 7) );
    
    int indices[7];
    double val[7];
    
    int * myGlobalElements = map.MyGlobalElements();
    for ( int i = 0; i < map.NumMyElements(); ++i) {
        /*
         * GlobalRow	- (In) Row number (in global coordinates) to put elements.
         * NumEntries	- (In) Number of entries.
         * Values	- (In) Values to enter.
         * Indices	- (In) Global column indices corresponding to values.
         */
        int globalRow = myGlobalElements[i];
        int numEntries = 0;
        
        IntVect iv = deserialize(globalRow, nx, ny, nz);
        
        // check left neighbor in x-direction
        IntVect xl(iv[0] - 1, iv[1], iv[2]);
        int gidx = serialize(nx, ny, nz, xl[0], xl[1], xl[2]);
        if ( isInside(xl, nx, ny, nz) ) {
            indices[numEntries] = gidx;
            val[numEntries]  = -1.0 / ( dx[0] * dx[0] );
            ++numEntries;
        }
        
        // check right neighbor in x-direction
        IntVect xr(iv[0] + 1, iv[1], iv[2]);
        gidx = serialize(nx, ny, nz, xr[0], xr[1], xr[2]);
        if ( isInside(xr, nx, ny, nz) ) {
            indices[numEntries] = gidx;
            val[numEntries]  = -1.0 / ( dx[0] * dx[0] );
            ++numEntries;
        }
        
        // check lower neighbor in y-direction
        IntVect yl(iv[0], iv[1] - 1, iv[2]);
        gidx = serialize(nx, ny, nz, yl[0], yl[1], yl[2]);
        if ( isInside(yl, nx, ny, nz) ) {
            indices[numEntries] = gidx;
            val[numEntries]  = -1.0 / ( dx[1] * dx[1] );
            ++numEntries;
        }
        
        // check upper neighbor in y-direction
        IntVect yr(iv[0], iv[1] + 1, iv[2]);
        gidx = serialize(nx, ny, nz, yr[0], yr[1], yr[2]);
        if ( isInside(yr, nx, ny, nz) ) {
            indices[numEntries] = gidx;
            val[numEntries]  = -1.0 / ( dx[1] * dx[1] );
            ++numEntries;
        }
        
        // check front neighbor in z-direction
        IntVect zl(iv[0], iv[1], iv[2] - 1);
        gidx = serialize(nx, ny, nz, zl[0], zl[1], zl[2]);
        if ( isInside(zl, nx, ny, nz) ) {
            indices[numEntries] = gidx;
            val[numEntries]  = -1.0 / ( dx[2] * dx[2] );
            ++numEntries;
        }
        
        // check back neighbor in z-direction
        IntVect zr(iv[0], iv[1], iv[2] + 1);
        gidx = serialize(nx, ny, nz, zr[0], zr[1], zr[2]);
        if ( isInside(zr, nx, ny, nz) ) {
            indices[numEntries] = gidx;
            val[numEntries]  = -1.0 / ( dx[2] * dx[2] );
            ++numEntries;
        }
        
        // check center
        gidx = globalRow;
        if ( isInside(gidx2lidx[globalRow], nx, ny, nz) ) {
            indices[numEntries] = globalRow;
            val[numEntries]  = 2.0 / ( dx[0] * dx[0] ) +
                               2.0 / ( dx[1] * dx[1] ) +
                               2.0 / ( dx[2] * dx[2] );
            
            ++numEntries;
        }
        
        int error = A->InsertGlobalValues(globalRow, numEntries, &val[0], &indices[0]);
        
        if ( error != 0 )
            std::cout << "Error in filling matrix" << std::endl;
    }
    
    A->FillComplete();
    
    /*
     * some printing
     */
    if ( Ippl::myNode() == 0 ) {
        std::cout << "Global info" << std::endl
                  << "Number of rows:      " << A->NumGlobalRows() << std::endl
                  << "Number of cols:      " << A->NumGlobalCols() << std::endl
                  << "Number of diagonals: " << A->NumGlobalDiagonals() << std::endl
                  << "Number of non-zeros: " << A->NumGlobalNonzeros() << std::endl
                  << std::endl;
    }
    
    Ippl::Comm->barrier();
    
    for (int i = 0; i < Ippl::getNodes(); ++i) {
        
        if ( i == Ippl::myNode() ) {
            std::cout << "Rank:                "
                      << i << std::endl
                      << "Number of rows:      "
                      << A->NumMyRows() << std::endl
                      << "Number of cols:      "
                      << A->NumMyCols() << std::endl
                      << "Number of diagonals: "
                      << A->NumMyDiagonals() << std::endl
                      << "Number of non-zeros: "
                      << A->NumMyNonzeros() << std::endl
                      << std::endl;
        }
        Ippl::Comm->barrier();
    }
    
//     std::cout << *b << std::endl;
//     std::cout << *A << std::endl;
    
    /*
     * unknowns
     */
    // each processor fill in its part of the right-hand side
    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp( new Epetra_Vector(map, false) ); // no initialization to zero
    x->PutScalar(0.0);
    
    
    
    /*
     * solve linear system Ax = b
     */
    Teuchos::RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator> > problem =
        Teuchos::rcp( new Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>(A, x, b) );
    bool set = problem->setProblem();
    if (set == false) {
        std::cerr << "ERROR: Belos::LinearProblem failed to set up correctly!" << std::endl;
        return;
    }
    
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp( new Teuchos::ParameterList );
    
//     params->set( "Block Size", 4 );
    params->set( "Maximum Iterations", 10000 );
    params->set("Convergence Tolerance", 1.0e-8);
    
    Belos::BlockCGSolMgr<double, Epetra_MultiVector, Epetra_Operator> solver(problem, params);
    
    Belos::ReturnType ret = solver.solve();
    
    // get the solution from the problem
    if ( ret == Belos::Converged ) {
        Teuchos::RCP<Epetra_MultiVector> sol = problem->getLHS();
        
        std::cout << *sol << std::endl;
        
        // print stuff
        if ( epetra_comm.MyPID() == 0 ) {
            std::cout << "Achieved tolerance: " << solver.achievedTol() << std::endl
                      << "Number of iterations: " << solver.getNumIters() << std::endl;
        }
        
    } else {
        if ( epetra_comm.MyPID() == 0 ) {
            std::cout << "Not converged. Achieved tolerance after " << solver.getNumIters() << " iterations is "
                      << solver.achievedTol() << "." << std::endl;
        }
    }
    
    
    /*
     * back to amrex
     */
    
    Teuchos::RCP<Epetra_MultiVector> sol = problem->getLHS();
    
    for (MFIter mfi(*theData[lev], false); mfi.isValid(); ++mfi) {
        const Box&          bx  = mfi.validbox();
        FArrayBox&          fab = (*theData[lev])[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    IntVect iv(i, j, k);
                    
                    int globalidx = serialize(nx, ny, nz, i, j, k);
                    
                    fab(iv) = (*sol)[0][globalidx];
                }
            }
        }
    }
}


int main(int argc, char* argv[]) {
    
    Ippl ippl(argc, argv);
    
    amrex::Initialize(argc, argv, true);
    
    
    ParmParse pp;
    TestParams parms;
  
    pp.get("nx", parms.nx);
    pp.get("ny", parms.ny);
    pp.get("nz", parms.nz);
    pp.get("max_grid_size", parms.max_grid_size);
    pp.get("nlevs", parms.nlevs);
    
    Vector<BoxArray> grids;
    Vector<Geometry> geom;
    Vector<std::unique_ptr<MultiFab> > theData = initData(parms, grids, geom);
    
    int lev = 0;
    amrex3trilinos(theData, grids, geom, lev);
    
    
    
    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);
    
    
    
    
    
    
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    return 0;
}