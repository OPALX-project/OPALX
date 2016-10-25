/* Author: Matthias Frey
 * 25. October 2016
 * 
 * 
 * 
 */

#include <iostream>

#include <BoxLib.H>
#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>

#include "Solver.h"
#include <fstream>

#include "writePlotFile.H"

double doSolve(int nr[3], int maxBoxSize, int nLevels) {
    // setup geometry
    IntVect lo(0, 0, 0);
    IntVect hi(nr[0] - 1, nr[1] - 1, nr[2] - 1);
    Box bx(lo, hi);
    
    RealBox domain;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain.setLo(i, 0.0);
        domain.setHi(i, 1.0);
    }
    
    // dirichlet boundary conditions
    int bc[BL_SPACEDIM] = {0, 0, 0};
    
    Array<Geometry> geom(nLevels);
    geom[0].define(bx, &domain, 0, bc);
    
    Array<BoxArray> ba(nLevels);
    
    ba[0].define(bx);
    ba[0].maxSize(maxBoxSize);
    
    Array<DistributionMapping> dmap(nLevels);
    dmap[0].define(ba[0], ParallelDescriptor::NProcs() /*nprocs*/);
    
    
    
    
    // ------------------------------------------------------------------------
    // Refined Meshes
    // ------------------------------------------------------------------------
    Array<int> refRatio(nLevels-1);
    for (int i = 0; i < refRatio.size(); ++i)
        refRatio[i] = 2;
    
    for (int lev = 1; lev < nLevels; ++lev) {
        geom[lev].define(BoxLib::refine(geom[lev - 1].Domain(),
                                        refRatio[lev - 1]),
                         &domain, 0, bc);
    }
    
    int fine = 1.0;
    for (int lev = 1; lev < nLevels; ++lev) {
        fine *= refRatio[lev - 1];
        
        IntVect refined_lo(0, 0, 0);
            
        IntVect refined_hi(nr[0] * fine - 1,
                           nr[1] * fine - 1,
                           nr[2] * fine - 1);
        
        Box refined_patch(refined_lo, refined_hi);
        ba[lev].define(refined_patch);
        
        ba[lev].maxSize(maxBoxSize); // / refRatio[lev - 1]);
        dmap[lev].define(ba[lev], ParallelDescriptor::NProcs() /*nprocs*/);
    }
    
    
    // ------------------------------------------------------------------------
    // Initialize MultiFabs
    // ------------------------------------------------------------------------
    
    
    PArray<MultiFab> rho(nLevels);
    PArray<MultiFab> phi(nLevels);
    PArray<MultiFab> efield(nLevels);
    PArray<MultiFab> phi_single(nLevels);
    
    for (int l = 0; l < nLevels; ++l) {
        rho.set(l, new MultiFab(ba[l], 1, 0));
        phi.set(l, new MultiFab(ba[l], 1, 1));
        phi_single.set(l, new MultiFab(ba[l], 1, 1));
        efield.set(l, new MultiFab(ba[l], BL_SPACEDIM, 1));
    }
    
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (multi level)
    // ------------------------------------------------------------------------
    int base_level = 0;
    int fine_level = nLevels - 1;
    
    // rho is equal to one everywhere
    for (int l = 0; l < nLevels; ++l)
        rho[l].setVal(-1.0);
    
    
    Real offset = 0.0;
    Solver sol;
    sol.solve_for_accel(rho, phi, efield, geom,
                        base_level, fine_level,
                        offset);
    
    // get solution on coarsest level by averaging down from finest.
    for (int i = fine_level-1; i >= 0; --i) {
        BoxLib::average_down(phi[i+1], phi[i], 0, 1, refRatio[i]);
    }
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (single level)
    // ------------------------------------------------------------------------
    base_level = 0;
    fine_level = 0;
    
    
    sol.solve_for_accel(rho, phi_single, efield, geom,
                        base_level, fine_level,
                        offset);
    
    
    // comparison
    MultiFab::Subtract(phi_single[0], phi[0], 0, 0, 1, 1);
    
//     std::string dir = "plt0000";
//     Real time = 0.0;
//     writePlotFile(dir, rho, phi, efield, refRatio, geom, time);
    
    return phi_single[0].norm2();
}

int main(int argc, char* argv[]) {
    
    BoxLib::Initialize(argc, argv, false);
    
    int nr[BL_SPACEDIM] = {
        std::atoi(argv[1]),
        std::atoi(argv[2]),
        std::atoi(argv[3])
    };
    
    
    int maxBoxSize = std::atoi(argv[4]);
    int nLevels = std::atoi(argv[5]);
    
    double l2error = doSolve(nr, maxBoxSize, nLevels);
    
    ParallelDescriptor::Barrier();
    
    if ( ParallelDescriptor::MyProc() == 0 ) {
        std::ofstream out("l2_error.dat", std::ios::app);
        out << nLevels << " " << l2error << std::endl;
        out.close();
    }
    
    BoxLib::Finalize();
    
    return 0;
}