#include <iostream>

#include <BoxLib.H>
#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>

#include "Solver.h"

int main(int argc, char* argv[]) {
    
    BoxLib::Initialize(argc, argv, false);
    
    int nLevels = 1;
    int nr[BL_SPACEDIM] = {64, 64, 64};
    int maxBoxSize = 16;
    
    
    // setup geometry
    IntVect lo(0, 0, 0);
    IntVect hi(nr[0] - 1, nr[1] - 1, nr[2] - 1);
    Box bx(lo, hi);
    
    RealBox domain;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain.setLo(i, 0.0);
        domain.setHi(i, 1.0);
    }
    
    // periodic boundary conditions
    int bc[BL_SPACEDIM] = {1, 1, 1};
    
    Array<Geometry> geom(nLevels);
    geom[0].define(bx, &domain, 0, bc);
    
    Array<BoxArray> ba(nLevels);
    
    ba[0].define(bx);
    ba[0].maxSize(maxBoxSize);
    
    Array<DistributionMapping> dmap(nLevels);
    dmap[0].define(ba[0], ParallelDescriptor::NProcs() /*nprocs*/);
    
    
    PArray<MultiFab> rho(nLevels);
    PArray<MultiFab> phi(nLevels);
    PArray<MultiFab> efield(nLevels);
    
    for (int l = 0; l < nLevels; ++l) {
        rho.set(l, new MultiFab(ba[l], 1, 0));
        phi.set(l, new MultiFab(ba[l], 1, 1));
        efield.set(l, new MultiFab(ba[l], BL_SPACEDIM, 1));
    }
    
    int base_level = 0;
    int fine_level = nLevels - 1;
    
    // rho is equal to one everywhere
    for (int l = 0; l < nLevels; ++l)
        rho[l].setVal(1.0);
    
    
    Real offset = 0.0;
    Solver sol;
    sol.solve_for_accel(rho, phi, efield, geom,
                        base_level, fine_level,
                        offset);
    
    
    for (int l = 0; l < nLevels; ++l) {
        std::cout << "[" << rho[l].min(0) << " " << rho[l].max(0) << "]" << std::endl
                  << "[" << phi[l].min(0) << " " << rho[l].max(0) << "]" << std::endl
                  << "[ (" << efield[l].min(0) << ", " << efield[l].min(1) << ", "
                  << efield[l].min(2) << ") , (" << efield[l].max(0) << ", "
                  << efield[l].max(1) << ", " << efield[l].max(2) << ") ]" << std::endl;
                  
    }
    
    BoxLib::Finalize();
    
    return 0;
}