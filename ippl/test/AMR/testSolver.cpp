/*!
 * @file testSolver.cpp
 * @author Matthias Frey
 * @date 13. - 14. October 2016, LBNL
 * @details Compute \f$\Delta\phi = -1\f$ and write plot files
 * that can be visualized by yt (python visualize.py)
 * or by AmrVis of the CCSE group at LBNL.
 * 
 * Domain:  [0, 1] x [0, 1] x [0, 1]
 * BC:      Dirichlet boundary conditions
 * #Levels: single level (currently)
 * 
 * Call:
 *  mpirun -np [#cores] testSolver [#gridpints x] [#gridpints y] [#gridpints z]
 *      [max. grid size] [#levels] [plotfile]
 * 
 * The refinement is done over the whole domain.
 * @brief Computes \f$\Delta\phi = -1\f$ using only grids.
 */

#include <iostream>

#include <BoxLib.H>
#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>

#include "Solver.h"

#include "writePlotFile.H"

#include "Ippl.h"

int main(int argc, char* argv[]) {
    
    if (argc != 7) {
        std::cerr << "mpirun -np [#cores] testSolver [#gridpints x] "
                  << "[#gridpints y] [#gridpints z] [max. grid size] "
                  << "[#levels] [plotfile (boolean)]" << std::endl;
        return -1;
    }
    
    Ippl ippl(argc, argv);
    
    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);
    
    BoxLib::Initialize(argc, argv, false);
    
    
    static IpplTimings::TimerRef solverTimer = IpplTimings::getTimer("solver");
    
    int nr[BL_SPACEDIM] = {
        std::atoi(argv[1]),
        std::atoi(argv[2]),
        std::atoi(argv[3])
    };
    
    
    int maxBoxSize = std::atoi(argv[4]);
    int nLevels = std::atoi(argv[5]);
    bool doWrite = std::atoi(argv[6]);
    
    
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
    
    for (int l = 0; l < nLevels; ++l) {
        rho.set(l, new MultiFab(ba[l], 1, 0));
        phi.set(l, new MultiFab(ba[l], 1, 1));
        efield.set(l, new MultiFab(ba[l], BL_SPACEDIM, 1));
    }
    
    
    
    int base_level = 0;
    int fine_level = nLevels - 1;
    
    // rho is equal to one everywhere
    for (int l = 0; l < nLevels; ++l)
        rho[l].setVal(-1.0);
    
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid
    // ------------------------------------------------------------------------
    Real offset = 0.0;
    Solver sol;
    
    IpplTimings::startTimer(solverTimer);
    sol.solve_for_accel(rho, phi, efield, geom,
                        base_level, fine_level,
                        offset);
    IpplTimings::stopTimer(solverTimer);
    
    // ------------------------------------------------------------------------
    // Write BoxLib plotfile
    // ------------------------------------------------------------------------
    
    if (doWrite) {
        std::string dir = "plt0000";
        Real time = 0.0;
        writePlotFile(dir, rho, phi, efield, refRatio, geom, time);
    }
    
    IpplTimings::stopTimer(mainTimer);
    
    IpplTimings::print();
    
    std::string tfn = std::string(argv[0])
                        + "-cores=" + std::to_string(Ippl::getNodes())
                        + "-levels=" + std::to_string(nLevels)
                        + "-gridx=" + std::to_string(nr[0])
                        + "-gridy=" + std::to_string(nr[1])
                        + "-gridz=" + std::to_string(nr[2])
                        + "-gridsize=" + std::to_string(maxBoxSize);
                        
    IpplTimings::print(tfn);
    
    return 0;
}