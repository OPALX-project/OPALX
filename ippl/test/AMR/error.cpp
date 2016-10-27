/* Author: Matthias Frey
 * 25. - 26. October 2016
 * 
 * Solve \Lap\ph = \rho
 * for different ways of \rho initializations.
 * 
 * doSolve: \rho = -1 everywhere (no particles)
 * doSolve1: \rho = -1 everywhere (initialized by particles, idea: we generate particles on the
 *           cell center and then do an AssignDensity-function call --> \rho = -1 everywhere
 * doSolve2: \rho is a Gaussian distribution initialized by particles
 */

#include <fstream>
#include <iostream>

#include <BoxLib.H>
#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>

#include "Solver.h"
#include "AmrPartBunch.h"
#include "Distribution.h"

#include "writePlotFile.H"

// solve \Lap\phi = - 1
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
        
        phi[l].setVal(0.0);
        phi_single[l].setVal(0.0);
        efield[l].setVal(0.0);
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
    
    
    std::string dir = "plt0000";
    Real time = 0.0;
    writePlotFile(dir, rho, phi, efield, refRatio, geom, time);
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (single level)
    // ------------------------------------------------------------------------
    base_level = 0;
    fine_level = 0;
    
    
    sol.solve_for_accel(rho, phi_single, efield, geom,
                        base_level, fine_level,
                        offset);
    
    
    // comparison
    std::cout << phi[0].sum() << " " << phi_single[0].sum() << std::endl;
    MultiFab::Subtract(phi_single[0], phi[0], 0, 0, 1, 1);
    
    return phi_single[0].norm2();
}

/* solve \Lap\phi = -1
 * where we use particles and assign density
 */
double doSolve1(int nr[3], int maxBoxSize, int nLevels) {
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
    
    Array<int> refRatio(nLevels - 1);
    for (int i = 0; i < refRatio.size(); ++i)
        refRatio[i] = 2;
    
    PartBunchBase* bunch = new AmrPartBunch(geom[0], dmap[0], ba[0]);
    
    int nParticles = std::pow(2 /* refinement ratio */, 3*(nLevels-1)) * nr[0] * nr[1] * nr[2];
    
    bunch->create(nParticles / ParallelDescriptor::NProcs());
    
    // each particle is in the center of a cell
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
        bunch->setQM(-1.0 / nParticles , i);
    }
    
    // ------------------------------------------------------------------------
    // Refined Meshes
    // ------------------------------------------------------------------------
    
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
    
    int cnt = 0;
    double ref = std::pow(2, nLevels - 1);
    
    Real dx = *(geom[nLevels - 1].CellSize()) * 0.5;
    
    double inx = 1.0 / (nr[0] * ref);
    double iny = 1.0 / (nr[1] * ref);
    double inz = 1.0 / (nr[2] * ref);
        
    for (int j = 0; j < ba[nLevels - 1].size(); ++j) {
        Box bx = ba[nLevels - 1].get(j);
        
        for (int k = bx.loVect()[0]; k <= bx.hiVect()[0]; ++k)
            for (int l = bx.loVect()[1]; l <= bx.hiVect()[1]; ++l)
                for (int m = bx.loVect()[2]; m <= bx.hiVect()[2]; ++m) {
                    bunch->setR(Vector_t(k * inx + dx,
                                         l * iny + dx,
                                         m * inz + dx),
                                cnt++);
                }
        
    }
    
//     for (unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
//         std::cout << bunch->getR(i)(0) << " " << bunch->getR(i)(1) << " "
//                   << bunch->getR(i)(2) << std::endl;
//     }
    
    
    dynamic_cast<AmrPartBunch*>(bunch)->Define (geom, dmap, ba, refRatio);
    
    bunch->myUpdate();
    
    bunch->gatherStatistics();
    
    // ------------------------------------------------------------------------
    // Initialize MultiFabs
    // ------------------------------------------------------------------------
    
    
    PArray<MultiFab> rho(nLevels);
    PArray<MultiFab> phi(nLevels);
    PArray<MultiFab> efield(nLevels);
    PArray<MultiFab> phi_single(nLevels);
    
    for (int l = 0; l < nLevels; ++l) {
        rho.set(l, new MultiFab(ba[l], 1, 0));
        rho[l].setVal(0.0);
        
        phi.set(l, new MultiFab(ba[l], 1, 1));
        phi[l].setVal(0.0);
        
        phi_single.set(l, new MultiFab(ba[l], 1, 1));
        phi_single[l].setVal(0.0);
        
        efield.set(l, new MultiFab(ba[l], BL_SPACEDIM, 1));
        efield[l].setVal(0.0);
    }
    
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (multi level)
    // ------------------------------------------------------------------------
    int base_level = 0;
    int fine_level = nLevels - 1;
    
    
    PArray<MultiFab> PartMF;
    PartMF.resize(nLevels,PArrayManage);
    
    for (int i = 0; i < nLevels; ++i) {
        PartMF.set(i, new MultiFab(ba[i], 1, 0));
        PartMF[i].setVal(0.0);
    }
    
    dynamic_cast<AmrPartBunch*>(bunch)->SetAllowParticlesNearBoundary(true);
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, PartMF, base_level, 1, fine_level);

    for (int lev = fine_level - 1 - base_level; lev >= 0; lev--)
        BoxLib::average_down(PartMF[lev+1],PartMF[lev], 0, 1, refRatio[lev]);

    for (int lev = 0; lev < nLevels; lev++)
        MultiFab::Add(rho[base_level+lev], PartMF[lev], 0, 0, 1, 0);
    
    Real offset = 0.0;
    Solver sol;
    sol.solve_for_accel(rho, phi, efield, geom,
                        base_level, fine_level,
                        offset);
    
    // get solution on coarsest level by averaging down from finest.
    for (int i = fine_level-1; i >= 0; --i) {
        BoxLib::average_down(phi[i+1], phi[i], 0, 1, refRatio[i]);
    }
    
    std::string dir = "plt0000";
    Real time = 0.0;
    writePlotFile(dir, rho, phi, efield, refRatio, geom, time);
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (single level)
    // ------------------------------------------------------------------------
    base_level = 0;
    fine_level = 0;
    
    
    sol.solve_for_accel(rho, phi_single, efield, geom,
                        base_level, fine_level,
                        offset);
    
    
    
    // comparison
    std::cout << phi[0].sum() << " " << phi_single[0].sum() << std::endl;
    MultiFab::Subtract(phi_single[0], phi[0], 0, 0, 1, 1);
    
    return phi_single[0].norm2();
}


double doSolve2(int nr[3], int maxBoxSize, int nLevels, int nParticles) {
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
    
    Array<int> refRatio(nLevels - 1);
    for (int i = 0; i < refRatio.size(); ++i)
        refRatio[i] = 2;
    
    PartBunchBase* bunch = new AmrPartBunch(geom[0], dmap[0], ba[0]);
    
    // ------------------------------------------------------------------------
    // Refined Meshes
    // ------------------------------------------------------------------------
    
    for (int lev = 1; lev < nLevels; ++lev) {
        geom[lev].define(BoxLib::refine(geom[lev - 1].Domain(),
                                        refRatio[lev - 1]),
                         &domain, 0, bc);
    }
    
    int fine = 1.0;
    for (int lev = 1; lev < nLevels; ++lev) {
        fine *= refRatio[lev - 1];
        
        if ( lev == 1) {
            IntVect refined_lo(0.25 * nr[0] * fine,
                               0.25 * nr[1] * fine,
                               0.25 * nr[2] * fine);
            
            IntVect refined_hi(0.75 * nr[0] * fine - 1,
                               0.75 * nr[1] * fine - 1,
                               0.75 * nr[2] * fine - 1);
        
            Box refined_patch(refined_lo, refined_hi);
            ba[lev].define(refined_patch);
            
        } else if ( lev == 2 ) {
            IntVect refined_lo(0.375 * nr[0] * fine,
                               0.375 * nr[1] * fine,
                               0.375 * nr[2] * fine);
            
            IntVect refined_hi(0.625 * nr[0] * fine - 1,
                               0.625 * nr[1] * fine - 1,
                               0.625 * nr[2] * fine - 1);
        
            Box refined_patch(refined_lo, refined_hi);
            ba[lev].define(refined_patch);
        }
        
        ba[lev].maxSize(maxBoxSize);
        dmap[lev].define(ba[lev], ParallelDescriptor::NProcs() /*nprocs*/);
    }
    
    int nloc = nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    dist.gaussian(0.5, 0.1, nloc, ParallelDescriptor::MyProc());
    dist.injectBeam(*bunch);
    
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
        bunch->setQM(-1.0, i);
    }
    
//     for (unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
//         std::cout << bunch->getR(i)(0) << " " << bunch->getR(i)(1) << " "
//                   << bunch->getR(i)(2) << std::endl;
//     }
    
    
    dynamic_cast<AmrPartBunch*>(bunch)->Define (geom, dmap, ba, refRatio);
    
    bunch->myUpdate();
    
    bunch->gatherStatistics();
    
    // ------------------------------------------------------------------------
    // Initialize MultiFabs
    // ------------------------------------------------------------------------
    
    
    PArray<MultiFab> rho(nLevels);
    PArray<MultiFab> phi(nLevels);
    PArray<MultiFab> efield(nLevels);
    PArray<MultiFab> phi_single(nLevels);
    
    for (int l = 0; l < nLevels; ++l) {
        rho.set(l, new MultiFab(ba[l], 1, 0));
        rho[l].setVal(0.0);
        
        phi.set(l, new MultiFab(ba[l], 1, 1));
        phi[l].setVal(0.0);
        
        phi_single.set(l, new MultiFab(ba[l], 1, 1));
        phi_single[l].setVal(0.0);
        
        efield.set(l, new MultiFab(ba[l], BL_SPACEDIM, 1));
        efield[l].setVal(0.0);
    }
    
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (multi level)
    // ------------------------------------------------------------------------
    int base_level = 0;
    int fine_level = nLevels - 1;
    
    
    PArray<MultiFab> PartMF;
    PartMF.resize(nLevels,PArrayManage);
    
    for (int i = 0; i < nLevels; ++i) {
        PartMF.set(i, new MultiFab(ba[i], 1, 0));
        PartMF[i].setVal(0.0);
    }
    
    dynamic_cast<AmrPartBunch*>(bunch)->SetAllowParticlesNearBoundary(true);
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, PartMF, base_level, 1, fine_level);

    for (int lev = fine_level - 1 - base_level; lev >= 0; lev--)
        BoxLib::average_down(PartMF[lev+1],PartMF[lev], 0, 1, refRatio[lev]);

    for (int lev = 0; lev < nLevels; lev++)
        MultiFab::Add(rho[base_level+lev], PartMF[lev], 0, 0, 1, 0);
    
    Real offset = 0.0;
    Solver sol;
    sol.solve_for_accel(rho, phi, efield, geom,
                        base_level, fine_level,
                        offset);
    
    // get solution on coarsest level by averaging down from finest.
    for (int i = fine_level-1; i >= 0; --i) {
        BoxLib::average_down(phi[i+1], phi[i], 0, 1, refRatio[i]);
    }
    
    std::string dir = "plt0000";
    Real time = 0.0;
    writePlotFile(dir, rho, phi, efield, refRatio, geom, time);
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (single level)
    // ------------------------------------------------------------------------
    base_level = 0;
    fine_level = 0;
    
    
    sol.solve_for_accel(rho, phi_single, efield, geom,
                        base_level, fine_level,
                        offset);
    
    
    
    // comparison
    std::cout << phi[0].sum() << " " << phi_single[0].sum() << std::endl;
    MultiFab::Subtract(phi_single[0], phi[0], 0, 0, 1, 1);
    
    return phi_single[0].norm2();
}

int main(int argc, char* argv[]) {
    
    if ( argc != 8 ) {
        std::cerr << "mpirun -np [#cores] testError [gridx] [gridy] [gridz] "
                  << "[max. grid] [#levels] "
                  << "[NOPARTICLES/UNIFORM/GAUSSIAN] "
                  << "[#particles]" << std::endl;
        return -1;
    }
    
    BoxLib::Initialize(argc, argv, false);
    
    int nr[BL_SPACEDIM] = {
        std::atoi(argv[1]),
        std::atoi(argv[2]),
        std::atoi(argv[3])
    };
    
    
    int maxBoxSize = std::atoi(argv[4]);
    int nLevels = std::atoi(argv[5]);
    int nParticles = std::atoi(argv[7]);
    
    double l2error = 0.0;
    bool solved = true;
    
    if ( std::strcmp(argv[6], "NOPARTICLES") == 0 )
        l2error = doSolve(nr, maxBoxSize, nLevels);
    else if ( std::strcmp(argv[6], "UNIFORM") == 0 )
        l2error = doSolve1(nr, maxBoxSize, nLevels);
    else if ( std::strcmp(argv[6], "GAUSSIAN") == 0 )
        l2error = doSolve2(nr, maxBoxSize, nLevels, nParticles);
    else {
        if ( ParallelDescriptor::MyProc() == 0 )
            std::cerr << "Not supported solver." << std::endl
                      << "Use:" << std::endl
                      << "- NOPARTICLES: rhs = -1 everywhere" << std::endl
                      << "- UNIFORM: As NOPARTICLES but with particles (single-core)"
                      << std::endl
                      << "- GAUSSIAN: rhs using particles" << std::endl;
        solved = false;
    }
    
    ParallelDescriptor::Barrier();
    
    if ( ParallelDescriptor::MyProc() == 0 && solved) {
        std::ofstream out("l2_error.dat", std::ios::app);
        out << nLevels << " " << l2error << std::endl;
        out.close();
    }
    
    BoxLib::Finalize();
    
    return 0;
}