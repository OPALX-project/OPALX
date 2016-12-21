#include "PoissonProblems.h"

#include "AmrOpal.h"

#include "Physics/Physics.h"

PoissonProblems::PoissonProblems(int nr[3], int maxGridSize, int nLevels,
                                 const std::vector<double>& lower,
                                 const std::vector<double>& upper)
    : maxGridSize_m(maxGridSize), nLevels_m(nLevels)
{
    nr_m[0] = nr[0];
    nr_m[1] = nr[1];
    nr_m[2] = nr[2];
    
    // setup geometry
    IntVect lo(0, 0, 0);
    IntVect hi(nr[0] - 1, nr[1] - 1, nr[2] - 1);
    Box bx(lo, hi);
    
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain_m.setLo(i, lower[i]);
        domain_m.setHi(i, upper[i]);
    }
    
    // dirichlet boundary conditions
    int bc[BL_SPACEDIM] = {0, 0, 0};
    
    geom_m.resize(nLevels);
    geom_m[0].define(bx, &domain_m, 0, bc);
    
    ba_m.resize(nLevels);
    ba_m[0].define(bx);
    ba_m[0].maxSize(maxGridSize);
    
    dmap_m.resize(nLevels);
    dmap_m[0].define(ba_m[0], ParallelDescriptor::NProcs() /*nprocs*/);
    
    // refinement ratios
    refRatio_m.resize(nLevels - 1);
    for (unsigned int i = 0; i < refRatio_m.size(); ++i)
        refRatio_m[i] = 2;
    
    // geometries of refined levels
    for (int lev = 1; lev < nLevels; ++lev) {
        geom_m[lev].define(BoxLib::refine(geom_m[lev - 1].Domain(),
                                          refRatio_m[lev - 1]),
                           &domain_m, 0, bc);
    }
}


double PoissonProblems::doSolveNoParticles() {
    
    // ------------------------------------------------------------------------
    // Refined Meshes
    // ------------------------------------------------------------------------
    refineWholeDomain_m();
    
    // ------------------------------------------------------------------------
    // Initialize MultiFabs (rho, phi, efield, ...)
    // ------------------------------------------------------------------------
    initMultiFabs_m();
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (multi level)
    // ------------------------------------------------------------------------
    int base_level = 0;
    int fine_level = nLevels_m - 1;
    
    // rho is equal to one everywhere
    for (int l = 0; l < nLevels_m; ++l)
        rho_m[l]->setVal(-1.0);
    
    
    Real offset = 0.0;
    Solver sol;
    sol.solve_for_accel(rho_m, phi_m, efield_m, geom_m,
                        base_level, fine_level,
                        offset);
    
    // get solution on coarsest level by averaging down from finest.
    for (int i = fine_level-1; i >= 0; --i) {
        BoxLib::average_down(*phi_m[i+1], *phi_m[i], 0, 1, refRatio_m[i]);
    }
    
    
    std::string dir = "plt0000";
    Real time = 0.0;
    writePlotFile(dir, rho_m, phi_m, efield_m, refRatio_m, geom_m, time);
    
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (single level)
    // ------------------------------------------------------------------------
    base_level = 0;
    fine_level = 0;
    
    
    sol.solve_for_accel(rho_m, phi_single_m, efield_m, geom_m,
                        base_level, fine_level,
                        offset);
    
    
    // comparison
    MultiFab::Subtract(*phi_single_m[0], *phi_m[0], 0, 0, 1, 1);
    
    return phi_single_m[0]->norm2();
}


double PoissonProblems::doSolveParticlesUniform() {
    
    PartBunchBase* bunch = new AmrPartBunch(geom_m[0], dmap_m[0], ba_m[0]);
    
    int nParticles = std::pow(2 /* refinement ratio */, 3*(nLevels_m - 1)) * nr_m[0] * nr_m[1] * nr_m[2];
    
    bunch->create(nParticles / ParallelDescriptor::NProcs());
    
    // each particle is in the center of a cell
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
        bunch->setQM(-1.0 / nParticles , i);
    
    
    // ------------------------------------------------------------------------
    // Refined Meshes
    // ------------------------------------------------------------------------
    refineWholeDomain_m();
    
    
    // ------------------------------------------------------------------------
    // Generate particles in the center of the finest cell
    // ------------------------------------------------------------------------
    int cnt = 0;
    double ref = std::pow(2, nLevels_m - 1);
    
    Real dx = *(geom_m[nLevels_m - 1].CellSize()) * 0.5;
    
    double inx = 1.0 / (nr_m[0] * ref);
    double iny = 1.0 / (nr_m[1] * ref);
    double inz = 1.0 / (nr_m[2] * ref);
        
    for (int j = 0; j < ba_m[nLevels_m - 1].size(); ++j) {
        Box bx = ba_m[nLevels_m - 1].get(j);
        
        for (int k = bx.loVect()[0]; k <= bx.hiVect()[0]; ++k)
            for (int l = bx.loVect()[1]; l <= bx.hiVect()[1]; ++l)
                for (int m = bx.loVect()[2]; m <= bx.hiVect()[2]; ++m) {
                    bunch->setR(Vector_t(k * inx + dx,
                                         l * iny + dx,
                                         m * inz + dx),
                                cnt++);
                }
    }
    
    
    // ------------------------------------------------------------------------
    // Redistribute particles on multi-level
    // ------------------------------------------------------------------------
    dynamic_cast<AmrPartBunch*>(bunch)->Define(geom_m, dmap_m, ba_m, refRatio_m);
    
    bunch->myUpdate();
    
    bunch->gatherStatistics();
    
    // ------------------------------------------------------------------------
    // Initialize MultiFabs (rho, phi, efield, ...)
    // ------------------------------------------------------------------------
    initMultiFabs_m();
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (multi level)
    // ------------------------------------------------------------------------
    int base_level = 0;
    int fine_level = nLevels_m - 1;
    
    
//     PArray<MultiFab> PartMF;
//     PartMF.resize(nLevels_m,PArrayManage);
    
//     for (int i = 0; i < nLevels_m; ++i) {
//         PartMF.set(i, new MultiFab(ba_m[i], 1, 0));
//         PartMF[i].setVal(0.0);
//     }
    
    dynamic_cast<AmrPartBunch*>(bunch)->SetAllowParticlesNearBoundary(true);
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rho_m, base_level, 1, fine_level);

//     for (int lev = fine_level - 1 - base_level; lev >= 0; lev--)
//         BoxLib::average_down(PartMF[lev+1],PartMF[lev], 0, 1, refRatio_m[lev]);

//     for (int lev = 0; lev < nLevels_m; lev++)
//         MultiFab::Add(*rho_m[base_level+lev], PartMF[lev], 0, 0, 1, 0);
    
    for (int i = 0; i < nLevels_m; ++i)
        std::cout << rho_m[i]->min(0) << " " << rho_m[i]->max(0) << std::endl;
    
    Real offset = 0.0;
    Solver sol;
    sol.solve_for_accel(rho_m, phi_m, efield_m, geom_m,
                        base_level, fine_level,
                        offset);
    
    // get solution on coarsest level by averaging down from finest.
    for (int i = fine_level-1; i >= 0; --i) {
        BoxLib::average_down(*phi_m[i+1], *phi_m[i], 0, 1, refRatio_m[i]);
    }
    
    std::string dir = "plt0000";
    Real time = 0.0;
    writePlotFile(dir, rho_m, phi_m, efield_m, refRatio_m, geom_m, time);
    
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (single level)
    // ------------------------------------------------------------------------
    base_level = 0;
    fine_level = 0;
    
    
    sol.solve_for_accel(rho_m, phi_single_m, efield_m, geom_m,
                        base_level, fine_level,
                        offset);
    
    
    
    // comparison
    MultiFab::Subtract(*phi_single_m[0], *phi_m[0], 0, 0, 1, 1);
    
    return phi_single_m[0]->norm2();
}



double PoissonProblems::doSolveParticlesGaussian(int nParticles, double mean, double stddev) {
    // ------------------------------------------------------------------------
    // Generate particles
    // ------------------------------------------------------------------------
    PartBunchBase* bunch = new AmrPartBunch(geom_m[0], dmap_m[0], ba_m[0]);
    
    int nloc = nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    dist.gaussian(mean, stddev, nloc, ParallelDescriptor::MyProc());
    dist.injectBeam(*bunch);
    
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
        bunch->setQM(1.0, i);
    
    
    // redistribute on multi-level
    dynamic_cast<AmrPartBunch*>(bunch)->Define (geom_m, dmap_m, ba_m, refRatio_m);
    
    bunch->myUpdate();
    
    bunch->gatherStatistics();
    
    
    // ------------------------------------------------------------------------
    // Do AMR with tagging
    // ------------------------------------------------------------------------
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr_m[i];
    
    ParmParse pp("amr");
    pp.add("max_grid_size", maxGridSize_m);
    
    AmrOpal myAmrOpal(&domain_m, nLevels_m - 1, nCells, 0 /* cartesian */, bunch);
    
    
    for (int i = 0; i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    dynamic_cast<AmrPartBunch*>(bunch)->Define (myAmrOpal.Geom(),
                                                myAmrOpal.DistributionMap(),
                                                myAmrOpal.boxArray(),
                                                refRatio_m);
    
    // Apply mapping and boxes to *this
    geom_m = myAmrOpal.Geom();
    dmap_m = myAmrOpal.DistributionMap();
    ba_m = myAmrOpal.boxArray();
    
    bunch->gatherStatistics();
    
    
    // ------------------------------------------------------------------------
    // Initialize MultiFabs (rho, phi, efield, ...)
    // ------------------------------------------------------------------------
    initMultiFabs_m();
    
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (multi level)
    // ------------------------------------------------------------------------
    int base_level = 0;
    int fine_level = nLevels_m - 1;
    
    dynamic_cast<AmrPartBunch*>(bunch)->SetAllowParticlesNearBoundary(true);
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rho_m, base_level, 1, fine_level);
    
    double constant = -1.0 / Physics::epsilon_0;
    for (int lev = 0; lev < nLevels_m; ++lev) {
        rho_m[lev]->mult(constant, 0, 1);
    }
        

//     for (int lev = fine_level - 1 - base_level; lev >= 0; lev--)
//         BoxLib::average_down(PartMF[lev+1], PartMF[lev], 0, 1, refRatio_m[lev]);
    
    // eps in C / (V * m)
//     double constant = -1.0 / Physics::epsilon_0;
    
//     for (int lev = 0; lev < nLevels_m; lev++) {
//         PartMF[lev].mult(constant, 0, 1);
//         
//         MultiFab::Add(*rho_m[base_level+lev], PartMF[lev], 0, 0, 1, 0);
//     }
    
    Real offset = 0.0;
    Solver sol;
    sol.solve_for_accel(rho_m, phi_m, efield_m, geom_m,
                        base_level, fine_level,
                        offset);
    
    // get solution on coarsest level by averaging down from finest.
    for (int i = fine_level-1; i >= 0; --i) {
        BoxLib::average_down(*phi_m[i+1], *phi_m[i], 0, 1, refRatio_m[i]);
    }
    
    std::string dir = "plt0000";
    Real time = 0.0;
    
//     for (int lev = 0; lev < nLevels_m; lev++) {
//         rho_m[lev]->mult(1.0 / constant, 0, 1);
//     }
    
    writePlotFile(dir, rho_m, phi_m, efield_m, refRatio_m, geom_m, time);
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (single level)
    // ------------------------------------------------------------------------
    base_level = 0;
    fine_level = 0;
    
    
//     for (int lev = 0; lev < nLevels_m; lev++) {
//         rho_m[lev]->mult(constant, 0, 1);
//     }
    
    sol.solve_for_accel(rho_m, phi_single_m, efield_m, geom_m,
                        base_level, fine_level,
                        offset);
    
    
    
    // comparison
    MultiFab::Subtract(*phi_single_m[0], *phi_m[0], 0, 0, 1, 1);
    
    return phi_single_m[0]->norm2();
    
}



double PoissonProblems::doSolveParticlesReal(int step, std::string h5file) {
    
    
    // ------------------------------------------------------------------------
    // Read in particle distribution
    // ------------------------------------------------------------------------
    PartBunchBase* bunch = new AmrPartBunch(geom_m[0], dmap_m[0], ba_m[0]);
    dynamic_cast<AmrPartBunch*>(bunch)->SetAllowParticlesNearBoundary(true);
    
    Distribution dist;
    dist.readH5(h5file, step);
    dist.injectBeam(*bunch);
    
//     for (uint i = 0; i < bunch->getLocalNum(); ++i)
//         std::cout << bunch->getR(i) << std::endl;
    
    // redistribute on multi-level
    dynamic_cast<AmrPartBunch*>(bunch)->Define (geom_m, dmap_m, ba_m, refRatio_m);
    
    bunch->gatherStatistics();
    
    // ------------------------------------------------------------------------
    // Do AMR with tagging
    // ------------------------------------------------------------------------
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr_m[i];
    
    ParmParse pp("amr");
    pp.add("max_grid_size", maxGridSize_m);
    
    AmrOpal myAmrOpal(&domain_m, nLevels_m - 1, nCells, 0 /* cartesian */, bunch);
    
    dynamic_cast<AmrPartBunch*>(bunch)->Define (myAmrOpal.Geom(),
                                                myAmrOpal.DistributionMap(),
                                                myAmrOpal.boxArray(),
                                                refRatio_m);
    
    for (int i = 0; i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    // Apply mapping and boxes to *this
    geom_m = myAmrOpal.Geom();
    dmap_m = myAmrOpal.DistributionMap();
    ba_m = myAmrOpal.boxArray();
    
    
    std::cout << "max. level = " << myAmrOpal.maxLevel() << std::endl
              << "finest level = " << myAmrOpal.finestLevel() << std::endl;
    
    myAmrOpal.info();
    
    
    for (int l = 0; l < nLevels_m; ++l) {
        for (int i = 0; i < 3; ++i)
            std::cout << myAmrOpal.Geom(l).CellSize()[i] << " ";
        std::cout << std::endl;
    }
    
    
//     myAmrOpal.writePlotFile("amr0000", 1.0);
    
    bunch->gatherStatistics();
    
    
    // ------------------------------------------------------------------------
    // Initialize MultiFabs (rho, phi, efield, ...)
    // ------------------------------------------------------------------------
    initMultiFabs_m();
    
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (multi level)
    // ------------------------------------------------------------------------
    int base_level = 0;
    int fine_level = nLevels_m - 1;
    
    
//     PArray<MultiFab> PartMF;
//     PartMF.resize(nLevels_m, PArrayManage);
    
//     for (int i = 0; i < nLevels_m; ++i) {
//         PartMF.set(0, new MultiFab(ba_m[0], 1, 1));
//         PartMF[0].setVal(0.0);
//     }
    
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rho_m, base_level, 1, fine_level);

//     for (int lev = fine_level - 1 - base_level; lev >= 0; lev--)
//         BoxLib::average_down(PartMF[lev+1], PartMF[lev], 0, 1, refRatio_m[lev]);

    // eps in C / (V * m)
//     double constant = -1.0 / Physics::epsilon_0;
//     
//     for (int lev = 0; lev < nLevels_m; lev++) {
//         PartMF[lev].mult(constant, 0, 1);
//         
//         MultiFab::Add(*rho_m[base_level+lev], PartMF[lev], 0, 0, 1, 0);
//     }
    
    Real offset = 0.0;
    Solver sol;
    sol.solve_for_accel(rho_m, phi_m, efield_m, geom_m,
                        base_level, fine_level,
                        offset);
    
    // get solution on coarsest level by averaging down from finest.
    for (int i = fine_level-1; i >= 0; --i) {
        BoxLib::average_down(*phi_m[i+1], *phi_m[i], 0, 1, refRatio_m[i]);
    }
    
    std::string dir = "plt0000";
    Real time = 0.0;
    
//     for (int lev = 0; lev < nLevels_m; lev++) {
//         rho_m[lev]->mult(1.0 / constant, 0, 1);
//     }
    
    writePlotFile(dir, rho_m, phi_m, efield_m, refRatio_m, geom_m, time);
    
    
//     for (int lev = 0; lev < nLevels_m; lev++) {
//         rho_m[lev]->mult(constant, 0, 1);
//     }
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (single level)
    // ------------------------------------------------------------------------
    base_level = 0;
    fine_level = 0;
    
    
    sol.solve_for_accel(rho_m, phi_single_m, efield_m, geom_m,
                        base_level, fine_level,
                        offset);
    
    
    
    // comparison
    MultiFab::Subtract(*phi_single_m[0], *phi_m[0], 0, 0, 1, 1);
    
    return phi_single_m[0]->norm2();
}

double PoissonProblems::doSolveMultiGaussians(int nParticles, double stddev) {
    PartBunchBase* bunch = new AmrPartBunch(geom_m[0], dmap_m[0], ba_m[0]);
    
    int nloc = nParticles / ParallelDescriptor::NProcs();
    
    bunch->create(3 * nloc);
    
    double mean[3] = { -0.5, 0.0, 0.5 };
    
    std::normal_distribution<double> zdist(0.0, stddev);
    
    for (int i = 0; i < 3; ++i) {
        std::mt19937_64 mt(ParallelDescriptor::NProcs() * (i + 1));
        std::normal_distribution<double> dist(mean[i], stddev);
        
        for (int j = 0; j < nloc; ++j)
            bunch->setR( Vector_t(dist(mt), dist(mt), zdist(mt)), i * nloc + j );
    }
    
    std::cout << bunch->getLocalNum() << std::endl;
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
        bunch->setQM(1.0, i);
    
    
    dynamic_cast<AmrPartBunch*>(bunch)->Define (geom_m, dmap_m, ba_m, refRatio_m);
    
    bunch->myUpdate();
    
    // ------------------------------------------------------------------------
    // Do AMR with tagging
    // ------------------------------------------------------------------------
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr_m[i];
    
    ParmParse pp("amr");
    pp.add("max_grid_size", maxGridSize_m);
    
    AmrOpal myAmrOpal(&domain_m, nLevels_m - 1, nCells, 0 /* cartesian */, bunch);
    
    
    for (int i = 0; i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    dynamic_cast<AmrPartBunch*>(bunch)->Define (myAmrOpal.Geom(),
                                                myAmrOpal.DistributionMap(),
                                                myAmrOpal.boxArray(),
                                                refRatio_m);
    
    // Apply mapping and boxes to *this
    geom_m = myAmrOpal.Geom();
    dmap_m = myAmrOpal.DistributionMap();
    ba_m = myAmrOpal.boxArray();
    
    
    std::cout << "max. level = " << myAmrOpal.maxLevel() << std::endl
              << "finest level = " << myAmrOpal.finestLevel() << std::endl;
    
    myAmrOpal.info();
    
    
    bunch->gatherStatistics();
    
    
    // ------------------------------------------------------------------------
    // Initialize MultiFabs (rho, phi, efield, ...)
    // ------------------------------------------------------------------------
    initMultiFabs_m();
    
    
    // ------------------------------------------------------------------------
    // Solve with MultiGrid (multi level)
    // ------------------------------------------------------------------------
    int base_level = 0;
    int fine_level = nLevels_m - 1;
    
    
//     PArray<MultiFab> PartMF;
//     PartMF.resize(nLevels_m, PArrayManage);
    
//     for (int i = 0; i < nLevels_m; ++i) {
//         PartMF.set(0, new MultiFab(ba_m[0], 1, 1));
//         PartMF[0].setVal(0.0);
//     }
    
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rho_m, base_level, 1, fine_level);

//     for (int lev = fine_level - 1 - base_level; lev >= 0; lev--)
//         BoxLib::average_down(PartMF[lev+1], PartMF[lev], 0, 1, refRatio_m[lev]);

//     double constant = -1.0 / Physics::epsilon_0;
    
//     for (int lev = 0; lev < nLevels_m; lev++) {
//         PartMF[lev].mult(constant, 0, 1);
//         
//         MultiFab::Add(*rho_m[base_level+lev], PartMF[lev], 0, 0, 1, 0);
//     }
    
    Real offset = 0.0;
    Solver sol;
    sol.solve_for_accel(rho_m, phi_m, efield_m, geom_m,
                        base_level, fine_level,
                        offset);
    
    // get solution on coarsest level by averaging down from finest.
    for (int i = fine_level-1; i >= 0; --i) {
        BoxLib::average_down(*phi_m[i+1], *phi_m[i], 0, 1, refRatio_m[i]);
    }
    
    std::string dir = "plt0000";
    Real time = 0.0;
    writePlotFile(dir, rho_m, phi_m, efield_m, refRatio_m, geom_m, time);
    
    return 0.0;
}

void PoissonProblems::refineWholeDomain_m() {
    int fine = 1.0;
    for (int lev = 1; lev < nLevels_m; ++lev) {
        fine *= refRatio_m[lev - 1];
        
        IntVect refined_lo(0, 0, 0);
            
        IntVect refined_hi(nr_m[0] * fine - 1,
                           nr_m[1] * fine - 1,
                           nr_m[2] * fine - 1);
        
        Box refined_patch(refined_lo, refined_hi);
        ba_m[lev].define(refined_patch);
        
        ba_m[lev].maxSize(maxGridSize_m);
        dmap_m[lev].define(ba_m[lev], ParallelDescriptor::NProcs() /*nprocs*/);
    }
}


void PoissonProblems::initMultiFabs_m() {
    rho_m.resize(nLevels_m);
    phi_m.resize(nLevels_m);
    efield_m.resize(nLevels_m);
    phi_single_m.resize(nLevels_m);   // single level solve
    
    for (int l = 0; l < nLevels_m; ++l) {
        rho_m[l] = std::unique_ptr<MultiFab>(new MultiFab(ba_m[l], 1, 0));
        rho_m[l]->setVal(0.0);
        
        phi_m[l] = std::unique_ptr<MultiFab>(new MultiFab(ba_m[l], 1, 1));
        phi_m[l]->setVal(0.0);
        
        phi_single_m[l] = std::unique_ptr<MultiFab>(new MultiFab(ba_m[l], 1, 1));
        phi_single_m[l]->setVal(0.0);
        
        efield_m[l] = std::unique_ptr<MultiFab>(new MultiFab(ba_m[l], BL_SPACEDIM, 1));
        efield_m[l]->setVal(0.0);
    }
}