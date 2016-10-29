#ifndef POISSONPROBLEMS_H
#define POISSONPROBLEMS_H

#include <BoxLib.H>
#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>

#include "Solver.h"
#include "AmrPartBunch.h"
#include "Distribution.h"

#include "writePlotFile.H"

class PoissonProblems {
    
public:
    PoissonProblems(int nr[3], int maxGridSize, int nLevels);
    
    
    // solve \Lap\phi = - 1
    double doSolveNoParticles();
    
    /* solve \Lap\phi = -1
     * where we use particles and assign density
     */
    double doSolveParticlesUniform();
    
    double doSolveParticlesGaussian(int nParticles);
    
    
    double doSolveParticlesReal(int step, std::string h5file);
    
private:
    void refineWholeDomain_m();
    void initMultiFabs_m();
    
private:
    int nr_m[3];
    int maxGridSize_m;
    int nLevels_m;
    Array<Geometry> geom_m;
    Array<DistributionMapping> dmap_m;
    Array<BoxArray> ba_m;
    Array<int> refRatio_m;
    
    PArray<MultiFab> rho_m;
    PArray<MultiFab> phi_m;
    PArray<MultiFab> efield_m;
    PArray<MultiFab> phi_single_m;
    
};

#endif