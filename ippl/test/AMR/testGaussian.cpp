// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 *
 ***************************************************************************/

/***************************************************************************

This test program sets up a simple sine-wave electric field in 3D,
  creates a population of particles with random q/m values (charge-to-mass
  ratio) and velocities, and then tracks their motions in the static
  electric field using nearest-grid-point interpolation.

Usage:
 mpirun -np 4 testGaussian [#gridpoints x] [#gridpoints y] [#gridpoints z]
                           [#particles] [#levels] [max. box size]

***************************************************************************/

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>

#include <ParmParse.H>

#include "PartBunch.h"
#include "AmrPartBunch.h"

#include "Distribution.h"
#include "Solver.h"
#include "AmrOpal.h"

#include "helper_functions.h"

#include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

void doSolve(AmrOpal& myAmrOpal, PartBunchBase* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             int nLevels,
             Inform& msg)
{
    static IpplTimings::TimerRef solvTimer = IpplTimings::getTimer("solv");
    // =======================================================================                                                                                                                                   
    // 4. prepare for multi-level solve                                                                                                                                                                          
    // =======================================================================

    rhs.resize(nLevels);
    phi.resize(nLevels);
    grad_phi.resize(nLevels);

    for (int lev = 0; lev < nLevels; ++lev) {
#ifdef UNIQUE_PTR
        //                                    # component # ghost cells                                                                                                                                          
        rhs[lev] = std::unique_ptr<MultiFab>(new MultiFab(myAmrOpal.boxArray()[lev],1          ,0));
        phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(myAmrOpal.boxArray()[lev],1          ,1));
        grad_phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(myAmrOpal.boxArray()[lev],BL_SPACEDIM, 1));

        rhs[lev]->setVal(0.0);
        phi[lev]->setVal(0.0);
        grad_phi[lev]->setVal(0.0);
#else
        //                                    # component # ghost cells                                                                                                                                          
        rhs.set(lev, new MultiFab(myAmrOpal.boxArray()[lev],1          ,0));
        phi.set(lev, new MultiFab(myAmrOpal.boxArray()[lev],1          ,1));
        grad_phi.set(lev, new MultiFab(myAmrOpal.boxArray()[lev],BL_SPACEDIM, 1));
        
        rhs[lev].setVal(0.0);
        phi[lev].setVal(0.0);
        grad_phi[lev].setVal(0.0);
#endif
    }
    

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();

    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rhs, base_level, 1, finest_level);
    
    Real totalCharge = 0.0;
    for (int i = 0; i <= finest_level; ++i) {
        Real vol = (*(geom[i].CellSize()) * *(geom[i].CellSize()) * *(geom[i].CellSize()) );
        
#ifdef UNIQUE_PTR
        Real sum = rhs[i]->sum(0) * vol;
#else
        Real sum = rhs[i].sum(0) * vol;
#endif
        totalCharge += sum;
    }
    msg << "Total Charge: " << totalCharge << " C" << endl
        << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << vol << " m^3" << endl;
    
    // eps in C / (V * m)
    double constant = -1.0 / Physics::epsilon_0;  // in [V m / C]
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(constant, 0, 1);
#endif
    }
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;

    // solve                                                                                                                                                                                                     
    Solver sol;
    IpplTimings::startTimer(solvTimer);
    sol.solve_for_accel(rhs,
                        phi,
                        grad_phi,
                        geom,
                        base_level,
                        finest_level,
                        offset);
    
    // for plotting unnormalize
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(1.0 / constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(1.0 / constant, 0, 1);
#endif
    }
    
    IpplTimings::stopTimer(solvTimer);
}

void doBoxLib(const Vektor<size_t, 3>& nr, size_t nParticles,
              int nLevels, size_t maxBoxSize, Inform& msg)
{
    static IpplTimings::TimerRef distTimer = IpplTimings::getTimer("dist");
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    double lower = -0.5; // m
    double upper =  0.5; // m
    
    RealBox domain;
    Array<BoxArray> ba;
    Array<Geometry> geom;
    Array<DistributionMapping> dmap;
    Array<int> rr;
    
    // in helper_functions.h
    init(domain, ba, dmap, geom, rr, nr, nLevels, maxBoxSize, lower, upper);
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    PartBunchBase* bunch = new AmrPartBunch(geom[0], dmap[0], ba[0]);
    
    
    // initialize a particle distribution
    unsigned long int nloc = nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    IpplTimings::startTimer(distTimer);
    double sig = 0.01;  // m
    dist.gaussian(0.0, sig, nloc, ParallelDescriptor::MyProc());
    
    // copy particles to the PartBunchBase object.
    dist.injectBeam(*bunch);
    
    
    
    
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
        bunch->setQM(Physics::q_e, i);  // in [C]
    
    // redistribute on single-level
    bunch->myUpdate();
    
    bunch->gatherStatistics();
    
    msg << "Bunch radius: " << sig << " m" << endl
        << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->getQM(0) << " C" << endl
        << "Total charge: " << nParticles * bunch->getQM(0) << " C" << endl;
    
    for (int i = 0; i < nLevels; ++i)
        msg << "#Cells per dim of level " << i << " for bunch : "
            << 2.0 * sig / *(geom[i].CellSize()) << endl;
    
    // ========================================================================
    // 2. tagging (i.e. create BoxArray's, DistributionMapping's of all
    //    other levels)
    // ========================================================================
    
    
    /*
     * create an Amr object
     */
    ParmParse pp("amr");
    pp.add("max_grid_size", int(maxBoxSize));
    
    Array<int> error_buf(nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr[i];
    
    AmrOpal myAmrOpal(&domain, nLevels - 1, nCells, 0 /* cartesian */, bunch);
    
    /*
     * do tagging
     */
    dynamic_cast<AmrPartBunch*>(bunch)->Define (myAmrOpal.Geom(),
                                                myAmrOpal.DistributionMap(),
                                                myAmrOpal.boxArray(),
                                                rr);
    
    
    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    container_t rhs;
    container_t phi;
    container_t grad_phi;
    
    std::string plotsolve = BoxLib::Concatenate("plt", 0, 4);
    doSolve(myAmrOpal, bunch, rhs, phi, grad_phi, geom, rr, nLevels, msg);
    
    msg << "Total field energy: " << totalFieldEnergy(grad_phi, rr) << endl;
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
#ifdef UNIQUE_PTR
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << grad_phi[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << grad_phi[i]->min(0) << endl;
#else
        msg << "Max. potential level " << i << ": "<< phi[i].max(0) << endl
            << "Min. potential level " << i << ": " << phi[i].min(0) << endl
            << "Max. ex-field level " << i << ": " << grad_phi[i].max(0) << endl
            << "Min. ex-field level " << i << ": " << grad_phi[i].min(0) << endl;
#endif
    }
    
    
    writePlotFile(plotsolve, rhs, phi, grad_phi, rr, geom, 0);
    
    dynamic_cast<AmrPartBunch*>(bunch)->python_format(0);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] testGaussian"
         << "[#gridpoints x] [#gridpoints y] [#gridpoints z] [#particles] "
         << "[#levels] [max. box size]";
    
    if ( argc < 7 ) {
        msg << call.str() << endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<size_t, 3> nr(std::atoi(argv[1]),
                         std::atoi(argv[2]),
                         std::atoi(argv[3]));
    
    
    size_t nParticles = std::atoi(argv[4]);
    
    
    msg << "Particle test running with" << endl
        << "- #particles = " << nParticles << endl
        << "- grid       = " << nr << endl;
        
    BoxLib::Initialize(argc,argv, false);
    size_t nLevels = std::atoi(argv[5]) + 1; // i.e. nLevels = 0 --> only single level
    size_t maxBoxSize = std::atoi(argv[6]);
    doBoxLib(nr, nParticles, nLevels, maxBoxSize, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    return 0;
}
