/*
 * @file testUnifSphereGrid.cpp
 * @author Matthias Frey
 * @date 20. Dec. 2016
 * @brief Solve the electrostatic potential for a cube
 *        with Dirichlet boundary condition.
 * @details In this example we put -1.0 on every grid
 *          point (cell-centered) within a sphere of radius
 *          R = 0.005 [m] for the right-hand side of the Poisson equation.
 *          The result can be compared
 *          with the iterative.cpp using the
 *          initMinusOneEverywhere() function for the right-hand
 *          side initialization. \n
 *          Domain: [-0.05 (m), 0.05 (m)]^3
 */
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

#include "helper_functions.h"

#include "Solver.h"
#include "AmrOpal.h"

#include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"
#include <random>

void initSphereOnGrid(container_t& rhs,
                      const Array<Geometry>& geom,
                      double a,
                      double R,
                      const Vektor<size_t, 3>& nr)
{
    // number of non-zero values
    int nnz = 0;
    
    // total charge
    long double qi = /*4.0 * Physics::pi * Physics::epsilon_0 **/ R * R;
    
#ifdef UNIQUE_PTR
    rhs[0]->setVal(0.0);
    for (MFIter mfi(*rhs[0 /*level*/]); mfi.isValid(); ++mfi) {
#else
    rhs[0].setVal(0.0);
    for (MFIter mfi(rhs[0 /*level*/]); mfi.isValid(); ++mfi) {
#endif
        const Box& bx = mfi.validbox();
#ifdef UNIQUE_PTR
        FArrayBox& rho = (*rhs[0])[mfi];
#else
        FArrayBox& rho = (rhs[0])[mfi];
#endif
        for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
            for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                    double x = 2.0 * a / double(nr[0] - 1) * i - a;
                    double y = 2.0 * a / double(nr[1] - 1) * j - a;
                    double z = 2.0 * a / double(nr[2] - 1) * k - a;
                    
                    if ( x * x + y * y + z * z <= R * R ) {
                        IntVect idx(i, j, k);
                        rho(idx, 0) = qi;//-1.0;;
                        ++nnz;
                    }
                }
            }
        }
    }
    
    std::cout << "#Non-zero grid points: " << nnz << std::endl;
    std::cout << "Total Charge: " << qi << std::endl;
    
#ifdef UNIQUE_PTR
        rhs[0]->mult(1.0 / double(nnz /**  * Physics::epsilon_0*/), 0, 1);
#else
        rhs[0].mult(1.0 / double(nnz /**  * Physics::epsilon_0*/), 0, 1);
#endif
}

void doSolve(AmrOpal& myAmrOpal, PartBunchBase* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Array<Geometry>& geom,
             const Array<int>& rr, int nLevels,
             Inform& msg,
             double a,
             double R,
             const Vektor<size_t, 3>& nr,
             bool grid = false)
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
    
    
    initSphereOnGrid(rhs, geom, a, R, nr);
    
    writePotential(rhs, *(geom[0].CellSize()), -a, "amr-rho_scalar-level-");
    
    // Check charge conservation
    Real totalCharge = 0.0;
    for (int i = 0; i <= finest_level; ++i) {
        Real invVol = (*(geom[i].CellSize()) * *(geom[i].CellSize()) * *(geom[i].CellSize()) );
#ifdef UNIQUE_PTR
        Real sum = rhs[i]->sum(0);// * invVol;
#else
        Real sum = rhs[i].sum(0);// * invVol;
#endif
        totalCharge += sum;
    }
    
    msg << "Total Charge (computed): " << totalCharge << " C" << endl
        << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    
    Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << *(geom[0].CellSize()) << "^3 = " << vol << " m^3" << endl;
    
    // eps in C / (V * m)
    double constant = -1.0 /* / Physics::epsilon_0*/;  // in [V m / C]
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
    
//     for (int i = 0; i <=finest_level; ++i) {
//         rhs[i]->mult(1.0 / constant, 0, 1);
//     }
    
    IpplTimings::stopTimer(solvTimer);
}

void doBoxLib(const Vektor<size_t, 3>& nr,
              int nLevels,
              size_t maxBoxSize,
              Inform& msg)
{
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    double lower = -0.05; // m
    double upper =  0.05; // m
    
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
    double R = 0.005; // radius of sphere [m]
    
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
    
    Array<int> nCells(BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; ++i)
        nCells[i] = nr[i];
    
    AmrOpal myAmrOpal(&domain, nLevels - 1, nCells, 0 /* cartesian */, bunch);
    
    container_t rhs;
    container_t phi;
    container_t grad_phi;
    
    std::string plotsolve = BoxLib::Concatenate("plt", 0, 4);
    doSolve(myAmrOpal, bunch, rhs, phi, grad_phi, geom, rr, nLevels, msg, upper, R, nr, true);
    
    
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
    
    writePotential(phi, *(geom[0].CellSize()), lower, "amr-phi_scalar-level-");
    writeElectricField(grad_phi, *(geom[0].CellSize()), lower);
    
    writePlotFile(plotsolve, rhs, phi, grad_phi, rr, geom, 0);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] testGaussian"
         << "[#gridpoints x] [#gridpoints y] [#gridpoints z] [#levels] [max. box size]";
    
    if ( argc < 6 ) {
        msg << call.str() << endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<size_t, 3> nr(std::atoi(argv[1]),
                         std::atoi(argv[2]),
                         std::atoi(argv[3]));
    
    
    msg << "Particle test running with" << endl
        << "- grid       = " << nr << endl;
        
    BoxLib::Initialize(argc,argv, false);
    size_t nLevels = std::atoi(argv[4]) + 1; // i.e. nLevels = 0 --> only single level
    size_t maxBoxSize = std::atoi(argv[5]);
    doBoxLib(nr, nLevels, maxBoxSize, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    return 0;
}
