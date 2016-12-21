/*
 * @file testGridSolve.cpp
 * @author Matthias Frey
 * @date 20. Dec. 2016
 * @brief Solve the electrostatic potential for a cube
 *        with Dirichlet boundary condition.
 * @details In this example we put -1.0 on every grid
 *          point (cell-centered) for the right-hand side
 *          of the Poisson equation. The result can be compared
 *          with the iterative.cpp using the
 *          initMinusOneEverywhere() function for the right-hand
 *          side initialization.\n
 *          Domain: [-0.05 (m), 0.05 (m)]^3
 */
#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>


#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>
#include <ParmParse.H>

#include "Solver.h"
#include "AmrOpal.h"

#include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"
#include <random>

#ifdef UNIQUE_PTR
typedef Array<std::unique_ptr<MultiFab> > container_t;
#else
#include <PArray.H>
typedef PArray<MultiFab> container_t;
#endif


void writePotential(const container_t& phi, double dx, double dlow, std::string filename) {
    // potential and efield a long x-direction (y = 0, z = 0)
    int lidx = 0;
#ifdef UNIQUE_PTR
    for (MFIter mfi(*phi[0 /*level*/]); mfi.isValid(); ++mfi) {
#else
    for (MFIter mfi(phi[0 /*level*/]); mfi.isValid(); ++mfi) {
#endif
        const Box& bx = mfi.validbox();
#ifdef UNIQUE_PTR
        const FArrayBox& lhs = (*phi[0])[mfi];
#else
        const FArrayBox& lhs = (phi[0])[mfi];
#endif
        /* Write the potential of all levels to file in the format: x, y, z, \phi
         */
        for (int proc = 0; proc < ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == ParallelDescriptor::MyProc() ) {
                std::string outfile = filename + std::to_string(0);
                std::ofstream out;
                
                if ( proc == 0 )
                    out.open(outfile);
                else
                    out.open(outfile, std::ios_base::app);
                
                int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
                int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
//                     for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
//                         for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                            IntVect ivec(i, j, k);
                            // add one in order to have same convention as PartBunch::computeSelfField()
                            out << i + 1 << " " << j + 1 << " " << k + 1 << " " << dlow + dx * i << " "
                                << lhs(ivec, 0)  << std::endl;
//                         }
//                     }
                }
                out.close();
            }
            ParallelDescriptor::Barrier();
        }
    }
}

void writeElectricField(const container_t& efield, double dx, double dlow) {
    // potential and efield a long x-direction (y = 0, z = 0)
    int lidx = 0;
#ifdef UNIQUE_PTR
    for (MFIter mfi(*efield[0 /*level*/]); mfi.isValid(); ++mfi) {
#else
    for (MFIter mfi(efield[0 /*level*/]); mfi.isValid(); ++mfi) {
#endif
        const Box& bx = mfi.validbox();
#ifdef UNIQUE_PTR
        const FArrayBox& lhs = (*efield[0])[mfi];
#else
        const FArrayBox& lhs = (efield[0])[mfi];
#endif
        /* Write the potential of all levels to file in the format: x, y, z, ex, ey, ez
         */
        for (int proc = 0; proc < ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == ParallelDescriptor::MyProc() ) {
                std::string outfile = "amr-efield_scalar-level-" + std::to_string(0);
                std::ofstream out;
                
                if ( proc == 0 )
                    out.open(outfile);
                else
                    out.open(outfile, std::ios_base::app);
                
                int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
                int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
//                     for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
//                         for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                            IntVect ivec(i, j, k);
                            // add one in order to have same convention as PartBunch::computeSelfField()
                            out << i + 1 << " " << j + 1 << " " << k + 1 << " " << dlow + dx * i << " "
                                << lhs(ivec, 0) << " " << lhs(ivec, 1) << " " << lhs(ivec, 2) << std::endl;
//                         }
//                     }
                }
                out.close();
            }
            ParallelDescriptor::Barrier();
        }
    }
}


void doSolve(const Array<BoxArray>& ba,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Array<Geometry>& geom,
             const Array<int>& rr, int nLevels,
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
        rhs[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev],1          ,0));
        phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev],1          ,1));
        grad_phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev],BL_SPACEDIM, 1));

        phi[lev]->setVal(0.0);
        grad_phi[lev]->setVal(0.0);
        
        rhs[lev]->setVal(-1.0);
#else
        //                                    # component # ghost cells                                                                                                                                          
        rhs.set(lev, new MultiFab(ba[lev],1          ,0));
        phi.set(lev, new MultiFab(ba[lev],1          ,1));
        grad_phi.set(lev, new MultiFab(ba[lev],BL_SPACEDIM, 1));
        phi[lev].setVal(0.0);
        grad_phi[lev].setVal(0.0);
        
        rhs[lev].setVal(-1.0);
#endif
    }
    
    writePotential(rhs, *(geom[0].CellSize()), 0, "amr-rho_scalar-level-");
    

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = 0;
    
    // Check charge conservation
    Real totalCharge = 0.0;
    for (int i = 0; i <= finest_level; ++i) {
#ifdef UNIQUE_PTR
        Real sum = rhs[i]->sum(0);
#else
        Real sum = rhs[i].sum(0);
#endif
        totalCharge += sum;
    }
    
    msg << "Total Charge (computed): " << totalCharge << " C" << endl
        << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    
    Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << *(geom[0].CellSize()) << "^3 = " << vol << " m^3" << endl;
    
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
    
    /*
     * nLevel is the number of levels allowed, i.e if nLevel = 1
     * we just run single-level
     */
    
    /*
     * set up the geometry
     */
    IntVect low(0, 0, 0);
    IntVect high(nr[0] - 1, nr[1] - 1, nr[2] - 1);    
    Box bx(low, high);
    
    // box [0, 0.05] x [0, 0.05] x [0, 0.05]
    RealBox domain;
    double a = 0.1;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain.setLo(i,  0); // m
        domain.setHi(i,  a); // m
    }
    
    // Dirichlet boundary conditions in all directions
    int bc[BL_SPACEDIM] = {0, 0, 0};
    
    
    Array<Geometry> geom;
    geom.resize(nLevels);
    
    // level 0 describes physical domain
    geom[0].define(bx, &domain, 0, bc);
    
    // Container for boxes at all levels
    Array<BoxArray> ba;
    ba.resize(nLevels);
    
    // box at level 0
    ba[0].define(bx);
    ba[0].maxSize(maxBoxSize);
    
    /*
     * distribution mapping
     */
    Array<DistributionMapping> dmap;
    dmap.resize(nLevels);
    dmap[0].define(ba[0], ParallelDescriptor::NProcs() /*nprocs*/);
    
    // refinement ratio
    Array<int> rr(nLevels - 1);
    for (unsigned int lev = 0; lev < rr.size(); ++lev)
        rr[lev] = 2;
    
    // geometries of refined levels
    for (int lev = 1; lev < nLevels; ++lev) {
        geom[lev].define(BoxLib::refine(geom[lev - 1].Domain(),
                                        rr[lev - 1]),
                         &domain, 0, bc);
    }
    
    msg << "Charge per grid point: 1.0 C" << endl
        << "Total charge: " << 1.0 * nr[0] * nr[1] * nr[2] << " C" << endl;
    
    container_t rhs;
    container_t phi;
    container_t grad_phi;
    
    std::string plotsolve = BoxLib::Concatenate("plt", 0, 4);
    doSolve(ba, rhs, phi, grad_phi, geom, rr, nLevels, msg);
    
    
    for (int i = 0; i < nLevels; ++i) {
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
    
    writePotential(phi, *(geom[0].CellSize()), 0, "amr-phi_scalar-level-");
    writeElectricField(grad_phi, *(geom[0].CellSize()), 0);
    
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
        
    BoxLib::Initialize(argc,argv, false);
    size_t nLevels = std::atoi(argv[4]) + 1; // i.e. nLevels = 0 --> only single level
    size_t maxBoxSize = std::atoi(argv[5]);
    doBoxLib(nr, nLevels, maxBoxSize, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    return 0;
}
