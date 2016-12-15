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

#include "PartBunch.h"
#include "AmrPartBunch.h"

#include "Solver.h"
#include "AmrOpal.h"

#include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"
#include <random>


typedef Array<std::unique_ptr<MultiFab> > container_t;

void initSphere(double r, PartBunchBase* bunch, int nParticles) {
    bunch->create(nParticles);
    
    std::mt19937_64 eng;
    std::uniform_real_distribution<> d(-r * 0.5, r * 0.5);
    
    std::string outfile = "amr-particles-level-" + std::to_string(0);
    std::ofstream out(outfile);
    long double qi = 4.0 * Physics::pi * Physics::epsilon_0 * r * r / double(nParticles);
    
    std::cout << "At R = " << r << ": " << qi * nParticles / (4.0 * Physics::pi * Physics::epsilon_0 * r) << std::endl;
    
    for (int i = 0; i < nParticles; ++i) {
        bunch->setR(Vector_t( d(eng), d(eng), d(eng)), i);    // m
        bunch->setQM(qi, i);   // C
        
        out << bunch->getR(i)[0] << std::endl;
    }
    out.close();
}


void writePotential(const container_t& phi, double dx, double dlow) {
    // potential and efield a long x-direction (y = 0, z = 0)
    int lidx = 0;
    for (MFIter mfi(*phi[0 /*level*/]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        FArrayBox& lhs = (*phi[0])[mfi];
        /* Write the potential of all levels to file in the format: x, y, z, \phi
         */
        for (int proc = 0; proc < ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == ParallelDescriptor::MyProc() ) {
                std::string outfile = "amr-phi_scalar-level-" + std::to_string(0);
                std::ofstream out;
                
                if ( proc == 0 )
                    out.open(outfile);
                else
                    out.open(outfile, std::ios_base::app);
                
//                 if ( !out.is_open() )
//                     throw OpalException("Error in TrilinosSolver::CopySolution",
//                                         "Couldn't open the file: " + outfile);
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
    for (MFIter mfi(*efield[0 /*level*/]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        FArrayBox& lhs = (*efield[0])[mfi];
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
                
//                 if ( !out.is_open() )
//                     throw OpalException("Error in TrilinosSolver::CopySolution",
//                                         "Couldn't open the file: " + outfile);
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


void doSolve(AmrOpal& myAmrOpal, PartBunchBase* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Array<Geometry>& geom,
             const Array<int>& rr, int nLevels)
{
    static IpplTimings::TimerRef solvTimer = IpplTimings::getTimer("solv");
    // =======================================================================                                                                                                                                   
    // 4. prepare for multi-level solve                                                                                                                                                                          
    // =======================================================================

    rhs.resize(nLevels);
    phi.resize(nLevels);
    grad_phi.resize(nLevels);

    for (int lev = 0; lev < nLevels; ++lev) {
        //                                    # component # ghost cells                                                                                                                                          
        rhs[lev] = std::unique_ptr<MultiFab>(new MultiFab(myAmrOpal.boxArray()[lev],1          ,0));
        phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(myAmrOpal.boxArray()[lev],1          ,1));
        grad_phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(myAmrOpal.boxArray()[lev],BL_SPACEDIM, 0));

        rhs[lev]->setVal(0.0);
        phi[lev]->setVal(0.0);
        grad_phi[lev]->setVal(0.0);
    }
    

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();

//     const Box& bx = rhs[0]->box(0);
//     
//     int cell[3] = {
//         (*bx.hiVect() - *bx.loVect()) / 2,
//         (*bx.hiVect() - *bx.loVect()) / 2,
//         (*bx.hiVect() - *bx.loVect()) / 2
//     };
//     
//     IntVect low(cell[0] - 1,
//                 cell[1] - 1,
//                 cell[2] - 1);
//     
//     IntVect high(cell[0] + 1,
//                  cell[1] + 1,
//                  cell[2] + 1);
//     
//     Box region(low, high);
    
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rhs, base_level, 1, finest_level);
    
    Real totalCharge = 0.0;
    for (int i = 0; i <= finest_level; ++i) {
        Real invVol = (*(geom[i].CellSize()) * *(geom[i].CellSize()) * *(geom[i].CellSize()) );
        Real sum = rhs[i]->sum(0) * invVol;
        std::cout << "level " << i << " sum = " << sum << std::endl;
        totalCharge += sum;
    }
    
    std::cout << rhs[0]->min(0) << " " << rhs[0]->max(0) << std::endl;
    
    std::cout << "Total Charge: " << totalCharge << " C" << std::endl;
    std::cout << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << std::endl;
    Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    std::cout << "Cell volume: " << *(geom[0].CellSize()) << "^3 = " << vol << " m^3" << std::endl;
    
    // eps in C / (V * m)
    double constant = -1.0/*Physics::q_e*/ / (4.0 * Physics::pi * Physics::epsilon_0);  // in [V m / C]
    for (int i = 0; i <=finest_level; ++i) {
        rhs[i]->mult(constant, 0, 1);       // in [V m]
    }
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;
//     if (geom[0].isAllPeriodic()) 
//     {
//         offset = dynamic_cast<AmrPartBunch*>(bunch)->sumParticleMass(0, 0);
//         offset /= geom[0].ProbSize();
//         std::cout << "offset = " << offset << std::endl;
//     }

    // solve                                                                                                                                                                                                     
    Solver sol;
    IpplTimings::startTimer(solvTimer);
    sol.solve_for_accel(rhs,            // [V m]
                        phi,            // [V m^3]
                        grad_phi,       // [V m^2]
                        geom,
                        base_level,
                        finest_level,
                        offset);
    
    for (int i = 0; i <=finest_level; ++i) {
        rhs[i]->mult(1.0 / constant, 0, 1);
    }
    
    IpplTimings::stopTimer(solvTimer);
}

void doBoxLib(const Vektor<size_t, 3>& nr,
              int nLevels, size_t maxBoxSize, Inform& msg, Inform& msg2all)
{
    static IpplTimings::TimerRef distTimer = IpplTimings::getTimer("dist");
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
    
    // box [-1,1]x[-1,1]x[-1,1]
    RealBox domain;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain.setLo(i, -1.0);
        domain.setHi(i,  1.0);
    }
    
    domain.setLo(0, -0.05); // m
    domain.setHi(0,  0.05); // m
    domain.setLo(1, -0.05); // m
    domain.setHi(1,  0.05); // m
    domain.setLo(2, -0.05); // m
    domain.setHi(2,  0.05); // m
    
    // periodic boundary conditions in all directions
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
    
    
    Array<int> rr(nLevels - 1);
    for (unsigned int lev = 0; lev < rr.size(); ++lev)
        rr[lev] = 2;
    
    // geometries of refined levels
    for (int lev = 1; lev < nLevels; ++lev) {
        geom[lev].define(BoxLib::refine(geom[lev - 1].Domain(),
                                        rr[lev - 1]),
                         &domain, 0, bc);
    }
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    PartBunchBase* bunch = new AmrPartBunch(geom[0], dmap[0], ba[0]);
    
    
    // initialize a particle distribution
    IpplTimings::startTimer(distTimer);
    
    double R = 0.01; // m
    int nParticles = 100000;
    initSphere(R, bunch, nParticles);
    
//     bunch->setR(Vector_t(0.0, 0.0, 0.0), 1);
//     bunch->setQM(0.0, 1);
    
    msg << "Particle location: " <<  bunch->getR(0) << " m" << endl
        << "Particle charge: " << bunch->getQM(0) << " C" << endl
        << "#Cells per dim for bunch: " << R / *(geom[0].CellSize()) << endl;
    
    // redistribute on single-level
    bunch->myUpdate();
    
//     bunch->gatherStatistics();
    
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
    doSolve(myAmrOpal, bunch, rhs, phi, grad_phi, geom, rr, nLevels);
    
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << grad_phi[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << grad_phi[i]->min(0) << endl;
    }
    
//     msg << "Potential: " << dynamic_cast<AmrPartBunch*>(bunch)->interpolate(0, *phi[0]) << endl;
    writePotential(phi, *(geom[0].CellSize()), -0.05);
    writeElectricField(grad_phi, *(geom[0].CellSize()), -0.05);
    
    writePlotFile(plotsolve, rhs, phi, grad_phi, rr, geom, 0);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    Inform msg2all("Solver", INFORM_ALL_NODES);
    

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
    doBoxLib(nr, nLevels, maxBoxSize, msg, msg2all);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    return 0;
}
