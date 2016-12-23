/*
 * @file testUnifSphere.cpp
 * @author Matthias Frey
 * @date 20. Dec. 2016
 * @brief Solve the electrostatic potential for a cube
 *        with Dirichlet boundary condition.
 * @details In this example we initialize 1e6 particles
 *          within a sphere of radius R = 0.005 [m].
 *          We then solve the Poisson equation.\n
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

#include "Solver.h"
#include "AmrOpal.h"

#include "helper_functions.h"

#include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"
#include <random>

void initSphere(double r, PartBunchBase* bunch, int nParticles) {
    bunch->create(nParticles);
    
    std::mt19937_64 eng;
    std::uniform_real_distribution<> ph(-1.0, 1.0);
    std::uniform_real_distribution<> th(0.0, 2.0 * Physics::pi);
    std::uniform_real_distribution<> u(0.0, 1.0);
    
    
    std::string outfile = "amr-particles-level-" + std::to_string(0);
    std::ofstream out(outfile);
    
    long double qi = 4.0 * Physics::pi * Physics::epsilon_0 * r * r / double(nParticles);
    
    for (int i = 0; i < nParticles; ++i) {
        // 17. Dec. 2016,
        // http://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
        // http://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume
        double phi = std::acos( ph(eng) );
        double theta = th(eng);
        double radius = r * std::cbrt( u(eng) );
        
        double x = radius * std::cos( theta ) * std::sin( phi );
        double y = radius * std::sin( theta ) * std::sin( phi );
        double z = radius * std::cos( phi );
        
        out << x << " " << y << " " << z << std::endl;
        
        bunch->setR(Vector_t( x, y, z ), i);    // m
        bunch->setQM(qi, i);   // C
    }
    out.close();
}

void doSolve(AmrOpal& myAmrOpal, PartBunchBase* bunch,
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
        initGridData(rhs, phi, grad_phi, myAmrOpal.boxArray()[lev], lev);
    }
    

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();
    
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rhs, base_level, 1, finest_level);
    
    writePotential(rhs, *(geom[0].CellSize()), -0.05, "amr-rho_scalar-level-");
    
    // Check charge conservation
    Real totalCharge = 0.0;
    for (int i = 0; i <= finest_level; ++i) {
        Real invVol = (*(geom[i].CellSize()) * *(geom[i].CellSize()) * *(geom[i].CellSize()) );
#ifdef UNIQUE_PTR
        Real sum = rhs[i]->sum(0) * invVol;
#else
        Real sum = rhs[i].sum(0) * invVol;
#endif
        totalCharge += sum;
    }
    
    msg << "Total Charge (computed): " << totalCharge << " C" << endl
        << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    
    Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << *(geom[0].CellSize()) << "^3 = " << vol << " m^3" << endl;
    
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
    sol.solve_for_accel(rhs,            // [V m]
                        phi,            // [V m^3]
                        grad_phi,       // [V m^2]
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
    int nParticles = 1e6;
    initSphere(R, bunch, nParticles);
    
    msg << "Bunch radius: " << R << " m" << endl
        << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->getQM(0) << " C" << endl
        << "Total charge: " << nParticles * bunch->getQM(0) << " C" << endl
        << "#Cells per dim for bunch: " << 2.0 * R / *(geom[0].CellSize()) << endl;
    
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
    
    Array<int> nCells(BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; ++i)
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
