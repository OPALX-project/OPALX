/*!
 * @file testGaussian.cpp
 * @author Matthias Frey
 * @date May 2017
 * @details Compute \f$\Delta\phi = -\rho\f$ where the charge
 * density is a Gaussian distribution. Write plot files
 * that can be visualized by yt (python visualize.py)
 * or by AmrVis of the CCSE group at LBNL.
 * 
 * Domain:  [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5]\n
 * BC:      Dirichlet boundary conditions\n
 * Charge/particle: elementary charge\n
 * Gaussian particle distribution N(0.0, 0.01)
 * 
 * Call:\n
 *  mpirun -np [#cores] testGaussian [#gridpoints x] [#gridpoints y] [#gridpoints z]
 *                                   [#particles] [#levels] [max. box size]
 *                                   [out: timing file name (optiona)]
 * 
 * @brief Computes \f$\Delta\phi = -\rho\f$
 */

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>
#include <iomanip>


#include <AMReX_ParmParse.H>

#include "AmrParticleBase.h"
#include "ParticleAmrLayout.h"
#include "PartBunchAmr.h"

#include "../Distribution.h"
#include "../Solver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

// #include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

#define Dim 3

using namespace amrex;

typedef Vektor<double, AMREX_SPACEDIM> Vector_t;

typedef ParticleAmrLayout<double,Dim> amrplayout_t;
typedef AmrParticleBase<amrplayout_t> amrbase_t;
typedef PartBunchAmr<amrplayout_t> amrbunch_t;

void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& efield,
             const Vector<Geometry>& geom,
             const Vector<int>& rr,
             int nLevels,
             Inform& msg)
{
    static IpplTimings::TimerRef allocTimer = IpplTimings::getTimer("alloc-memory-grid");
    
    static IpplTimings::TimerRef assignTimer = IpplTimings::getTimer("assign-charge");
    
    // =======================================================================                                                                                                                                   
    // 4. prepare for multi-level solve                                                                                                                                                                          
    // =======================================================================
    
    rhs.resize(nLevels);
    phi.resize(nLevels);
    efield.resize(nLevels);
    
    IpplTimings::startTimer(allocTimer);
    
    for (int lev = 0; lev < nLevels; ++lev) {
        initGridData(rhs, phi, efield, myAmrOpal.boxArray()[lev],
                     myAmrOpal.DistributionMap(lev), lev);
    }
    
    IpplTimings::stopTimer(allocTimer);

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();
    
    Vector<std::unique_ptr<MultiFab> > partMF(nLevels);
    for (int lev = 0; lev < nLevels; lev++) {
        const BoxArray& ba = rhs[lev]->boxArray();
        const DistributionMapping& dmap = rhs[lev]->DistributionMap();
        partMF[lev].reset(new MultiFab(ba, dmap, 1, 2));
        partMF[lev]->setVal(0.0, 2);
    }

    IpplTimings::startTimer(assignTimer);
    bunch->AssignDensityFort(bunch->qm, partMF, base_level, 1, finest_level);
    IpplTimings::stopTimer(assignTimer);
    
    for (int lev = 0; lev < nLevels; ++lev) {
        MultiFab::Copy(*rhs[lev], *partMF[lev], 0, 0, 1, 0);
    }
    
    double totCharge = totalCharge(rhs, finest_level, geom);
    msg << "Total Charge: " << totCharge << " C" << endl
        << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << vol << " m^3" << endl;
    
    // eps in C / (V * m)
    double constant = -1.0 / Physics::epsilon_0;  // in [V m / C]
    for (int i = 0; i <=finest_level; ++i)
        rhs[i]->mult(constant, 0, 1);       // in [V m]
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;

    // solve                                                                                                                                                                                                     
    Solver sol;
    sol.solve_for_accel(rhs,
                        phi,
                        efield,
                        geom,
                        base_level,
                        finest_level,
                        offset);
    
    // for plotting unnormalize
    for (int i = 0; i <=finest_level; ++i)
        rhs[i]->mult(1.0 / constant, 0, 1);       // in [V m]
}

void doAMReX(const Vektor<size_t, 3>& nr, size_t nParticles,
              int nLevels, size_t maxBoxSize, Inform& msg)
{
    static IpplTimings::TimerRef distTimer = IpplTimings::getTimer("init-distribution");
    static IpplTimings::TimerRef regridTimer = IpplTimings::getTimer("regrid");
    static IpplTimings::TimerRef redistTimer = IpplTimings::getTimer("particle-redistr");
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    std::array<double, AMREX_SPACEDIM> lower = {{-0.0525, -0.0525, -0.0525}}; // m
    std::array<double, AMREX_SPACEDIM> upper = {{ 0.0525,  0.0525,  0.0525}}; // m
    
    RealBox domain;
    
    // in helper_functions.h
    init(domain, nr, lower, upper);
    
    
    /*
     * create an Amr object
     */
    ParmParse pp("amr");
    pp.add("max_grid_size", int(maxBoxSize));
    
    Vector<int> error_buf(nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    Vector<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr[i];
    
    AmrOpal myAmrOpal(&domain, nLevels - 1, nCells, 0 /* cartesian */);
    
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const Vector<BoxArray>& ba = myAmrOpal.boxArray();
    const Vector<DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    const Vector<Geometry>& geom = myAmrOpal.Geom();
    
    Vector<int> rr(nLevels);
    for (int i = 0; i < nLevels; ++i)
        rr[i] = 2;
    
    
    //create a new layout using ParticleAmrLayout class
    amrplayout_t* PL = new amrplayout_t(geom, dmap, ba, rr);
    
    //create a particle bunch
    PartBunchAmr<amrplayout_t>* bunch = new PartBunchAmr<amrplayout_t>();
    bunch->initialize(PL);
    bunch->initializeAmr();
    
    
    // initialize a particle distribution
    unsigned long int nloc = nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    IpplTimings::startTimer(distTimer);
    double sig = 0.01;  // m
    dist.gaussian(0.0, sig, nloc, ParallelDescriptor::MyProc());
    
    // copy particles to the PartBunchBase object.
    dist.injectBeam(*bunch);
    IpplTimings::stopTimer(distTimer);
    
    
    
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
        bunch->qm[i] = Physics::q_e;  // in [C]
    
    // redistribute on single-level
    IpplTimings::startTimer(redistTimer);
    bunch->update();
    IpplTimings::stopTimer(redistTimer);
    
    Vector_t rmin, rmax;
    bounds(bunch->R, rmin, rmax);
    msg << "rmin = " << rmin << endl << "rmax = " << rmax << endl;
    
    msg << "Bunch radius: " << sig << " m" << endl
        << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << nParticles * bunch->qm[0] << " C" << endl;
    
    
    myAmrOpal.setBunch(bunch);
        
    const Vector<Geometry>& geoms = myAmrOpal.Geom();
    for (int i = 0; i < nLevels; ++i)
        msg << "#Cells per dim of level " << i << " for bunch : "
            << 2.0 * sig / *(geoms[i].CellSize()) << endl;
    
    // ========================================================================
    // 2. tagging (i.e. create BoxArray's, DistributionMapping's of all
    //    other levels)
    // ========================================================================
    
    /*
     * do tagging
     */
//     dynamic_cast<AmrPartBunch*>(bunch)->Define (myAmrOpal.Geom(),
//                                                 myAmrOpal.DistributionMap(),
//                                                 myAmrOpal.boxArray(),
//                                                 rr);
//     
//     
    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    IpplTimings::startTimer(regridTimer);
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    IpplTimings::stopTimer(regridTimer);
    
    container_t rhs, phi, efield;
    
    std::string plotsolve = amrex::Concatenate("plt", 0, 4);
    
    doSolve(myAmrOpal, bunch, rhs, phi, efield, geoms, rr, nLevels, msg);
    
    msg << "Total field energy: " << totalFieldEnergy(efield, rr) << endl;
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(0) << endl;
    }
    
    
//     writePlotFile(plotsolve, rhs, phi, efield, rr, geoms, 0);
    
//     bunch->>python_format(0);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z] [#particles] "
         << "[#levels] [max. box size] [out: timing file name (optiona)]";
    
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
        
    amrex::Initialize(argc,argv, false);
    size_t nLevels = std::atoi(argv[5]) + 1; // i.e. nLevels = 0 --> only single level
    size_t maxBoxSize = std::atoi(argv[6]);
    doAMReX(nr, nParticles, nLevels, maxBoxSize, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    std::stringstream timefile;
    timefile << std::string(argv[0]) << "-timing-cores-"
             << std::setfill('0') << std::setw(6) << Ippl::getNodes()
             << "-threads-1.dat";
    
    if ( argc == 8 ) {
        timefile.str(std::string());
        timefile << std::string(argv[7]);
    }
    
    IpplTimings::print(timefile.str());
    
    return 0;
}
