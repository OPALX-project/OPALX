/*!
 * @file testPerformance.cpp
 * @author Matthias Frey
 * @date February 2017
 * 
 * Domain:  [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5]\n
 * BC:      Perodic boundary conditions\n
 * Charge/particle: elementary charge\n
 * Gaussian particle distribution N(0.0, 0.01)
 * 
 * @details Do some timing measurements on the different parts of a
 *          tracking code, i.e. computation of force, computation of Poisson
 *          equation, etc.\n
 *          In every time step each particle is randomly perturbed by a small \f$\delta x\f$.
 * 
 * Call:\n
 *  mpirun -np [#cores] testPerformance [#gridpoints x] [#gridpoints y] [#gridpoints z]
 *                                     [#particles] [#levels] [max. box size] [#steps]
 * 
 * @brief Perturbes particles randomly in space for several time steps.
 */

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>
#include <iomanip>


#include <ParmParse.H>


#include "../Distribution.h"
#include "../Solver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

// #include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

typedef Vektor<double, BL_SPACEDIM> Vector_t;

void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             int nLevels,
             Inform& msg,
             IpplTimings::TimerRef& assignTimer)
{
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

    IpplTimings::startTimer(assignTimer);
    bunch->AssignDensity(bunch->qm, false, rhs, base_level, finest_level);
    IpplTimings::stopTimer(assignTimer);
    
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
    sol.solve_for_accel(rhs,
                        phi,
                        grad_phi,
                        geom,
                        base_level,
                        finest_level,
                        offset,
                        false);
    
    // for plotting unnormalize
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(1.0 / constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(1.0 / constant, 0, 1);
#endif
    }
}

void doBoxLib(const Vektor<size_t, 3>& nr, size_t nParticles,
              int nLevels, size_t maxBoxSize, int nSteps, Inform& msg)
{
    static IpplTimings::TimerRef regridTimer = IpplTimings::getTimer("tracking-regrid");
    static IpplTimings::TimerRef solveTimer = IpplTimings::getTimer("tracking-solve");
    static IpplTimings::TimerRef stepTimer = IpplTimings::getTimer("tracking-step");
    static IpplTimings::TimerRef totalTimer = IpplTimings::getTimer("tracking-total");
    static IpplTimings::TimerRef statisticsTimer = IpplTimings::getTimer("tracking-statistics");
    static IpplTimings::TimerRef assignTimer = IpplTimings::getTimer("assign-charge");
    static IpplTimings::TimerRef gravityTimer = IpplTimings::getTimer("gravity-eval");
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    std::array<double, BL_SPACEDIM> lower = {{-0.5, -0.5, -0.5}}; // m
    std::array<double, BL_SPACEDIM> upper = {{ 0.5,  0.5,  0.5}}; // m
    
    RealBox domain;
    
    // in helper_functions.h
    init(domain, nr, lower, upper);
    
    
    /*
     * create an Amr object
     */
    ParmParse pp("amr");
    pp.add("max_grid_size", int(maxBoxSize));
    
    Array<int> error_buf(nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
//     ParmParse pgeom("geometry");
//     Array<int> is_per = { 1, 1, 1};
//     pgeom.addarr("is_periodic", is_per);
    
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr[i];
    
    AmrOpal myAmrOpal(&domain, nLevels - 1, nCells, 0 /* cartesian */);
    
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const Array<BoxArray>& ba = myAmrOpal.boxArray();
    const Array<DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    const Array<Geometry>& geom = myAmrOpal.Geom();
    
    Array<int> rr(nLevels);
    for (int i = 0; i < nLevels; ++i)
        rr[i] = 2;
    
    
    amrplayout_t* playout = new amrplayout_t(geom, dmap, ba, rr);
    
    std::unique_ptr<amrbunch_t> bunch( new amrbunch_t() );
    bunch->initialize(playout);
    bunch->initializeAmr(); // add attributes: level, grid
    
    bunch->setAllowParticlesNearBoundary(true);
    
    // initialize a particle distribution
    unsigned long int nloc = nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    dist.gaussian(0.0, 0.07, nloc, ParallelDescriptor::MyProc());
    
    // copy particles to the PartBunchBase object.
    dist.injectBeam( *(bunch.get()) );
    
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
        bunch->qm[i] = 1.0e-19;  // in [C]
    }
    
    // redistribute on single-level
    bunch->update();
    
    msg << "Single-level statistics" << endl;
    bunch->gatherStatistics();
    
    msg << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << nParticles * bunch->qm[0] << " C" << endl;
    
    
    myAmrOpal.setBunch(bunch.get());
    
    const Array<Geometry>& geoms = myAmrOpal.Geom();
    
    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    msg << "Multi-level statistics" << endl;
    bunch->gatherStatistics();
    
    container_t rhs;
    container_t phi;
    container_t grad_phi;
    
//     std::string plotsolve = BoxLib::Concatenate("plt", 0, 4);
    
    /* units:
     * 
     * [P] = \gamma\beta
     * [R] = m
     * 
     */
    std::mt19937_64 eng(0);
    std::uniform_real_distribution<> perturbation(-0.01, 0.01);
    
    
    std::string statistics = "particle-statistics-ncores-" + std::to_string(Ippl::getNodes()) + ".dat";
    IpplTimings::startTimer(totalTimer);
    for (int t = 0; t < nSteps; ++t) {
        
        //        bunch->python_format(t);
        
        IpplTimings::startTimer(solveTimer);
        doSolve(myAmrOpal, bunch.get(), rhs, phi, grad_phi, geoms, rr, nLevels, msg, assignTimer);
        IpplTimings::stopTimer(solveTimer);
        
        IpplTimings::startTimer(gravityTimer);
        bunch->GetGravity(bunch->E, grad_phi);
        IpplTimings::stopTimer(gravityTimer);
        
        IpplTimings::startTimer(stepTimer);
        for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
            bunch->R[i] += Vector_t(perturbation(eng),
                                    perturbation(eng),
                                    perturbation(eng)
                                   );
        }
        IpplTimings::stopTimer(stepTimer);
        
        IpplTimings::startTimer(regridTimer);
        if ( myAmrOpal.maxLevel() > 0 ) {
        for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
            myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
        } else
            bunch->update();
        IpplTimings::stopTimer(regridTimer);
        
        IpplTimings::startTimer(statisticsTimer);
        bunch->dumpStatistics(statistics);
        IpplTimings::stopTimer(statisticsTimer);
    }
    IpplTimings::stopTimer(totalTimer);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z] [#particles] "
         << "[#levels] [max. box size] [#steps]";
    
    if ( argc < 8 ) {
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
    int nSteps = std::atoi(argv[7]);
    doBoxLib(nr, nParticles, nLevels, maxBoxSize, nSteps, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    std::stringstream timefile;
    timefile << std::string(argv[0]) << "-timing-cores-"
             << std::setfill('0') << std::setw(6) << Ippl::getNodes()
             << "-threads-1.dat";
    
    IpplTimings::print(timefile.str());
    
    return 0;
}
