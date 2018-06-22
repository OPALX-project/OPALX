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
 *                                     [#particles] [#levels] [max. box size] [#steps] [reuse]
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


#include <AMReX_ParmParse.H>

#include <AMReX_Interpolater.H>


#include "../Distribution.h"
#include "../Solver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

// #include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

using namespace amrex;

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

typedef Vektor<double, AMREX_SPACEDIM> Vector_t;

void precondition(AmrOpal& myAmrOpal,
                  container_t& rhs,
                  container_t& phi,
                  container_t& grad_phi,
                  int nLevels,
                  const Vector<Geometry>& geom)
{
    for (int lev = 0; lev < nLevels; ++lev) {
        initGridData(rhs, grad_phi, myAmrOpal.boxArray()[lev],
                     myAmrOpal.DistributionMap(lev), lev);
    }
    
    for (int i = myAmrOpal.finestLevel()-1; i >= 0; --i) {
        MultiFab tmp(phi[i]->boxArray(), phi[i]->DistributionMap(), 1, 0);
        tmp.setVal(0.0);
        amrex::average_down(*(phi[i+1].get()), tmp, 0, 1, myAmrOpal.refRatio(i));
        MultiFab::Add(*(phi[i].get()), tmp, 0, 0, 1, 0);
    }
    
    for (int lev = 1; lev < nLevels; ++lev) {
        phi[lev].reset(new MultiFab(myAmrOpal.boxArray()[lev],myAmrOpal.DistributionMap(lev),1,1));
        phi[lev]->setVal(0.0);
    }
    
    
    
        
    PCInterp mapper;
    Vector<BCRec> bc(1);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        bc[0].setLo(i, INT_DIR);
        bc[0].setHi(i, INT_DIR);
    }
    
    IntVect fine_ratio(2, 2, 2);
    for (int lev = 0; lev < nLevels - 1; ++lev) {
        /*
         * fmfi: MultiFab iterator for fine grids
         * cmfi: MultiFab iterator for coarse grids
         */
        for (MFIter fmfi(*(phi[lev + 1].get()), false);
             fmfi.isValid(); ++fmfi) {
            
            const Box& bx = fmfi.validbox();
            FArrayBox& fab = (*(phi[lev + 1].get()))[fmfi];
            
            
            FArrayBox finefab(bx, 1);
            FArrayBox crsefab(mapper.CoarseBox(finefab.box(), fine_ratio), 1);
            
            for (MFIter cmfi(*(phi[lev].get()), false);cmfi.isValid(); ++cmfi) {
                crsefab.copy((*(phi[lev].get()))[cmfi]);
            }
            
            mapper.interp(crsefab,
                          0, // comp
                          finefab,
                          0, // comp
                          1, // ncomp
                          finefab.box(),
                          fine_ratio,
                          geom[lev],
                          geom[lev + 1],
                          bc,
                          0,
                          0);
            
            fab.copy(finefab);
        }
    }
}

void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Vector<Geometry>& geom,
             const Vector<int>& rr,
             int nLevels,
             bool reuse,
             Inform& msg,
             IpplTimings::TimerRef& assignTimer,
             IpplTimings::TimerRef& actualSolveTimer)
{
    // =======================================================================                                                                                                                                   
    // 4. prepare for multi-level solve                                                                                                                                                                          
    // =======================================================================
    
    if ( reuse ) {
        precondition(myAmrOpal, rhs, phi, grad_phi, nLevels, geom);
        
    } else {
        for (int lev = 0; lev < nLevels; ++lev) {
            initGridData(rhs, phi, grad_phi, myAmrOpal.boxArray()[lev],
                         myAmrOpal.DistributionMap(lev), lev);
        }
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
        rhs[i]->mult(constant, 0, 1);       // in [V m]
    }
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;

    // solve                                                                                                                                                                                                     
    Solver sol;
    IpplTimings::startTimer(actualSolveTimer);
    sol.solve_for_accel(rhs,
                        phi,
                        grad_phi,
                        geom,
                        base_level,
                        finest_level,
                        offset,
                        false);
    IpplTimings::stopTimer(actualSolveTimer);
    
    // for plotting unnormalize
    for (int i = 0; i <=finest_level; ++i) {
        rhs[i]->mult(1.0 / constant, 0, 1);       // in [V m]
    }
}

void doBoxLib(const Vektor<size_t, 3>& nr, size_t nParticles,
              int nLevels, size_t maxBoxSize, int nSteps, bool reuse, Inform& msg)
{
    static IpplTimings::TimerRef regridTimer = IpplTimings::getTimer("tracking-regrid");
    static IpplTimings::TimerRef solveTimer = IpplTimings::getTimer("tracking-solve");
    static IpplTimings::TimerRef stepTimer = IpplTimings::getTimer("tracking-step");
    static IpplTimings::TimerRef totalTimer = IpplTimings::getTimer("tracking-total");
    static IpplTimings::TimerRef statisticsTimer = IpplTimings::getTimer("tracking-statistics");
    static IpplTimings::TimerRef assignTimer = IpplTimings::getTimer("assign-charge");
    static IpplTimings::TimerRef actualSolveTimer = IpplTimings::getTimer("actual-solve");
    static IpplTimings::TimerRef gravityTimer = IpplTimings::getTimer("gravity-eval");
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    std::array<double, AMREX_SPACEDIM> lower = {{-0.5, -0.5, -0.5}}; // m
    std::array<double, AMREX_SPACEDIM> upper = {{ 0.5,  0.5,  0.5}}; // m
    
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
    
//     ParmParse pgeom("geometry");
//     Vector<int> is_per = { 1, 1, 1};
//     pgeom.addarr("is_periodic", is_per);
    
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
    
    const Vector<Geometry>& geoms = myAmrOpal.Geom();
    
    myAmrOpal.setTagging(AmrOpal::kCenteredRegion);
    
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
    rhs.resize(nLevels);
    phi.resize(nLevels);
    grad_phi.resize(nLevels);
    
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
        doSolve(myAmrOpal, bunch.get(), rhs, phi, grad_phi, geoms,
                rr, nLevels, (reuse && bool(t)), msg, assignTimer, actualSolveTimer);
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
        
        bunch->gatherStatistics();
        
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
         << "[#levels] [max. box size] [#steps] [reuse solution of previous step: true / false]";
    
    if ( argc < 9 ) {
        msg << call.str() << endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<size_t, 3> nr(std::atoi(argv[1]),
                         std::atoi(argv[2]),
                         std::atoi(argv[3]));
    
    
    size_t nParticles = std::atoi(argv[4]);
    
    
        
    amrex::Initialize(argc,argv, false);
    size_t nLevels = std::atoi(argv[5]) + 1; // i.e. nLevels = 0 --> only single level
    size_t maxBoxSize = std::atoi(argv[6]);
    int nSteps = std::atoi(argv[7]);
    std::string reuse_str = argv[8];
    bool reuse = false;
    
    if ( reuse_str.find('1') != std::string::npos ||
         reuse_str.find("true") != std::string::npos )
        reuse = true;
        
    
    
    msg << "Particle test running with" << endl
        << "- #particles      = " << nParticles << endl
        << "- grid            = " << nr << endl
        << "- max grid        = " << maxBoxSize << endl
        << "- #steps          = " << nSteps << endl
        << "- #level          = " << nLevels << endl
        << "- reuse prev. sol = " << ((reuse) ? "True" : "False") << endl;
    
    
    doBoxLib(nr, nParticles, nLevels, maxBoxSize, nSteps, reuse, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    std::stringstream timefile;
    timefile << std::string(argv[0]) << "-timing-cores-"
             << std::setfill('0') << std::setw(6) << Ippl::getNodes()
             << "-threads-1.dat";
    
    IpplTimings::print(timefile.str());
    
    return 0;
}
