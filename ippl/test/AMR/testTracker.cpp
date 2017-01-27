/*!
 * @file testTracker.cpp
 * @author Matthias Frey
 * @date November 2016
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
 *  mpirun -np [#cores] testTracker [#gridpoints x] [#gridpoints y] [#gridpoints z]
 *                                  [#particles] [#levels] [max. box size]
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

Vector_t kick(const Vector_t& R,
              Vector_t P,
              const Vector_t& Ef,
              const Vector_t& Bf,
              const double& dt,
              const double& mass,
              const double& charge)
{
    double const gamma = std::sqrt(1.0 + dot(P, P));
    const Vector_t t = 0.5 * dt * charge * Physics::c * Physics::c / ( gamma * mass) * Bf;
    const Vector_t w = P + cross(P, t);
    const Vector_t s = 2.0 / (1.0 + dot(t, t)) * t;
    P += cross(w, s);
    
    return P + 0.5 * dt * charge * Physics::c / mass * Ef;
}

Vector_t push(Vector_t R, const Vector_t& P, const double& dt) {
    return R + 0.5 * dt * P / std::sqrt(1.0 + dot(P, P));
}

void doSolve(AmrOpal& myAmrOpal, PartBunchBase* bunch,
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
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rhs, base_level, 1, finest_level);
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
              int nLevels, size_t maxBoxSize, Inform& msg)
{
    static IpplTimings::TimerRef regridTimer = IpplTimings::getTimer("tracking-regrid");
    static IpplTimings::TimerRef solveTimer = IpplTimings::getTimer("tracking-solve");
    static IpplTimings::TimerRef stepTimer = IpplTimings::getTimer("tracking-step");
    static IpplTimings::TimerRef totalTimer = IpplTimings::getTimer("tracking-total");
    static IpplTimings::TimerRef statisticsTimer = IpplTimings::getTimer("tracking-statistics");
    static IpplTimings::TimerRef assignTimer = IpplTimings::getTimer("assign-charge");
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    std::array<double, BL_SPACEDIM> lower = {{-0.15, -0.15, -0.15}}; // m
    std::array<double, BL_SPACEDIM> upper = {{ 0.15,  0.15,  0.15}}; // m
    
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
    
    PartBunchBase* bunch = new AmrPartBunch(geom, dmap, ba, rr);
    
    
    // initialize a particle distribution
    unsigned long int nloc = nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    double sig = 0.02;  // m
    dist.gaussian(0.0, sig, nloc, ParallelDescriptor::MyProc());
    
    // copy particles to the PartBunchBase object.
    dist.injectBeam(*bunch);
    
    
    
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
        bunch->setQM(1.0e-8, i);  // in [C]
    }
    
    // redistribute on single-level
    bunch->myUpdate();
    
    msg << "Single-level statistics" << endl;
    bunch->gatherStatistics();
    
    msg << "Bunch radius: " << sig << " m" << endl
        << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->getQM(0) << " C" << endl
        << "Total charge: " << nParticles * bunch->getQM(0) << " C" << endl;
    
    
    myAmrOpal.setBunch(dynamic_cast<AmrPartBunch*>(bunch));
        
    const Array<Geometry>& geoms = myAmrOpal.Geom();
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
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    msg << "Multi-level statistics" << endl;
    bunch->gatherStatistics();
    
    container_t rhs;
    container_t phi;
    container_t grad_phi;
    
    std::string plotsolve = BoxLib::Concatenate("plt", 0, 4);
    
    double mass = 0.983272; // mass
    double dt = 0.005;
    
    /* units:
     * 
     * [P] = \gamma\beta
     * [R] = m
     * 
     */
    
    std::string statistics = "particle-statistics-ncores-" + std::to_string(Ippl::getNodes()) + ".dat";
    IpplTimings::startTimer(totalTimer);
    for (int t = 0; t < 100; ++t) {
        
        IpplTimings::startTimer(solveTimer);
        doSolve(myAmrOpal, bunch, rhs, phi, grad_phi, geoms, rr, nLevels, msg, assignTimer);
        IpplTimings::stopTimer(solveTimer);
        
        IpplTimings::startTimer(stepTimer);
        for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
            int level = dynamic_cast<AmrPartBunch*>(bunch)->getLevel(i);
            
            
            Vector_t R = push(bunch->getR(i),
                              bunch->getP(i),
                              dt);
            
            bunch->setR(R, i);
            
            Vector_t externalE = dynamic_cast<AmrPartBunch*>(bunch)->interpolate(i, grad_phi[level]);
            Vector_t externalB = Vector_t(0.0, 0.0, 0.0);
            
            double gamma = std::sqrt(1.0 + dot(bunch->getP(i), bunch->getP(i)));
            double beta  = std::sqrt( 1.0 - 1.0 / ( gamma * gamma ) );
            
            externalB[1] =  gamma * beta * externalE[2];
            externalB[2] = -gamma * beta * externalE[1];
            
            externalE[1] *= gamma;
            externalE[2] *= gamma;
            
            
            
            Vector_t P = kick(bunch->getR(i),
                              bunch->getP(i),
                              externalE,
                              externalB,
                              dt,
                              mass * 1.0e9,
                              bunch->getQM(i));
            
            bunch->setP(P, i);
            
            R = push(bunch->getR(i),
                     bunch->getP(i),
                     dt);
            
            bunch->setR(R, i);
        }
        IpplTimings::stopTimer(stepTimer);
        
        IpplTimings::startTimer(regridTimer);
        for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
            myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
        IpplTimings::stopTimer(regridTimer);
        
        IpplTimings::startTimer(statisticsTimer);
        bunch->dumpStatistics(statistics);
        IpplTimings::stopTimer(statisticsTimer);
        
//         writePlotFile(plotsolve, rhs, phi, grad_phi, rr, geoms, t);
//         dynamic_cast<AmrPartBunch*>(bunch)->python_format(t);
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
    
    std::stringstream timefile;
    timefile << std::string(argv[0]) << "-timing-cores-"
             << std::setfill('0') << std::setw(6) << Ippl::getNodes()
             << "-threads-1.dat";
    
    IpplTimings::print(timefile.str());
    
    return 0;
}
