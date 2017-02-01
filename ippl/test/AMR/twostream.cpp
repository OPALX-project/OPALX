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

void doSolve(AmrOpal& myAmrOpal, PartBunchBase* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             int nLevels)
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

    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rhs, base_level, 1, finest_level);
    
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

void doTwoStream(Vektor<std::size_t, 3> nr,
                 std::size_t nLevels,
                 std::size_t maxBoxSize,
                 double alpha,
                 double dt,
                 double epsilon,
                 std::size_t nIter,
                 Inform& msg)
{
    std::array<double, BL_SPACEDIM> lower = {{0.0, 0.0, 0.0}}; // m
    std::array<double, BL_SPACEDIM> upper = {{4.0 * Physics::pi,
                                              4.0 * Physics::pi,
                                              4.0 * Physics::pi}
    }; // m
    
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
    
    ParmParse pgeom("geometry");
    Array<int> is_per = { 1, 1, 1};
    pgeom.addarr("is_periodic", is_per);
    
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
    for (std::size_t i = 0; i < nLevels; ++i)
        rr[i] = 2;
    
    PartBunchBase* bunch = new AmrPartBunch(geom, dmap, ba, rr);
    
    
    // initialize a particle distribution
    Distribution dist;
    dist.twostream(Vektor<double, 3>(lower[0], lower[1], lower[2]),
                   Vektor<double, 3>(upper[0], upper[1], upper[2]),
                   Vektor<std::size_t, 3>(4, 4, 8/*32*/),
                   Vektor<std::size_t, 3>(8, 8, 16/*128*/),
                   Vektor<double, 3>(3/*6.0*/, 3/*6.0*/, 3/*6.0*/), alpha);
    
    // copy particles to the PartBunchBase object.
    dist.injectBeam(*bunch);
    
    // redistribute on single-level
    bunch->myUpdate();
    
//     dynamic_cast<AmrPartBunch*>(bunch)->SetAllowParticlesNearBoundary(true);
    
//     writePlotFile(plotsolve, rhs, phi, grad_phi, rr, geoms, 0);
    
    myAmrOpal.setBunch(dynamic_cast<AmrPartBunch*>(bunch));
    
    const Array<Geometry>& geoms = myAmrOpal.Geom();
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    msg << "Multi-level statistics" << endl;
    bunch->gatherStatistics();
    
    for (std::size_t i = 0; i < nIter; ++i) {
        
        dynamic_cast<AmrPartBunch*>(bunch)->python_format(i);
        
        for (std::size_t j = 0; j < bunch->getLocalNum(); ++j) {
            int level = dynamic_cast<AmrPartBunch*>(bunch)->getLevel(j);
            bunch->setR( bunch->getR(j) + dt * bunch->getP(j), j );
        }
            
        for (int j = 0; j <= myAmrOpal.finestLevel() && j < myAmrOpal.maxLevel(); ++j)
            myAmrOpal.regrid(j /*lbase*/, 0.0 /*time*/);
        
        container_t rhs;
        container_t phi;
        container_t grad_phi;
        doSolve(myAmrOpal, bunch, rhs, phi, grad_phi, geoms, rr, nLevels);
        
        for (std::size_t j = 0; j < bunch->getLocalNum(); ++j) {
            int level = dynamic_cast<AmrPartBunch*>(bunch)->getLevel(j);
            Vector_t Ef = dynamic_cast<AmrPartBunch*>(bunch)->interpolate(j, grad_phi[level]);
            
            Ef *= Physics::epsilon_0; /* not used by Ulmer */
            bunch->setP( bunch->getP(j) + dt * bunch->getQM(j) / bunch->getMass(j) * Ef, j);
        }
        
    }
    
    
    
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    BoxLib::Initialize(argc,argv, false);
    
    Inform msg("TwoStream");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z] [#levels]"
         << " [max. box size] [alpha] [timstep] [epsilon] [#iterations]"
         << " [out: timing file name (optiona)]";
    
    if ( argc < 9 ) {
        msg << call.str() << endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<std::size_t, 3> nr(std::atoi(argv[1]),
                              std::atoi(argv[2]),
                              std::atoi(argv[3]));
    
    
    std::size_t nLevels    = std::atoi( argv[4] ) + 1; // i.e. nLevels = 0 --> only single level
    std::size_t maxBoxSize = std::atoi( argv[5] );
    double alpha           = std::atof( argv[6] );
    double dt              = std::atof( argv[7] );
    double epsilon         = std::atof( argv[8] );
    std::size_t nIter      = std::atof( argv[9] );
    
    msg << "Particle test running with" << endl
        << "- #level     = " << nLevels << endl
        << "- grid       = " << nr << endl
        << "- max. size  = " << maxBoxSize << endl
        << "- amplitude  = " << alpha << endl
        << "- epsilon    = " << epsilon << endl
        << "- time step  = " << dt << endl
        << "- #steps     = " << nIter << endl;
    
    doTwoStream(nr, nLevels, maxBoxSize,
                alpha, dt, epsilon, nIter, msg);
    
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