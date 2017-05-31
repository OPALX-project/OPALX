/*!
 * @file testDomainTransformSolve.cpp
 * @author Matthias Frey
 * @date May 2017
 * @details Compute \f$\Delta\phi = -\rho\f$ where the charge
 * density is a Gaussian distribution. Write plot files
 * that can be visualized by yt (python visualize.py)
 * or by AmrVis of the CCSE group at LBNL.
 * The calculation is done on a computation domain [-1, 1]^3.
 * The particles are scaled into the domain and then scaled back. The
 * potential and electric field are scaled appropriately.
 * 
 * Domain:  [-1.0, 1.0] x [-1.0, 1.0] x [-1.0, 1.0]\n
 * BC:      Dirichlet boundary conditions\n
 * Charge/particle: elementary charge\n
 * Gaussian particle distribution N(0.0, 0.01)
 * 
 * Call:\n
 *  mpirun -np [#cores] testNewTracker [#gridpoints x] [#gridpoints y] [#gridpoints z]
 *                                     [#particles] [#levels] [max. box size] [#steps]
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

double domainMapping(amrbase_t& PData,
                     const Vector_t& lold,
                     const Vector_t& uold,
                     const Vector_t& lnew,
                     const Vector_t& unew)
{
    // [lold, uold] --> [lnew, unew]
    Vector_t invdiff = 1.0 / ( uold - lold );
    Vector_t slope = ( unew - lnew ) * invdiff;
    Vector_t intercept = ( uold * lnew - lold * unew ) * invdiff;
    
    for (unsigned int i = 0; i < PData.getLocalNum(); ++i)
        PData.R[i] = slope * PData.R[i] + intercept;
    
    Vector_t dnew = unew - lnew;
    Vector_t dold = uold - lold;
    
    return  std::sqrt( dot(dnew, dnew) / dot(dold, dold) ); //( unew - lnew ) / ( uold - lold );
}

void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& efield,
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
    efield.resize(nLevels);
    
    for (int lev = 0; lev < nLevels; ++lev) {
        initGridData(rhs, phi, efield,
                     myAmrOpal.boxArray()[lev], myAmrOpal.DistributionMap(lev), lev);
    }

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();
    
   container_t partMF(nLevels);
   for (int lev = 0; lev < nLevels; lev++) {
        const amrex::BoxArray& ba = myAmrOpal.boxArray()[lev];
        const amrex::DistributionMapping& dmap = myAmrOpal.DistributionMap(lev);
        partMF[lev].reset(new amrex::MultiFab(ba, dmap, 1, 2));
        partMF[lev]->setVal(0.0, 2);
   }
    
    IpplTimings::startTimer(assignTimer);
    bunch->AssignDensityFort(bunch->qm, partMF, base_level, 1, finest_level);
    IpplTimings::stopTimer(assignTimer);
    
    for (int lev = 0; lev < nLevels; ++lev) {
        amrex::MultiFab::Copy(*rhs[lev], *partMF[lev], 0, 0, 1, 0);
    }
    
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
                        efield,
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

void doAMReX(const Vektor<size_t, 3>& nr, size_t nParticles,
              int nLevels, size_t maxBoxSize, int nSteps, Inform& msg)
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
    
    std::array<double, BL_SPACEDIM> lower = {{-1.0, -1.0, -1.0}}; // m
    std::array<double, BL_SPACEDIM> upper = {{ 1.0,  1.0,  1.0}}; // m
    
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
    Array<int> is_per = { 0, 0, 0};
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
    dist.gaussian(0.0, 0.01, nloc, ParallelDescriptor::MyProc());
    
    // copy particles to the PartBunchBase object.
    dist.injectBeam( *bunch );
    
    
    
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
        bunch->qm[i] = Physics::q_e;  // in [C]
    
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
    container_t efield;
    
    std::string plotsolve = amrex::Concatenate("plt", 0, 4);
    
    bunch->python_format(0);
    
    
    // map particles
    Vector_t rmin, rmax;
    bounds(bunch->R, rmin, rmax);
    msg << "rmin = " << rmin << endl << "rmax = " << rmax << endl;
    
    Vector_t low(-1.0, -1.0, -1.0);
    Vector_t up = - low;
    
    msg << endl << "Transformed positions" << endl << endl;
    
    double factor = domainMapping(*bunch, 1.05 * rmin, 1.05 * rmax, low , up);
    
    bunch->update();
    
    Vector_t tmp_min, tmp_max;
    bounds(bunch->R, tmp_min, tmp_max);
    msg << "rmin = " << tmp_min << endl << "rmax = " << tmp_max << endl;
    
    doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, geoms, rr, nLevels, msg, assignTimer);
    
    msg << endl << "Back to normal positions" << endl << endl;
    
    domainMapping(*bunch, low , up, 1.05 * rmin, 1.05 * rmax);
    
    bunch->update();
    
    bounds(bunch->R, tmp_min, tmp_max);
    
    msg << "rmin = " << tmp_min << endl << "rmax = " << tmp_max << endl;
    
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        if ( efield[i]->contains_nan(false) )
            msg << "Efield: Nan" << endl;
    }
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        
        msg << i << ": " << phi[i]->nGrow() << " "
                  << phi[i]->min(0, 1) << " " << phi[i]->max(0, 1) << endl;
        
        if ( phi[i]->contains_nan(false) )
            msg << "Pot: Nan" << endl;
    }
    
    // scale solution
    double factor2 = factor * factor;
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        phi[i]->mult(factor, 0, 1);
        
        efield[i]->mult(factor2, 0, 3, 1);
    }
    
    msg << "Total field energy: " << totalFieldEnergy(efield, rr) << endl;
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(0) << endl;
    }
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
        
    amrex::Initialize(argc,argv, true);
    size_t nLevels = std::atoi(argv[5]) + 1; // i.e. nLevels = 0 --> only single level
    size_t maxBoxSize = std::atoi(argv[6]);
    int nSteps = std::atoi(argv[7]);
    doAMReX(nr, nParticles, nLevels, maxBoxSize, nSteps, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    std::stringstream timefile;
    timefile << std::string(argv[0]) << "-timing-cores-"
             << std::setfill('0') << std::setw(6) << Ippl::getNodes()
             << "-threads-1.dat";
    
    IpplTimings::print(timefile.str());
    
    return 0;
}
