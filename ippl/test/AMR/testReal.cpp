// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 *
 ***************************************************************************/

/***************************************************************************

This test program sets up a simple sine-wave electric field in 3D,
  creates a population of particles with random q/m values (charge-to-mass
  ratio) and velocities, and then tracks their motions in the static
  electric field using nearest-grid-point interpolation.

Usage:
 mpirun -np 4 testReal [#gridpoints x] [#gridpoints y] [#gridpoints z]
                       [#levels] [max. box size] [h5 file] [start:by:end]

***************************************************************************/

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>
#include <tuple>

#include <ParmParse.H>

#include "PartBunch.h"
#include "AmrPartBunch.h"
#include "helper_functions.h"

#include "Distribution.h"
#include "Solver.h"
#include "AmrOpal.h"

#include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

void doSolve(AmrOpal& myAmrOpal, PartBunchBase* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             int nLevels,
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
    
    double totCharge = totalCharge(rhs, finest_level, geom);
    
    msg << "Total Charge: " << totCharge << " C" << endl
        << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << vol << " m^3" << endl;
    
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
    sol.solve_for_accel(rhs,
                        phi,
                        grad_phi,
                        geom,
                        base_level,
                        finest_level,
                        offset);
    
    // for plotting unnormalize
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(1.0 / constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(1.0 / constant, 0, 1);
#endif
    }
    
    IpplTimings::stopTimer(solvTimer);
}

void doBoxLib(const Vektor<size_t, 3>& nr,
              int nLevels, size_t maxBoxSize,
              std::string h5file,
	      const std::tuple<int, int, int>& h5steps,
	      Inform& msg)
{
    static IpplTimings::TimerRef distTimer = IpplTimings::getTimer("dist");
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    std::array<double, BL_SPACEDIM> lower = {{-0.15, -0.15, -0.15}}; // m
    std::array<double, BL_SPACEDIM> upper = {{0.45, 0.15, 0.15}}; // m
    
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
    Distribution dist;
    IpplTimings::startTimer(distTimer);

    int begin = 0, end = 0, by = 0;
    unsigned int nBunches = 0;

    std::vector<double> xshift = {0.037,
				  0.05814,
				  0.0852,
				  0.11325,
				  0.13933,
				  0.16177,
				  0.18013,
				  0.19567,
				  0.21167,
				  0.2313,
				  0.25567,
				  0.28326,
				  0.31034,
				  0.33308,
				  0.35041};
    
    std::tie(begin, end, by) = h5steps;
    for (int i = std::min(begin, end); i <= std::max(begin, end); i += by) {
	msg << "Reading step " << i << endl;
	dist.readH5(h5file, i /* step */);
	// copy particles to the PartBunchBase object.
	dist.injectBeam(*bunch, false, {{xshift[nBunches], 0.0, 0.0}});
	msg << "Injected step " << i << endl;

	if ( ++nBunches == xshift.size() )
	    break;
    }

    msg << "#Bunches: " << nBunches << endl;

    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
	bunch->setQM(2.1717e-16, i);
    }
    
    // redistribute on single-level
    bunch->myUpdate();
    
    bunch->gatherStatistics();
    
    int nParticles = bunch->getTotalNum();
    msg << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->getQM(0) << " C" << endl
        << "Total charge: " << nParticles * bunch->getQM(0) << " C" << endl;
    
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
    doSolve(myAmrOpal, bunch, rhs, phi, grad_phi, geom, rr, nLevels, msg);
    
    msg << "Total field energy: " << totalFieldEnergy(grad_phi, rr) << endl;
    
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
    
    
    writePlotFile(plotsolve, rhs, phi, grad_phi, rr, geom, 0);
    
//     dynamic_cast<AmrPartBunch*>(bunch)->python_format(0);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z] "
         << "[#levels] [max. box size] [h5file] [start:by:end]";
    
    if ( argc < 8 ) {
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
    std::string h5file = argv[6];
    std::string stepping = argv[7];

    std::string::size_type semipos1 = stepping.find_first_of(':');
    std::string::size_type semipos2 = stepping.find_first_of(':', semipos1 + 1);
    int begin = std::atoi( stepping.substr(0, semipos1).c_str() );
    int by = std::atoi( stepping.substr(semipos1 + 1, semipos2 - semipos1 - 1).c_str() );
    int end = std::atoi( stepping.substr(semipos2 + 1).c_str() );
    
    auto h5steps = std::make_tuple(begin, end, by);

    if ( begin < 0 || end < 0 ) {
	msg << "Negative values in range." << endl;
	return -1;
    } else if ( (begin > end && by >= 0) || (begin < end && by < 0) ) {
	msg << "Please check [start:by:end]. Neither upward nor "
	    << "downward stepping possible." << endl;
	return -1;
    } else if ( by && (end - begin) % by) {
	msg << "Stepping not a multiplicative of range." << endl;
	return -1;
    }
    
    msg << "Particle test running with" << endl
        << "- grid:         " << nr << endl
        << "- #level:       " << nLevels << endl
        << "- H5:           " << h5file << endl
	<< "- start:by:end: " << begin << ":" << by << ":" << end << endl;
    
    doBoxLib(nr, nLevels, maxBoxSize, h5file, h5steps, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    return 0;
}
