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

 mpirun -np 4 testAmrPartBunch IPPL 32 32 32 100 10
 
 mpirun -np 4 testAmrPartBunch BOXLIB 32 32 32 100 10 0

***************************************************************************/

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

#include "Distribution.h"
#include "Solver.h"
#include "AmrOpal.h"

#include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

// #define SOLVER


// typedef std::vector<std::unique_ptr<MultiFab> > container_t;
typedef Array<std::unique_ptr<MultiFab> > container_t;

double dt = 1.0;          // size of timestep


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
        grad_phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(myAmrOpal.boxArray()[lev],BL_SPACEDIM,1));

        rhs[lev]->setVal(0.0);
        phi[lev]->setVal(0.0);
        grad_phi[lev]->setVal(0.0);
    }
    

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();

    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rhs, base_level, 1, finest_level);
    
    // eps in C / (V * m)
    double constant = -1.0 / Physics::epsilon_0;
    for (int i = 0; i <= finest_level; ++i)
        rhs[i]->mult(constant, 0, 1);
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;
    if (geom[0].isAllPeriodic())
    {
        for (int lev = 0; lev < nLevels; lev++)
            offset = dynamic_cast<AmrPartBunch*>(bunch)->sumParticleMass(0,lev);
        offset /= geom[0].ProbSize();
    }

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


// ============================================================================
// IPPL
// ============================================================================
void doIppl(const Vektor<size_t, 3>& nr, size_t nParticles,
            size_t nTimeSteps, Inform& msg, Inform& msg2all)
{

    static IpplTimings::TimerRef distTimer = IpplTimings::getTimer("dist");
    static IpplTimings::TimerRef tracTimer = IpplTimings::getTimer("trac");    

    e_dim_tag decomp[Dim];
    unsigned serialDim = 2;

    msg << "Serial dimension is " << serialDim  << endl;

    Mesh_t *mesh;
    FieldLayout_t *FL;

    NDIndex<Dim> domain;
    for (unsigned i=0; i<Dim; i++)
        domain[i] = domain[i] = Index(nr[i] + 1);

    for (unsigned d=0; d < Dim; ++d)
        decomp[d] = (d == serialDim) ? SERIAL : PARALLEL;

    // create mesh and layout objects for this problem domain
    mesh          = new Mesh_t(domain);
    FL            = new FieldLayout_t(*mesh, decomp);
    playout_t* PL = new playout_t(*FL, *mesh);

    /*
     * In case of periodic BC's define
     * the domain with hr and rmin
     */
    Vector_t hr(1.0 / nr(0), 1.0 / nr(1), 1.0 / nr(2));
    Vector_t rmin(0.0);
    Vector_t rmax(1.0);
    
    PartBunchBase* bunch = new PartBunch<playout_t>(PL,hr,rmin,rmax,decomp);
        
    /*
     * initialize a particle distribution
     */
    IpplTimings::startTimer(distTimer);
    unsigned long int nloc = nParticles / Ippl::getNodes();
    Distribution dist;
    dist.uniform(0.2, 0.8, nloc, Ippl::myNode());
    dist.injectBeam(*bunch);
    // bunch->print();    
    bunch->myUpdate();
    IpplTimings::stopTimer(distTimer);

    double q = 1.0/nParticles;

    // random initialization for charge-to-mass ratio
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
        bunch->setQM(q, i);

    msg << "particles created and initial conditions assigned " << endl;

    // redistribute particles based on spatial layout

//     msg << "initial update and initial mesh done .... Q= " << sum(bunch->getQM()) << endl;
    msg << dynamic_cast<PartBunch<playout_t>*>(bunch)->getMesh() << endl;
    msg << dynamic_cast<PartBunch<playout_t>*>(bunch)->getFieldLayout() << endl;

    msg << "scatter test done delta= " <<  bunch->scatter() << endl;

    bunch->initFields();
    msg << "bunch->initField() done " << endl;

    // begin main timestep loop
    msg << "Starting iterations ..." << endl;
    IpplTimings::startTimer(tracTimer);
    for (unsigned int it=0; it<nTimeSteps; it++) {
        bunch->gatherStatistics();
//         // advance the particle positions
//         // basic leapfrogging timestep scheme.  velocities are offset
//         // by half a timestep from the positions.
        for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
            bunch->setR(bunch->getR(i) + dt * bunch->getP(i), i);

        // update particle distribution across processors
        bunch->myUpdate();

//         // gather the local value of the E field
//         bunch->gatherCIC();
// 
//         // advance the particle velocities
//         for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
//             bunch->setP(bunch->getP(i) + dt * bunch->getQM(i) * bunch->getE(i), i);
//         
//         msg << "Finished iteration " << it << " - min/max r and h " << bunch->getRMin()
//             << bunch->getRMax() << bunch->getHr() << endl;
    }
    Ippl::Comm->barrier();
    IpplTimings::stopTimer(tracTimer);
    delete bunch;
    msg << "Particle test testAmrPartBunch: End." << endl;
}



// ============================================================================
// BOXLIB
// ============================================================================
void doBoxLib(const Vektor<size_t, 3>& nr, size_t nParticles,
              int nLevels, size_t maxBoxSize,
              size_t nTimeSteps, Inform& msg, Inform& msg2all)
{
    static IpplTimings::TimerRef distTimer = IpplTimings::getTimer("dist");    
    static IpplTimings::TimerRef tracTimer = IpplTimings::getTimer("trac");
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
    
    domain.setLo(0, -0.025);
    domain.setHi(0,  0.025);
    domain.setLo(1, -0.025);
    domain.setHi(1,  0.025);
    domain.setLo(2, -0.025);
    domain.setHi(2,  0.025);
    
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
    msg << "Max. Grid Size = " << maxBoxSize << endl;
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
    unsigned long int nloc = nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    IpplTimings::startTimer(distTimer);
    dist.gaussian(0.0, 0.001, nloc, ParallelDescriptor::MyProc());
//     dist.uniform(-0.02, 0.02, nloc, ParallelDescriptor::MyProc());
    
//     dist.readH5("/home/matthias/Accelerated.h5", 0);
    
    
    // copy particles to the PartBunchBase object.
    dist.injectBeam(*bunch);
    
//     bunch->setR(Vector_t(0.003, 0.001, 0), 0);
    
    msg << "****************************************************" << endl
        << "            BEAM INJECTED (single level)            " << endl
        << "****************************************************" << endl;
    
    bunch->gatherStatistics();
    // redistribute on single-level
    bunch->myUpdate();
    IpplTimings::stopTimer(distTimer);
    
    msg << "****************************************************" << endl
        << "         BEAM REDISTRIBUTED (single level)          " << endl
        << "****************************************************" << endl;
    
    bunch->gatherStatistics();
    
    
    // ========================================================================
    // 2. tagging (i.e. create BoxArray's, DistributionMapping's of all
    //    other levels)
    // ========================================================================
    
    
    /*
     * create an Amr object
     */
    
    
//     dynamic_cast<AmrPartBunch*>(bunch)->Define (geom, dmap, ba, rr);
    
    ParmParse pp("amr");
    pp.add("max_grid_size", int(maxBoxSize));
    
    Array<int> error_buf(nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr[i];
    
    AmrOpal myAmrOpal(&domain, nLevels - 1, nCells, 0 /* cartesian */, bunch);
    
    for (int i = 0; i < nLevels; ++i)
        msg << "Max. grid size level" << i << ": "
            << myAmrOpal.maxGridSize(i) << endl;
    
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
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i) {
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    }
    
    msg << "****************************************************" << endl
        << "          BEAM REDISTRIBUTED (multi level)          " << endl
        << "****************************************************" << endl;
    
    bunch->gatherStatistics();
    
    
    
    msg << "****************************************************" << endl
        << "                   BEAM GEOMETRY                    " << endl
        << "****************************************************" << endl;
    for (int i = 0; i < nLevels; ++i)
        msg << dynamic_cast<AmrPartBunch*>(bunch)->GetParGDB()->Geom(i) << endl;
    
    
    msg << "****************************************************" << endl
        << "                   BEAM BOXARRAY                    " << endl
        << "****************************************************" << endl;
    for (int i = 0; i < nLevels; ++i)
        msg << dynamic_cast<AmrPartBunch*>(bunch)->GetParGDB()->boxArray(i) << endl;
    
    Inform amr("AMR");
    amr << "Max. level   = " << myAmrOpal.maxLevel() << endl
        << "Finest level = " << myAmrOpal.finestLevel() << endl;
    for (int i = 0; i < nLevels - 1; ++i)
        amr << "Max. ref. ratio level " << i << ": "
            << myAmrOpal.MaxRefRatio(i) << endl;
    for (int i = 0; i < nLevels; ++i)
        amr << "Max. grid size level" << i << ": "
            << myAmrOpal.maxGridSize(i) << endl;
    
    for (int i = 0; i < nLevels; ++i)
        amr << "BoxArray level" << i << ": "
            << myAmrOpal.boxArray(i) << endl;
    
    
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
        bunch->setP(Vector_t(1.0, 0.0, 0.0), i);         
    
    
    container_t rhs;
    container_t phi;
    container_t grad_phi;
    
    std::string plotsolve = BoxLib::Concatenate("solve_", 0, 4);
    doSolve(myAmrOpal, bunch, rhs, phi, grad_phi, geom, rr, nLevels);
    
    writePlotFile(plotsolve, rhs, phi, grad_phi, rr, geom, 0);
    
    // begin main timestep loop
    msg << "Starting iterations ..." << endl;
    IpplTimings::startTimer(tracTimer);
    for (unsigned int it=0; it<nTimeSteps; it++) {
        bunch->gatherStatistics();
        
//         
//         myAmrOpal.writePlotFile(plotfilename, it);
        myAmrOpal.assignDensity();
        
        
        std::string plotfilename = BoxLib::Concatenate("amr_", it, 4);
        myAmrOpal.writePlotFile(plotfilename, it);
        
        dynamic_cast<AmrPartBunch*>(bunch)->python_format(it);
        
        // update time step according to levels
//         dt = 0.9 * *( myAmrOpal.Geom(myAmrOpal.finestLevel() - 1).CellSize() );
        
        // advance the particle positions
//         for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
//             bunch->setR(bunch->getR(0) + dt * bunch->getP(0), 0);
        
        dist.setDistribution(*bunch, "/home/matthias/Accelerated.h5", it + 1);
        
        bunch->myUpdate();
        
        
        // update Amr object
        // update particle distribution across processors
        for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i) {
            myAmrOpal.regrid(i /*lbase*/, Real(it) /*time*/);
            msg << "DONE " << i << "th regridding." << endl;
        }
        
        
//         doSolve(myAmrOpal, bunch, rhs, phi, grad_phi, geom, rr, nLevels);
    }
    ParallelDescriptor::Barrier();
    IpplTimings::stopTimer(tracTimer);

    delete bunch;
    msg << "Particle test testAmrPartBunch: End." << endl;
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg(argv[1]);
    Inform msg2all(argv[1], INFORM_ALL_NODES);
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] testAmrPartBunch [IPPL or BOXLIB] "
         << "[#gridpoints x] [#gridpoints y] [#gridpoints z] [#particles] [#timesteps] ";
    
    if ( argc < 7 ) {
        msg << call.str() << endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<size_t, 3> nr(std::atoi(argv[2]),
                         std::atoi(argv[3]),
                         std::atoi(argv[4]));
    
    
    size_t nParticles = std::atoi(argv[5]);
    size_t nTimeSteps = std::atoi(argv[6]);
    
    
    msg << "Particle test running with" << endl
        << "- #timesteps = " << nTimeSteps << endl
        << "- #particles = " << nParticles << endl
        << "- grid       = " << nr << endl;
    
    
    if ( std::strcmp(argv[1], "IPPL") ) {
        if ( argc != 9 ) {
            msg << call.str() << "[#levels] [max. box size]" << endl;
            return -1;
        }
        
        BoxLib::Initialize(argc,argv, false);
        size_t nLevels = std::atoi(argv[7]) + 1; // i.e. nLevels = 0 --> only single level
        size_t maxBoxSize = std::atoi(argv[8]);
        doBoxLib(nr, nParticles, nLevels, maxBoxSize, nTimeSteps, msg, msg2all);
    } else
        doIppl(nr, nParticles, nTimeSteps, msg, msg2all);
    

    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    std::string tfn = std::string(argv[1]) + "-cores=" + std::to_string(Ippl::getNodes());
    IpplTimings::print(tfn);

    return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $
 ***************************************************************************/
