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

#include <cmath>

// #define SOLVER


double dt = 1.0;          // size of timestep


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
              size_t nMaxLevels, size_t maxBoxSize,
              size_t nTimeSteps, Inform& msg, Inform& msg2all)
{
    static IpplTimings::TimerRef distTimer = IpplTimings::getTimer("dist");    
    static IpplTimings::TimerRef tracTimer = IpplTimings::getTimer("trac");
    static IpplTimings::TimerRef solvTimer = IpplTimings::getTimer("solv");
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    /*
     * set up the geometry
     */
    IntVect low(0, 0, 0);
    IntVect high(nr[0] - 1, nr[1] - 1, nr[2] - 1);    
    Box bx(low, high);
    
    // box [-1,1]x[-1,1]x[-1,1]
    RealBox domain;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain.setLo(i, 0.0);
        domain.setHi(i, 1.0);
    }
    
    
    // periodic boundary conditions in all directions
    int bc[BL_SPACEDIM] = {1, 1, 1};
    
    
    Geometry geom;
    
    // level 0 describes physical domain
    geom.define(bx, &domain, 0, bc);
    
    // Container for boxes at all levels
    BoxArray ba;
    
    // box at level 0
    ba.define(bx);
    msg << "Max. Grid Size = " << maxBoxSize << endl;
    ba.maxSize(maxBoxSize);
    
    /*
     * distribution mapping
     */
    DistributionMapping dmap;
    dmap.define(ba, ParallelDescriptor::NProcs() /*nprocs*/);
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    PartBunchBase* bunch = new AmrPartBunch(geom, dmap, ba);
    
    
    // initialize a particle distribution
    unsigned long int nloc = nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    IpplTimings::startTimer(distTimer);
    dist.uniform(0.3, 0.7, nloc, ParallelDescriptor::MyProc());
    // copy particles to the PartBunchBase object.
    dist.injectBeam(*bunch);
    
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
    //     bunch->print();
    
    double q = 1.0 / nParticles;

    // random initialization for charge-to-mass ratio
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
        bunch->setQM(q, i);
    
    
    // ========================================================================
    // 2. tagging (i.e. create BoxArray's, DistributionMapping's of all
    //    other levels)
    // ========================================================================
    
    
    /*
     * create an Amr object
     */
    
    ParmParse pp("amr");
    pp.add("max_grid_size", int(maxBoxSize));
    
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr[i];
    
    // nLevel = nMaxLevels + 1
    int nLevels = nMaxLevels + 1;
    AmrOpal myAmrOpal(&domain, int(nMaxLevels), nCells, 0 /* cartesian */, bunch);
    
    for (int i = 0; i < nLevels; ++i)
        msg << "Max. grid size level" << i << ": "
            << myAmrOpal.maxGridSize(i) << endl;
    
    /*
     * do tagging
     */
    Array<int> rr(nMaxLevels);
    for (int lev = 1; lev < nLevels; ++lev)
        rr[lev-1] = 2;
    dynamic_cast<AmrPartBunch*>(bunch)->Define (myAmrOpal.Geom(),
                                                myAmrOpal.DistributionMap(),
                                                myAmrOpal.boxArray(),
                                                rr);
    
    
    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    for (int i = 0; i < myAmrOpal.maxLevel(); ++i) {
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
        msg << "DONE " << i << "th regridding." << endl;
    }
    
    
    myAmrOpal.info();
//     myAmrOpal.averageDown();
//     myAmrOpal.finestDensityAssign();
//     bunch->myUpdate();
    
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
    
    
    
//     // print some information
//     myAmrOpal.writePlotFile("amr0000");
//     myAmrOpal.free();
    
    Inform amr("AMR");
    amr << "Max. level   = " << myAmrOpal.maxLevel() << endl
        << "Finest level = " << myAmrOpal.finestLevel() << endl;
    for (int i = 0; i < int(nMaxLevels); ++i)
        amr << "Max. ref. ratio level " << i << ": "
            << myAmrOpal.MaxRefRatio(i) << endl;
    for (int i = 0; i < nLevels; ++i)
        amr << "Max. grid size level" << i << ": "
            << myAmrOpal.maxGridSize(i) << endl;
    
    for (int i = 0; i < nLevels; ++i)
        amr << "BoxArray level" << i << ": "
            << myAmrOpal.boxArray(i) << endl;
    
    
//     // write particle file
//     dynamic_cast<AmrPartBunch*>(bunch)->Checkpoint(".", "particles0000", true);
    
    
    // =======================================================================                                                                                                                                   
    // 4. prepare for multi-level solve                                                                                                                                                                          
    // =======================================================================
#ifdef SOLVER
    PArray<MultiFab> rhs;
    PArray<MultiFab> phi;
    PArray<MultiFab> grad_phi;

    rhs.resize(nLevels,PArrayManage);
    phi.resize(nLevels,PArrayManage);
    grad_phi.resize(nLevels,PArrayManage);

    for (int lev = 0; lev < nLevels; ++lev) {
        //                                    # component # ghost cells                                                                                                                                          
        rhs.set     (lev,new MultiFab(ba[lev],1          ,0));
        phi.set     (lev,new MultiFab(ba[lev],1          ,1));
        grad_phi.set(lev,new MultiFab(ba[lev],BL_SPACEDIM,1));

        rhs[lev].setVal(0.0);
        phi[lev].setVal(0.0);
        grad_phi[lev].setVal(0.0);
    }
    

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = nMaxLevels;

    PArray<MultiFab> PartMF;
    PartMF.resize(nLevels,PArrayManage);
    PartMF.set(0,new MultiFab(ba[0],1,1));
    PartMF[0].setVal(0.0);
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, PartMF, base_level, 1, finest_level);

    for (int lev = finest_level - 1 - base_level; lev >= 0; lev--)
        BoxLib::average_down(PartMF[lev+1],PartMF[lev],0,1,rr[lev]);

    for (int lev = 0; lev < nLevels; lev++)
        MultiFab::Add(rhs[base_level+lev], PartMF[lev], 0, 0, 1, 0);
    
    
    
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
    sol.solve_for_accel(rhs,phi,grad_phi,geom,base_level,finest_level,offset);
    IpplTimings::stopTimer(solvTimer);
#endif
    // ========================================================================
    // do some operations on the data
    // ========================================================================
    
    
    for (int i = 0; i < myAmrOpal.finestLevel() + 1; ++i) {
        
        if ( myAmrOpal.boxArray(i).empty() )
            break;
        
        MultiFab mf(myAmrOpal.boxArray(i), 1, 1);
        dynamic_cast<AmrPartBunch*>(bunch)->AssignDensitySingleLevel(0, /*attribute*/ mf, i /*level*/);
        Real charge = dynamic_cast<AmrPartBunch*>(bunch)->sumParticleMass(0 /*attribute*/, i /*level*/);
        
        double invVol = 1.0 / (nr[0] * nr[1] * nr[2]);
        
        /* refinement factor */
        invVol = (i > 0) ? invVol / std::pow( ( 2 << (i - 1) ), 3) : invVol;
        
        std::cout << "MultiFab sum: " << mf.sum() * invVol << std::endl
                  << "Charge sum: " << charge << std::endl;

        
    }
    
    std::string plotfilename = BoxLib::Concatenate("amr", 0, 4);
    myAmrOpal.writePlotFile(plotfilename, 0);
    
    
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
        bunch->setP(Vector_t(1.0, 0.0, 0.0), i);         
    
    // begin main timestep loop
    msg << "Starting iterations ..." << endl;
    IpplTimings::startTimer(tracTimer);
    for (unsigned int it=0; it<nTimeSteps; it++) {
        bunch->gatherStatistics();
        
        std::string plotfilename = BoxLib::Concatenate("amr_", it, 4);
        myAmrOpal.writePlotFile(plotfilename, it);
        
        // update time step according to levels
        dt = 0.5 * *( myAmrOpal.Geom(myAmrOpal.finestLevel() - 1).CellSize() );
        
        // advance the particle positions
        for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
            bunch->setR(bunch->getR(i) + dt * bunch->getP(i), i);
        
        
        // update Amr object
        // update particle distribution across processors
        for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i) {
            myAmrOpal.regrid(i /*lbase*/, Real(it) /*time*/);
            msg << "DONE " << i << "th regridding." << endl;
        }
        
        //
        for (int i = 0; i < myAmrOpal.finestLevel() + 1; ++i) {
            MultiFab mf(myAmrOpal.boxArray(i), 1, 1);
            dynamic_cast<AmrPartBunch*>(bunch)->AssignDensitySingleLevel(0, /*attribute*/ mf, i /*level*/);
            Real charge = dynamic_cast<AmrPartBunch*>(bunch)->sumParticleMass(0 /*attribute*/, i /*level*/);
            
            double invVol = 1.0 / ( nr[0] * nr[1] * nr[2] );
            
            /* refinement factor */
            invVol = (i > 0) ? invVol / std::pow( ( 2 << (i - 1) ), 3) : invVol;
            
            std::cout << "MultiFab sum (not normalized): " << mf.sum() << std::endl;
            std::cout << "MultiFab sum: " << mf.sum() * invVol << std::endl
                      << "Charge sum: " << charge << std::endl;

        
        }
        //
        
//         myAmrOpal.averageDown();
        
        amr << "Max. level   = " << myAmrOpal.maxLevel() << endl
            << "Finest level = " << myAmrOpal.finestLevel() << endl;
        for (int i = 0; i < int(nMaxLevels); ++i)
            amr << "Max. ref. ratio level " << i << ": "
                << myAmrOpal.MaxRefRatio(i) << endl;
        for (int i = 0; i < nLevels; ++i)
            amr << "Max. grid size level" << i << ": "
                << myAmrOpal.maxGridSize(i) << endl;
        
        for (int i = 0; i < nLevels; ++i)
            amr << "BoxArray level" << i << ": "
                << myAmrOpal.boxArray(i) << endl;

        msg << "Done timestep " << it << " using dt = " << dt << endl;
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
        size_t nMaxLevels = std::atoi(argv[7]);
        size_t maxBoxSize = std::atoi(argv[8]);
        doBoxLib(nr, nParticles, nMaxLevels, maxBoxSize, nTimeSteps, msg, msg2all);
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
