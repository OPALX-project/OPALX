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

#define SOLVER 1


const double dt = 1.0;          // size of timestep


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


void doBoxLib(const Vektor<size_t, 3>& nr, size_t nParticles,
              size_t nLevels, size_t maxBoxSize,
              size_t nTimeSteps, Inform& msg, Inform& msg2all)
{
    static IpplTimings::TimerRef distTimer = IpplTimings::getTimer("dist");    
    static IpplTimings::TimerRef tracTimer = IpplTimings::getTimer("trac");    
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
    int bc[BL_SPACEDIM] = {0, 0, 0};
    
    Array<Geometry> geom(nLevels);
    
    // level 0 describes physical domain
    geom[0].define(bx, &domain, 0, bc);
    
    // Container for boxes at all levels
    Array<BoxArray> ba(nLevels);
    
    // box at level 0
    ba[0].define(bx);
    ba[0].maxSize(maxBoxSize);
    
    /*
     * distribution mapping
     */
    Array<DistributionMapping> dmap(nLevels);
    dmap[0].define(ba[0], ParallelDescriptor::NProcs() /*nprocs*/);
    
    
    /*
     * refinement ratio
     */
    Array<int> rr(nLevels - 1);
    for (size_t lev = 1; lev < nLevels; ++lev)
        rr[lev - 1] = 2;
    
    
    for (size_t lev = 1; lev < nLevels; ++lev) {
        geom[lev].define(BoxLib::refine(geom[lev - 1].Domain(),
                                        rr[lev - 1]),
                         &domain, 0, bc);
    }
    
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    PartBunchBase* bunch = new AmrPartBunch(geom, dmap, ba, rr);
    
    
    // initialize a particle distribution
    unsigned long int nloc = nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    IpplTimings::startTimer(distTimer);
    dist.uniform(0.2, 0.8, nloc, ParallelDescriptor::MyProc());
    // copy particles to the PartBunchBase object.
    dist.injectBeam(*bunch);
    bunch->gatherStatistics();
    // redistribute on single-level
    bunch->myUpdate();
    IpplTimings::stopTimer(distTimer);
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
    
    ParmParse pp("amr");
    pp.add("max_level", int(nLevels - 1));
    pp.add("ref_ratio", int(2)); //FIXME
    pp.add("max_grid_size", int(maxBoxSize));
    
    std::vector<int> tmp(3); tmp[0] = nr[0]; tmp[1] = nr[1]; tmp[2] = nr[2];
    pp.addarr("n_cell", tmp);//IntVect(nr[0], nr[1], nr[2]));
    
    AmrOpal myAmrOpal;
    
    myAmrOpal.SetBoxArray(0, ba[0]);
    
    MultiFab nPartPerCell;
    TagBoxArray tags;
    Real time = 0.0;
    int tag_level = 0;
    
    myAmrOpal.ErrorEst(tag_level, nPartPerCell, tags, time);
    
    ///@todo Next higher boxes are hard-coded.
    int fine = 1.0;
    for (int lev = 1; lev < int(nLevels); ++lev) {
        fine *= rr[lev - 1];
        
        if ( lev == 1) {
            IntVect refined_lo(0.25 * nr[0] * fine,
                                0.25 * nr[1] * fine,
                                0.25 * nr[2] * fine);
            
            IntVect refined_hi(0.75 * nr[0] * fine - 1,
                                0.75 * nr[1] * fine - 1,
                                0.75 * nr[2] * fine - 1);
        Box refined_patch(refined_lo, refined_hi);
        ba[lev].define(refined_patch);
        } else if ( lev == 2) {
            IntVect refined_lo(0.375 * nr[0] * fine,
                                0.375* nr[1] * fine,
                                0.375 * nr[2] * fine);
            
            IntVect refined_hi(0.625 * nr[0] * fine - 1,
                               0.625 * nr[1] * fine - 1,
                                0.625* nr[2] * fine - 1);

            Box refined_patch(refined_lo, refined_hi);
            ba[lev].define(refined_patch);
        }
    }
    
    // break BoxArrays into maxBoxSize^3
    for (size_t lev = 1; lev < nLevels; ++lev) {
        ba[lev].maxSize(maxBoxSize);
        dmap[lev].define(ba[lev], ParallelDescriptor::NProcs() /*nprocs*/);
    }
    
    for (size_t lev = 1; lev < nLevels; ++lev) {
        dynamic_cast<AmrPartBunch*>(bunch)->SetParticleBoxArray(lev, /*dmap[lev],*/ ba[lev]);
        dynamic_cast<AmrPartBunch*>(bunch)->SetParticleDistributionMap(lev, dmap[lev]);
    }
    
    
    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    
    bunch->myUpdate();
    
    
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

    for (size_t lev = 0; lev < nLevels; ++lev) {
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
    int finest_level = nLevels - 1;

    PArray<MultiFab> PartMF;
    PartMF.resize(nLevels,PArrayManage);
    PartMF.set(0,new MultiFab(ba[0],1,1));
    PartMF[0].setVal(0.0);
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, PartMF, base_level, 1, finest_level);
    //MyPC->AssignDensitySingleLevel(0, PartMF[0],0,1,0);                                                                                                                                                        

    for (int lev = finest_level - 1 - base_level; lev >= 0; lev--)
        BoxLib::average_down(PartMF[lev+1],PartMF[lev],0,1,rr[lev]);

    for (size_t lev = 0; lev < nLevels; lev++)
        MultiFab::Add(rhs[base_level+lev], PartMF[lev], 0, 0, 1, 0);
    
    
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;
//     if (geom[0].isAllPeriodic())
//     {
//         for (size_t lev = 0; lev < nLevels; lev++)
//             offset = dynamic_cast<AmrPartBunch*>(bunch)->sumParticleMass(0,lev);
//         offset /= geom[0].ProbSize();
//     }

    // solve                                                                                                                                                                                                     
    Solver sol;
    sol.solve_for_accel(rhs,phi,grad_phi,geom,base_level,finest_level,offset);
#endif
    // ========================================================================
    // do some operations on the data
    // ========================================================================
    
    MultiFab mf(ba[0], 1, 1);
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensitySingleLevel(0, /*attribute*/ mf, 0 /*level*/);
    
    
//     std::cout << mf.max(0) << std::endl << mf.min(0) << std::endl;
    
    Real charge = dynamic_cast<AmrPartBunch*>(bunch)->sumParticleMass(0 /*attribute*/, 0 /*level*/);
    
    double invVol = 1.0 / (nr[0] * nr[1] * nr[2]);
    
    std::cout << "MultiFab sum: " << mf.sum() * invVol << std::endl
              << "Charge sum: " << charge << std::endl;
    
    
    
    // begin main timestep loop
    msg << "Starting iterations ..." << endl;
    IpplTimings::startTimer(tracTimer);
    for (unsigned int it=0; it<nTimeSteps; it++) {
        bunch->gatherStatistics();
        // advance the particle positions
        for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
            bunch->setR(bunch->getR(i) + dt * bunch->getP(i), i);

        // update particle distribution across processors
        bunch->myUpdate();
        
        dynamic_cast<AmrPartBunch*>(bunch)->WritePlotFile(".", "boxlib-timestep-" + std::to_string(it));

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
        size_t nLevels = std::atoi(argv[7]);
        size_t maxBoxSize = std::atoi(argv[8]);
        doBoxLib(nr, nParticles, nLevels, maxBoxSize, nTimeSteps, msg, msg2all);
    } else
        doIppl(nr, nParticles, nTimeSteps, msg, msg2all);
    

    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    std::string tfn = std::string(argv[1]) + std::string("-cores=") + std::to_string(Ippl::getNodes());
    IpplTimings::print(tfn);

    return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $
 ***************************************************************************/
