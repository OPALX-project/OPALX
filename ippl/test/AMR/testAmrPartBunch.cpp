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

#include "PartBunch.h"
#include "AmrPartBunch.h"

#include "Distribution.h"


const double dt = 1.0;          // size of timestep


void doIppl(const Vektor<size_t, 3>& nr, size_t nParticles,
            size_t nTimeSteps, Inform& msg, Inform& msg2all)
{

    e_dim_tag decomp[Dim];
    unsigned serialDim = 2;

//     msg << "Serial dimension is " << serialDim  << endl;

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
    unsigned long int nloc = nParticles / Ippl::getNodes();
    Distribution dist;
    dist.uniform(0.0, 1.0, nloc, Ippl::myNode());
    dist.injectBeam(*bunch);
    
    bunch->print();
    
    bunch->myUpdate();
    
    double q = 1.0/nParticles;

    // random initialization for charge-to-mass ratio
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
        bunch->setQM(q, i);

//     msg << "particles created and initial conditions assigned " << endl;

    // redistribute particles based on spatial layout

//     msg << "initial update and initial mesh done .... Q= " << sum(bunch->getQM()) << endl;
//     msg << dynamic_cast<PartBunch<playout_t>*>(bunch)->getMesh() << endl;
//     msg << dynamic_cast<PartBunch<playout_t>*>(bunch)->getFieldLayout() << endl;

//     msg << "scatter test done delta= " <<  bunch->scatter() << endl;

    bunch->initFields();
//     msg << "bunch->initField() done " << endl;

    // begin main timestep loop
    msg << "Starting iterations ..." << endl;
    for (unsigned int it=0; it<nTimeSteps; it++) {
        bunch->gatherStatistics();
        // advance the particle positions
        // basic leapfrogging timestep scheme.  velocities are offset
        // by half a timestep from the positions.
        for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
            bunch->setR(bunch->getR(i) + dt * bunch->getP(i), i);

        // update particle distribution across processors
        bunch->myUpdate();

        // gather the local value of the E field
        bunch->gatherCIC();

        // advance the particle velocities
        for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
            bunch->setP(bunch->getP(i) + dt * bunch->getQM(i) * bunch->getE(i), i);
        
        msg << "Finished iteration " << it << " - min/max r and h " << bunch->getRMin()
            << bunch->getRMax() << bunch->getHr() << endl;
    }
    Ippl::Comm->barrier();
    
    delete bunch;
    msg << "Particle test testAmrPartBunch: End." << endl;
}


void doBoxLib(const Vektor<size_t, 3>& nr, size_t nParticles,
              size_t maxLevel, size_t maxBoxSize,
              size_t nTimeSteps, Inform& msg, Inform& msg2all)
{
    /*
     * set up the geometry
     */
    IntVect low(0, 0, 0);
    IntVect high(nr[0] - 1, nr[1] - 1, nr[2] - 1);    
    Box bx(low, high);
    
    BoxArray ba(bx);
    
    ba.maxSize(maxBoxSize);
    
    // box [-1,1]x[-1,1]x[-1,1]
    RealBox domain;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain.setLo(i, 0.0);
        domain.setHi(i, 1.0);
    }
    
    
    // periodic boundary conditions in all directions
    int bc[BL_SPACEDIM] = {1, 1, 1};
    
    Geometry geom;
    geom.define(bx, &domain, 0, bc);
    
    DistributionMapping dmap;
    dmap.define(ba, ParallelDescriptor::NProcs() /*nprocs*/);
    
//     msg << geom << endl;
    
    /*
     * initialize particle bunch
     */
    
    PartBunchBase* bunch = new AmrPartBunch(geom, dmap, ba);
    
    
    /*
     * initialize a particle distribution
     */
    unsigned long int nloc = nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    dist.uniform(0.0, 1.0, nloc, ParallelDescriptor::MyProc());
    dist.injectBeam(*bunch);
    
    bunch->myUpdate();
    bunch->print();
    
    double q = 1.0 / nParticles;

    // random initialization for charge-to-mass ratio
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
        bunch->setQM(q, i);
    
    
    
    
    MultiFab mf(ba, 1, 1);
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensitySingleLevel(0, /*attribute*/ mf, 0 /*level*/);
    
    
//     std::cout << mf.max(0) << std::endl << mf.min(0) << std::endl;
    
    Real charge = dynamic_cast<AmrPartBunch*>(bunch)->sumParticleMass(0 /*attribute*/, 0 /*level*/);
    
    double invVol = 1.0 / (nr[0] * nr[1] * nr[2]);
    
    std::cout << "MultiFab sum: " << mf.sum() * invVol << std::endl
              << "Charge sum: " << charge << std::endl;
    
    
    
    // begin main timestep loop
    msg << "Starting iterations ..." << endl;
    for (unsigned int it=0; it<nTimeSteps; it++) {
        bunch->gatherStatistics();
        // advance the particle positions
        // basic leapfrogging timestep scheme.  velocities are offset
        // by half a timestep from the positions.
        for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
            bunch->setR(bunch->getR(i) + dt * bunch->getP(i), i);

        // update particle distribution across processors
        bunch->myUpdate();

        // gather the local value of the E field
        bunch->gatherCIC();

        // advance the particle velocities
        for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
            bunch->setP(bunch->getP(i) + dt * bunch->getQM(i) * bunch->getE(i), i);
        
        msg << "Finished iteration " << it << " - min/max r and h " << bunch->getRMin()
            << bunch->getRMax() << bunch->getHr() << endl;
    }
    ParallelDescriptor::Barrier();
    
    delete bunch;
    msg << "Particle test testAmrPartBunch: End." << endl;
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg(argv[1]);
    Inform msg2all(argv[1], INFORM_ALL_NODES);
    
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
            msg << call.str() << "[max. level] [max. box size]" << endl;
            return -1;
        }
        
        BoxLib::Initialize(argc,argv, false);
        size_t maxLevel = std::atoi(argv[7]);
        size_t maxBoxSize = std::atoi(argv[8]);
        doBoxLib(nr, nParticles, maxLevel, maxBoxSize, nTimeSteps, msg, msg2all);
    } else
        doIppl(nr, nParticles, nTimeSteps, msg, msg2all);
    
    return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $
 ***************************************************************************/
