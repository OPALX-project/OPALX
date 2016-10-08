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


#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>

#include "PartBunch.h"
#include "AmrPartBunch.h"


const double dt = 1.0;          // size of timestep


void doIppl(const Vektor<size_t, 3>& nr, size_t nParticles,
            size_t nTimeSteps, Inform& msg, Inform& msg2all)
{

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
    Vector_t hr(1.0);
    Vector_t rmin(0.0);
    Vector_t rmax(nr);
    
    PartBunchBase* bunch = new PartBunch<playout_t>(PL,hr,rmin,rmax,decomp);

    // initialize the particle object: do all initialization on one node,
    // and distribute to others

    unsigned long int nloc = nParticles / Ippl::getNodes();

    bunch->create(nloc);
    for (unsigned long int i = 0; i< nloc; i++) {
        for (int d = 0; d<3; d++)
            bunch->getR(i)(d) =  IpplRandom() * nr[d];
    }

    double q = 1.0/nParticles;

    // random initialization for charge-to-mass ratio
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
        bunch->getQM(i) = q;

    msg << "particles created and initial conditions assigned " << endl;

    // redistribute particles based on spatial layout
    bunch->myUpdate();

//     msg << "initial update and initial mesh done .... Q= " << sum(bunch->getQM()) << endl;
    msg << dynamic_cast<PartBunch<playout_t>*>(bunch)->getMesh() << endl;
    msg << dynamic_cast<PartBunch<playout_t>*>(bunch)->getFieldLayout() << endl;

    msg << "scatter test done delta= " <<  bunch->scatter() << endl;

    bunch->initFields();
    msg << "bunch->initField() done " << endl;

    // begin main timestep loop
    msg << "Starting iterations ..." << endl;
    for (unsigned int it=0; it<nTimeSteps; it++) {
        bunch->gatherStatistics();
        // advance the particle positions
        // basic leapfrogging timestep scheme.  velocities are offset
        // by half a timestep from the positions.
        for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
            bunch->getR(i) += dt * bunch->getP(i);

        // update particle distribution across processors
        bunch->myUpdate();

        // gather the local value of the E field
        bunch->gatherCIC();

        // advance the particle velocities
        for (unsigned int i = 0; i < bunch->getLocalNum(); ++i)
            bunch->getP(i) += dt * bunch->getQM(i) * bunch->getE(i);
        
        msg << "Finished iteration " << it << " - min/max r and h " << bunch->getRMin()
            << bunch->getRMax() << bunch->getHr() << endl;
    }
    Ippl::Comm->barrier();
    
    delete bunch;
    msg << "Particle test testAmrPartBunch: End." << endl;
}


void doBoxLib(const Vektor<size_t, 3>& nr, size_t nParticles,
              size_t maxLevel, size_t nTimeSteps, Inform& msg,
              Inform& msg2all)
{
    /*
     * set up the geometry
     */
    IntVect low(0, 0, 0);
    IntVect high(nr[0] - 1, nr[1] - 1, nr[2] - 1);    
    Box bx(low, high);
    
    BoxArray ba(bx);
    
    ba.maxSize(16);
    
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
    
    msg << geom << endl;
    
    /*
     * initialize particle bunch
     */
    
    PartBunchBase* bunch = new AmrPartBunch(geom, dmap, ba);
    
    RealBox particle_box;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        particle_box.setLo(i, 0.1);
        particle_box.setHi(i, 0.9);
    }
    
    dynamic_cast<AmrPartBunch*>(bunch)->InitRandom(nParticles,
                      42 /*seed*/,
                      1 /* particle mass*/,
                      false /*serialize*/,
                      particle_box);
    
    
    bunch->myUpdate();
    
    for (size_t i = 0; i < bunch->getLocalNum(); ++i)
        msg2all << bunch->getR(i) << endl;
    
    
    MultiFab mf(ba, 1, 1);
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensitySingleLevel(mf, 0);
    
    
    std::cout << mf.max(0) << std::endl << mf.min(0) << std::endl;
    
    Real charge = dynamic_cast<AmrPartBunch*>(bunch)->sumParticleMass(0);
    
    double invVol = 1.0 / (nr[0] * nr[1] * nr[2]);
    
    std::cout << "MultiFab sum: " << mf.sum() * invVol << std::endl
              << "Charge sum: " << charge << std::endl;
    
    delete bunch;
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg(argv[1]);
    Inform msg2all(argv[1], INFORM_ALL_NODES);
    
    
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
        BoxLib::Initialize(argc,argv, false);
        size_t maxLevel = std::atoi(argv[7]);
        doBoxLib(nr, nParticles, maxLevel, nTimeSteps, msg, msg2all);
    } else
        doIppl(nr, nParticles, nTimeSteps, msg, msg2all);
    
    return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $
 ***************************************************************************/
