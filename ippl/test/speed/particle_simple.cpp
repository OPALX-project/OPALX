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
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

/***************************************************************************
 *
 * Benchmark to measure time required for particle create/delete
 *
 * Usage: progname  <num particles>  <vnodes/processor>
 *
 ***************************************************************************/

#include "Ippl.h"
#include <stdlib.h>

// constants for this benchmark
const unsigned int Iterations =  100;
const unsigned int Dim        =    2;
const unsigned int XSize      =  100;
const unsigned int YSize      =  100;

// typedefs for this benchmark
typedef double                            PPos_t;
typedef double                            PData_t;
typedef ParticleSpatialLayout<PPos_t,Dim> PLayout_t;
typedef GenParticle<PLayout_t,PData_t>    Particle_t;


//////////////////////////////////////////////////////////////////////

void usage(int argc, char *argv[]) {
  Inform msg("Error");
  msg << "Usage: " << argv[0] << " <num particles> <vnodes/processor>";
  msg << endl;
  exit(1);
}


//////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[] )
{
  int i, N, vpp, npp;
  Ippl ippl(argc,argv);
  Inform msg(argv[0], INFORM_ALL_NODES);
  Inform dbg("pdebug");

  // get parameters for the benchmark
  if (argc != 3)
    usage(argc,argv);
  if ( (N=atoi(argv[1])) < 1)
    usage(argc,argv);
  if ( (vpp=atoi(argv[2])) < 1)
    usage(argc,argv);
  vpp *= Ippl::getNodes();
  npp = N / Ippl::getNodes();

  // create a particle object, and initialize it
  FieldLayout<Dim> flayout(Index(XSize),Index(YSize),PARALLEL,PARALLEL,vpp);
  Particle_t  P(new PLayout_t(flayout));
  P.getBConds()[0] = ParticlePeriodicBCond;
  P.getBConds()[1] = ParticlePeriodicBCond;
  P.getBConds()[2] = ParticlePeriodicBCond;
  P.getBConds()[3] = ParticlePeriodicBCond;
  msg << "Created new particle object = " << P << endl;

  // perform benchmark: for each iteration, create a number of particles
  // (the same number on each node), assign random values to them, then
  // do an update (which will require particles to be swapped among nodes).
  // After this, delete about 1/2 the particles, and just do a delete
  // with no swap.
  // Measure the total time required to do each step for all the iterations,
  // and report the average time per particle per iteration.
  Timer tippl_create;
  Timer tippl_assign;
  Timer tippl_update;
  Timer tippl_delete;
  double usec_times[4], glob_times[4];
  usec_times[0] = glob_times[0] = 0.0;
  usec_times[1] = glob_times[1] = 0.0;
  usec_times[2] = glob_times[2] = 0.0;
  usec_times[3] = glob_times[3] = 0.0;

  msg << "Starting " << Iterations << " iterations ..." << endl;
  for (i=0; i < Iterations; ++i) {
    // create N particles
    tippl_create.clear();
    tippl_create.start();
    P.create(npp);
    tippl_create.stop();
    usec_times[0] += tippl_create.clock_time();

    // initialize them with random positions and data
    tippl_assign.clear();
    tippl_assign.start();
    assign(P.R(0), IpplRandom * (double)XSize);
    assign(P.R(1), IpplRandom * (double)YSize);
    assign(P.data, 10.0 * IpplRandom);
    tippl_assign.stop();
    usec_times[1] += tippl_assign.clock_time();

    // do an update, which will swap the particle positions and update the
    // count of the total number of particles on each node
    tippl_update.clear();
    tippl_update.start();
    P.update();
    tippl_update.stop();
    usec_times[2] += tippl_update.clock_time();

    if (i % 10 == 0) {
      dbg << "iter " << i << ": " << usec_times[0] << ", " << usec_times[1];
      dbg << ", " << usec_times[2] << ", " << usec_times[3] << ", " << endl;
      dbg << "  "; Message::printDebug(dbg); dbg << endl;
      dbg << "  "; P.printDebug(dbg);        dbg << endl;
    }

    // delete some of the particles, and then just do an update to change
    // the particle numbers (do not swap)
    int currp = P.getLocalNum();
    int maxdp = 0;
    while (maxdp < (currp/2 - 1)) {
      if (maxdp % 3 != 0)
	P.destroy(1, maxdp);
      maxdp++;
    }
    P.setUpdateFlag(PLayout_t::SWAP, false);
    tippl_delete.clear();
    tippl_delete.start();
    P.update();
    tippl_delete.stop();
    usec_times[3] += tippl_delete.clock_time();

    // delete all the other particles, so we start out fresh
    P.destroy(P.getLocalNum(), 0);
    P.update();
    P.setUpdateFlag(PLayout_t::SWAP, true);
  }
  msg << "Iterations complete.\n" << endl;

  // at the end, get the total times, and take the average over all nodes
  msg << "Node " << Ippl::myNode() << " times (c,a,u,d): ";
  msg << usec_times[0] << ", " << usec_times[1] << ", " << usec_times[2];
  msg << ", " << usec_times[3] << endl;

  reduce(usec_times, usec_times + 4, glob_times, OpAddAssign());

  double scalefactor = 1.e6 / ((double)(Iterations * N * Ippl::getNodes()));
  for (i=0; i < 4; ++i)
    usec_times[i] = glob_times[i] * scalefactor;

  Inform m2("results");
  m2 << "Create CPU time/iter/particle/node (usec)= " << usec_times[0] << endl;
  m2 << "Assign CPU time/iter/particle/node (usec)= " << usec_times[1] << endl;
  m2 << "Update CPU time/iter/particle/node (usec)= " << usec_times[2] << endl;
  m2 << "Delete CPU time/iter/particle/node (usec)= " << usec_times[3] << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: particle_simple.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:40 $
 * IPPL_VERSION_ID: $Id: particle_simple.cpp,v 1.1.1.1 2003/01/23 07:40:40 adelmann Exp $ 
 ***************************************************************************/
