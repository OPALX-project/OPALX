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
 * Particle Boundary Conditions test ... initializes a set of 2D particles,
 * changes their position to move beyond the boundary, and prints out the
 * resulting positions.
 *
 * This version includes the commands to instrument the code for connection
 * to another application or visualization component.
 ***************************************************************************/

#include "Ippl.h"

// dimension of our positions
const unsigned Dim = 2;

// data type for our positions
typedef double Real;

// typedef for our particle layout type
typedef ParticleSpatialLayout<Real,Dim> playout_t;



int main(int argc, char *argv[])
{
  
  int i, j;
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  testmsg << "Particle test: Spatial layout: Begin." << endl;

  // create layout objects
  FieldLayout<Dim> FL(Index(100), Index(100), PARALLEL, PARALLEL, 8);

  // create an empty particle object, with 2D position and 2
  // integer attributes, using a uniform layout
  playout_t* PL = new playout_t(FL);
  GenArrayParticle<playout_t,int,2> P(PL);

  // initialize the particle object:
  // The basic method is this: on just one node, call create(N)
  // and then set the values for the N new particles.  Then, on
  // every node, call update to distribute the particles where
  // they are supposed to go.  This can be done several times
  // if the number of particles is very large and cannot all fit
  // on one processor.
  if (P.singleInitNode()) {
    testmsg << "Start of initialization ..." << endl;
    P.create(20);	// makes new particles at end of current array
    // put ten particles near right edge, and ten near left edge
    for (j=0; j < 10; j++) { // initialize new values
      P.R[j] = Vektor<Real,Dim>(5, j * 10);
      P.data[0][j] = Ippl::myNode();
      P.data[1][j] = 5;
    }
    for (j=10; j < 20; j++) { // initialize new values
      P.R[j] = Vektor<Real,Dim>(95, (j-10) * 10 + 2);
      P.data[0][j] = Ippl::myNode();
      P.data[1][j] = 95;
    }
  }
  testmsg << "Performing initial update ..." << endl;
  P.update();         // do update to move particles to proper loc
  testmsg << "End of initialization" << i << ": localnum=";
  testmsg << P.getLocalNum() << ", totalnum=";
  testmsg << P.getTotalNum() << endl;

  // establish connection for particles
  DataConnect *regcon = P.connect("pbconds");
  P.data[0].connect("data[0]", regcon);
  P.data[1].connect("data[1]", regcon);

  // hand off control to other agency for user interaction, if possible
  testmsg << "Doing interaction ..." << endl;
  P.interact();

  // Just print out the contents of node 0's data
  testmsg << "------------------" << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testmsg << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testmsg << P.R[i] << endl;
  }
  testmsg << "Updating connection ..." << endl;
  P.updateConnection();		// also does update for attributes

  // hand off control to other agency for user interaction, if possible
  testmsg << "Doing interaction ..." << endl;
  P.interact();

  //
  // test 1: 
  // add offset to X positions, moving particles beyond the border
  // bconds: none
  //

  testmsg << "Testing null boundary conditions ... moving by (10,0):" << endl;
  assign(P.R(0), P.R(0) + 10);
  P.update();
  testmsg << "------------------" << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testmsg << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testmsg << P.R[i] << endl;
  }
  testmsg << "Updating connection ..." << endl;
  P.updateConnection();		// also does update for attributes

  // hand off control to other agency for user interaction, if possible
  testmsg << "Doing interaction ..." << endl;
  P.interact();

  //
  // test 2: 
  // add offset to X & Y positions, moving particles beyond the border
  // bconds: periodic in X & Y
  //

  testmsg << "Testing periodic boundary conditions ... moving by (5,-1)";
  testmsg << endl;
  P.getBConds()[0] = ParticlePeriodicBCond;
  P.getBConds()[1] = ParticlePeriodicBCond;
  P.getBConds()[2] = ParticlePeriodicBCond;
  P.getBConds()[3] = ParticlePeriodicBCond;
  //  Vektor<Real,Dim> offset(5,-1);
  //  assign(P.R, P.R + offset);
  assign(P.R(0), P.R(0) + 5);
  assign(P.R(1), P.R(1) - 1);
  P.update();
  testmsg << "------------------" << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testmsg << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testmsg << P.R[i] << endl;
  }
  testmsg << "Updating connection ..." << endl;
  P.updateConnection();		// also does update for attributes

  // hand off control to other agency for user interaction, if possible
  testmsg << "Doing interaction ..." << endl;
  P.interact();

  //
  // test 3:
  // add offset to X & Y, moving particles beyond the border
  // bconds: reflective in X, periodic in Y
  //

  testmsg << "Testing reflective (X) and periodic (Y) boundary conditions ";
  testmsg << " ... moving by (-15, 10)" << endl;
  P.getBConds()[0] = ParticleReflectiveBCond;
  P.getBConds()[1] = ParticleReflectiveBCond;
  //  Vektor<Real,Dim> offset(5,-1);
  //  assign(P.R, P.R + offset);
  assign(P.R(0), P.R(0) - 15);
  assign(P.R(1), P.R(1) + 10);
  P.update();
  testmsg << "------------------" << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testmsg << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testmsg << P.R[i] << endl;
  }
  testmsg << "Updating connection ..." << endl;
  P.updateConnection();		// also does update for attributes

  // hand off control to other agency for user interaction, if possible
  testmsg << "Doing interaction ..." << endl;
  P.interact();

  //
  // test 4:
  // create new particles, and only update the numbers, do not do swapping
  //

  testmsg << "Testing flags - only do num update." << endl;
  P.create(5);
  P.R(0) = 145.0;
  P.R(1) = -40.0;
  P.setUpdateFlag(playout_t::SWAP, false);
  P.setUpdateFlag(playout_t::BCONDS, false);
  P.update();
  testmsg << "------------------" << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testmsg << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testmsg << P.R[i] << endl;
  }
  testmsg << "Now do complete update." << endl;
  P.setUpdateFlag(playout_t::ALL, true);
  P.update();
  testmsg << "------------------" << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testmsg << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testmsg << P.R[i] << endl;
  }

  testmsg << "Particle test: Spatial layout: End." << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: pbconds_vis.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: pbconds_vis.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/
