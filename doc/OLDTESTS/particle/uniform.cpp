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

#include "Ippl.h"


// dimension of our positions
const unsigned Dim = 2;

// define our particle layout type
typedef ParticleUniformLayout<double,Dim> playout_t;


int main(int argc, char *argv[])
{
  int i, j;
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  testmsg << "Particle test: Uniform layout: Begin." << endl;

  // create an empty particle object, with 2D position and 2
  // integer attributes, using a default uniform layout
  GenArrayParticle< playout_t, int, 2> P(new playout_t());

  // initialize the particle object:
  // The basic method is this: on just one node, call create(N)
  // and then set the values for the N new particles.  Then, on
  // every node, call update to distribute the particles where
  // they are supposed to go.
  testmsg << "Start of initialization ..." << endl;
  if (P.singleInitNode()) {
    P.create(20);		// makes new particles at end of current array
    for (j=0; j < 20; j++) {
      P.R[j] = Vektor<double,2>((double)(Ippl::myNode()), (double)j);
      P.data[0][j] = Ippl::myNode();
      P.data[1][j] = j;
    }
  }
  testmsg << "Performing initial update ..." << endl;
  P.update();   		// do update to move particles to proper loc
  testmsg << "End of initialization: local num = ";
  testmsg << P.getLocalNum() << ", total num = " << P.getTotalNum() << endl;

  // print out info about particles so far
  testmsg << endl << "Particle container contents = " << P << endl << endl;

  // Just print out the contents of node 0's data
  testmsg << "------------------" << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testmsg << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testmsg << P.R[i] << endl;
  }

  testmsg << "Particle test: Uniform layout: End." << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: uniform.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:03:48 $
 * IPPL_VERSION_ID: $Id: uniform.cpp,v 1.1.1.1 2000/11/30 21:03:48 adelmann Exp $ 
 ***************************************************************************/
