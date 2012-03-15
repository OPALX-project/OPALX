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

// number of particles to create
const unsigned Num = 10;

// data type for our positions
typedef double Real;

// typedef for our particle layout type
typedef ParticleSpatialLayout<Real,Dim> playout_t;


int main(int argc, char *argv[])
{
  int i, j;
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  testmsg << "Particle test: Particle expressions: Begin." << endl;

  // create layout objects
  FieldLayout<Dim> FL(Index(100), Index(200), PARALLEL, PARALLEL, 8);

  // create an empty particle object, with 2D position and 2
  // integer attributes, using a spatial layout
  playout_t* PL = new playout_t(FL);
  GenArrayParticle<playout_t, int, 3> P(PL);

  // initialize the particle object:
  // The basic method is this: on just one node, call create(N)
  // and then set the values for the N new particles.  Then, on
  // every node, call update to distribute the particles where
  // they are supposed to go.  This can be done several times
  // if the number of particles is very large an cannot all fit
  // on one processor.
  if (P.singleInitNode()) {
    testmsg << "Start of initialization ..." << endl;
    P.create(Num);	// makes new particles at end of current array
    for (j=0; j < Num; j++) { // initialize new values
      P.R[j] = Vektor<double,2>((double)(j*10), (double)(j*20));
      P.data[0][j] = Ippl::myNode() + 10;
      P.data[1][j] = j;
      P.data[2][j] = -1;
    }
  }
  testmsg << "Performing initial update ..." << endl;
  P.update();         // do update to move particles to proper loc
  testmsg << "Local number of particles = " << P.getLocalNum() << endl;
  testmsg << "Initial positions ..." << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "  Local particle " << i << ": R=" << P.R[i] << endl;
  }

  testmsg << "Testing sum(data[1] = 0...9) ..." << endl;
  testmsg << "  result = " << sum(P.data[1]) << endl;

  testmsg << "\nTesting multipply(10 + data[1] = 0...9) ..." << endl;
  testmsg << "  result = " << prod(P.data[1] + 1) << endl;

  testmsg << "\nTesting sum(data[1] + Random) ..." << endl;
  testmsg << "  result = " << sum(P.data[1] + IpplRandom) << endl;

  testmsg << "\nTesting R(0) = IpplRandom ..." << endl;
#ifdef __MWERKS__
  // ??? Really need to find if there's a way to make CW4 accept as op=()
  assign(P.R(0), IpplRandom * 100);
#else
  P.R(0) = IpplRandom * 100;
#endif // __MWERKS__
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "  Local particle " << i << ": R=" << P.R[i] << endl;
  }

  testmsg << "\nTesting data[2] = data[0] + data[1] ... " << endl;
  P.data[2] = P.data[0] + P.data[1];

  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "  Local particle " << i << ": data[2]=" << P.data[2][i];
    testmsg << ", data[0]=" << P.data[0][i];
    testmsg << ", data[1]=" << P.data[1][i] << endl;
  }

  testmsg << "\nTesting R(1) = R(0) ..." << endl;
#ifdef __MWERKS__
  // ??? Really need to find if there's a way to make CW4 accept as op=()
  assign(P.R(1), P.R(0));
#else
  P.R(1) = P.R(0);
#endif // __MWERKS__

  for (i=0; i < P.getLocalNum(); i++)
    testmsg << "  Local particle " << i << ": R=" << P.R[i] << endl;

  testmsg << " \nTesting data[2] = data[0] + 2 + 3 + log(4.0) ... " << endl;
  P.data[2] = P.data[0] + 2.0 + 3.0 + log(4.0);
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "  Local particle " << i << ": data[2]=" << P.data[2][i];
    testmsg << ", data[0]=" << P.data[0][i] << endl;
  }

  testmsg << " \nTesting R(0) = log(R(1) + 100.0) + data[1] ... " << endl;
#ifdef __MWERKS__
  // ??? Really need to find if there's a way to make CW4 accept as op=()
  assign(P.R(0), log(P.R(1) + 100.0) + P.data[1]);
#else
  P.R(0) = log(P.R(1) + 100.0) + P.data[1];
#endif // __MWERKS__
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "  Local particle " << i << ": R(0)=" << P.R[i](0);
    testmsg << ", R(1)=" << P.R[i](1);
    testmsg << ", data[1]=" << P.data[1][i] << endl;
  }

  testmsg << "\nTesting R(0) = fmod(R(0) + 10.0, 10.0) ..." << endl;

#ifdef __MWERKS__
  // ??? Really need to find if there's a way to make CW4 accept as op=()
  assign(P.R(0), fmod(P.R(0) + 10.0, 10.0));
#else
  P.R(0) = fmod(P.R(0) + 10.0, 10.0);
#endif // __MWERKS__
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "  Local particle " << i << ": R = " << P.R[i] << endl;
  }

  testmsg << " \nTesting data[2] = -50.0 ... " << endl;
  P.data[2] = -50;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "  Local particle " << i << ": data[2]=" << P.data[2][i] <<endl;
  }

  testmsg << " \nTesting data[2] = where(data[1] > 5, data[2], R(0)) ...";
  testmsg << endl;
  P.data[2] = where( gt(P.data[1],5), P.data[2], P.R(0) );
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "  Local particle " << i << ": data[2]=" << P.data[2][i];
    testmsg << ", R(0)=" << P.R[i](0);
    testmsg << ", data[1]=" << P.data[1][i] << endl;
  }

  testmsg << " \nTesting data[2] = where(data[1] <= 5, 1000, R(1)) ...";
  testmsg << endl;
  P.data[2] = where( le(P.data[1],5), 1000, P.R(1));
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "  Local particle " << i << ": data[2]=" << P.data[2][i];
    testmsg << ", R(1)=" << P.R[i](1) << endl;
  }

  testmsg << "\nParticle test: Particle expressions: End." << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: parexp.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:02:43 $
 * IPPL_VERSION_ID: $Id: parexp.cpp,v 1.1.1.1 2000/11/30 21:02:43 adelmann Exp $ 
 ***************************************************************************/
