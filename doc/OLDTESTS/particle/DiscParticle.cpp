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
// lengths of axes
const int nx=100, ny=200;
// typedef for our particle layout type
typedef ParticleSpatialLayout< double, Dim, Cartesian<Dim> > playout_t;

int main(int argc, char *argv[])
{
  int i, j;
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  testmsg << "DiscParticle test: Spatial layout: Begin." << endl;

  // set mesh spacing
  double meshspacing = 0.5;
  double** delta = new double*[Dim];
  delta[0] = new double[nx];
  delta[1] = new double[ny];
  int indx;
  for (indx=0; indx<nx; indx++) delta[0][indx] = meshspacing;
  for (indx=0; indx<ny; indx++) delta[1][indx] = meshspacing;

  // set mesh origin
  Vektor<double,Dim> orig(-3.0, 5.0);

  // create the Mesh
  Index I(nx+1), J(ny+1);
  Cartesian<Dim> myMesh(I,J,delta,orig);

  for (int dimen=0; dimen<Dim; dimen++)
    delete [] delta[dimen];
  delete [] delta;

  // A FieldLayout which describes the spatial region we wish to move
  // the particles around in
  CenteredFieldLayout<Dim,Cartesian<Dim>,Cell> FL(myMesh, PARALLEL, PARALLEL);

  // create an empty particle object, with 2D position and 2
  // integer attributes
  playout_t* PL = new playout_t(FL,myMesh);
  playout_t* PL2 = new playout_t(FL,myMesh);
  GenArrayParticle<playout_t, int, 2> P(PL);
  GenArrayParticle<playout_t, int, 2> P2(PL2);

  // print out info about particles so far
  testmsg << endl << "Particle container contents = " << P << endl << endl;

  // initialize the particle object:
  // The basic method is this: on just one node, call create(N)
  // and then set the values for the N new particles.  Then, on
  // every node, call update to distribute the particles where
  // they are supposed to go.  This can be done several times
  // if the number of particles is very large and cannot all fit
  // on one processor.
  testmsg << "Start of initialization ..." << endl;
  if (P.singleInitNode()) {
    P.create(10);	      // makes new particles at end of current array
    for (j=0; j < 10; j++) {
      P.R[j] = myMesh.get_origin() + meshspacing *
	 Vektor<double,Dim>(double(j*(nx/20-1)),double(j*(ny/20-1)));
      P.data[0][j] = Ippl::myNode();
      P.data[1][j] = j;
    }
  }
  testmsg << "Performing initial update ..." << endl;
  P.update();         // do update to move particles to proper loc
  testmsg << "End of initialization: local num = ";
  testmsg << P.getLocalNum() << ", total num = " << P.getTotalNum() << endl;

  // print out info about particles so far
  testmsg << endl << "Particle container contents = " << P << endl << endl;

  // create an output file, via DataConnect connection
  testmsg << endl << "Creating output file 'outputdata' ..." << endl;
  DataConnect *dp = new FileDataConnect("outputdata.config");
  P.connect("outputdata", dp, DataSource::OUTPUT);
  P.updateConnection();

  // Just print out the contents of node 0's data
  testmsg << "------------------" << endl;
  testmsg << "Local Particles = " << P.getLocalNum() << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testmsg << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testmsg << P.R[i] << endl;
  }

  // test out bit-reverse coordinate initialization
  Vektor<double,Dim> upper(double(nx/2),double(ny/2));
  //  unsigned InitSeed = P.getTotalNum() * Ippl::myNode();
  //  Vektor<RNGBitReverse,Dim> BR(RNGBitReverse(2,InitSeed),
  //                               RNGBitReverse(3,InitSeed));
  //  assign(P.R(0), BR(0) * upper(0));
  //  assign(P.R(1), BR(1) * upper(1));
  P.R(0) = myMesh.get_origin()(0) + meshspacing * IpplRandom * upper(0);
  P.R(1) = myMesh.get_origin()(1) + meshspacing * IpplRandom * upper(1);
  P.update();
  P.updateConnection();

  // Just print out the contents of node 0's data
  testmsg << "------------------" << endl;
  testmsg << "Local Particles = " << P.getLocalNum() << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testmsg << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testmsg << P.R[i] << endl;
  }

  // test out random coordinate initialization
  Vektor<RNGSimple,Dim> RND(RNGSimple(Ippl::myNode()),
                            RNGSimple(Ippl::myNode()+Ippl::getNodes()));
  //  assign(P.R(0), RND(0) * upper(0));
  //  assign(P.R(1), RND(1) * upper(1));
  P.R(0) = myMesh.get_origin()(0) + meshspacing * RND(0) * upper(0);
  P.R(1) = myMesh.get_origin()(1) + meshspacing * RND(1) * upper(1);
  P.update();
  P.updateConnection();

  // Just print out the contents of node 0's data
  testmsg << "------------------" << endl;
  testmsg << "Local Particles = " << P.getLocalNum() << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testmsg << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testmsg << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testmsg << P.R[i] << endl;
  }

  // Read particles in from DiscField file
  testmsg << endl;
  testmsg << "Reading particles from DiscParticle file 'outputdata' ..."<<endl;
  DataConnect *dp2 = new FileDataConnect("outputdata.config");
  P2.connect("outputdata", dp2, DataSource::INPUT);
  for (int rec=0; rec < 3; ++rec) {
    testmsg << "------------------ for record " << rec << ":" << endl;
    P2.updateConnection();
    testmsg << "Local Particles read in = " << P2.getLocalNum() << endl;
    for (i=0; i < P2.getLocalNum(); i++) {
      testmsg << "LID=" << i << ", GID=" << P2.ID[i] << ", data[0]=";
      testmsg << P2.data[0][i] << ", data[1]=" << P2.data[1][i] << ", R=";
      testmsg << P2.R[i] << endl;
    }
  }

  // Create and write a set of particles in a single attrib
  testmsg << endl;
  testmsg << "Writing single-attribute DiscParticle ..." << endl;
  ParticleAttrib<double> attr;
  attr.create(5);
  attr.connect("outputattr", dp2, DataSource::OUTPUT);
  for (int a=0; a < 4; ++a) {
    attr = IpplRandom;
    attr.updateConnection();
  } 

  // Now read the attributes back in
  testmsg << endl;
  testmsg << "Reading attribute back in ..." << endl;
  attr.disconnect(dp2);
  attr.connect("outputattr", dp2, DataSource::INPUT);
  for (int b=0; b < 4; ++b) {
    attr.updateConnection();
    testmsg << "------------------ for record " << b << ":" << endl;
    testmsg << "Local Particles read in = " << attr.size() << endl;
    for (i=0; i < attr.size(); i++) {
      testmsg << "attr[" << i << "] = " << attr[i] << endl;
    }
  } 

  testmsg << endl << "------------------" << endl;
  testmsg << "Particle test: Spatial layout: End." << endl;

  delete dp;
  delete dp2;

  return 0;
}

/***************************************************************************
 * $RCSfile: DiscParticle.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:02:31 $
 * IPPL_VERSION_ID: $Id: DiscParticle.cpp,v 1.1.1.1 2000/11/30 21:02:31 adelmann Exp $ 
 ***************************************************************************/

