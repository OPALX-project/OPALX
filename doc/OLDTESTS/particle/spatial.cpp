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
#include "Utility/RNGBitReverse.h"

// dimension of our positions
const unsigned Dim = 2;

// lengths of axes
const int nx=100, ny=200;

// number of particles
const int np=20;

// typedef for our particle layout type
typedef ParticleSpatialLayout< double, Dim, Cartesian<Dim> > playout_t;


int main(int argc, char *argv[])
{
  TAU_PROFILE("main", "int (int, char**)", TAU_DEFAULT);

  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);
  Inform testinfo(argv[0],INFORM_ALL_NODES);

  testmsg << "Particle test: Spatial layout: Begin." << endl;

  int i, j, nploc;
  unsigned vnodes;

  if (Ippl::getNodes() > 1) {
    // read in number of vnodes from node 0, broadcast to others
    Message* mess;
    int parent = 0;
    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG0);
    if (Ippl::myNode() == parent) {
      cout << "Enter number of vnodes (0 to exit): ";
      cin >> vnodes;
      mess = new Message;
      mess->put(vnodes);
      Ippl::Comm->broadcast_others(mess, tag);
    }
    else {
      mess = Ippl::Comm->receive_block(parent, tag);
      mess->get(vnodes);
      delete mess;
    }
  }
  else {
    cout << "Enter number of vnodes (0 to exit): ";
    cin >> vnodes;
  }

  if (vnodes == 0) {
    testmsg << "Number of vnodes = 0, exiting ... " << endl;
    testmsg << "Particle test: Spatial layout: End." << endl;
    return 0;
  }
    
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
  CenteredFieldLayout<Dim,Cartesian<Dim>,Cell>
    FL(myMesh, PARALLEL, PARALLEL, vnodes);

  // A ParticleLayout constructed from the FieldLayout and Mesh
  playout_t* PL = new playout_t(FL,myMesh);
  // create an empty particle object, with 2D position and 2
  // integer attributes
  GenArrayParticle<playout_t, int, 2> P(PL);

  // turn off update flag for optimized destroy method
  // P.setUpdateFlag(playout_t::OPTDESTROY,false);

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

  nploc = np / Ippl::getNodes();
  P.create(nploc);	      // makes new particles at end of current array
  for (j=0; j < nploc; j++) {
    P.R[j] = myMesh.get_origin() + meshspacing *
      Vektor<double,Dim>(static_cast<double>(j+nploc*Ippl::myNode())/np*nx,
                         static_cast<double>(j+nploc*Ippl::myNode())/np*ny);
    P.data[0][j] = Ippl::myNode();
    P.data[1][j] = j;
  }
  testinfo << "Local Particles = " << P.getLocalNum() << endl;

  testmsg << "Performing initial update ..." << endl;
  P.update();         // do update to move particles to proper loc
  testmsg << "End of initialization: local num = ";
  testmsg << P.getLocalNum() << ", total num = " << P.getTotalNum() << endl;

  // print out info about particles so far
  testmsg << endl << "Particle container contents = " << P << endl << endl;

  // print out the contents of each node's data
  //  testmsg << "------------------" << endl;
  testinfo << "Local Particles = " << P.getLocalNum() << endl;

  for (i=0; i < P.getLocalNum(); i++) {
    testinfo << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testinfo << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testinfo << P.R[i] << endl;
  }

  // test out bit-reverse coordinate initialization
  testmsg << "Testing bit-reverse coordinate initialization ..." << endl;
  Vektor<double,Dim> upper(static_cast<double>(nx),static_cast<double>(ny));
  unsigned InitSeed = P.getTotalNum() * Ippl::myNode();
  Vektor<RNGBitReverseSequence,Dim> BR(RNGBitReverseSequence(2,InitSeed),
                                       RNGBitReverseSequence(3,InitSeed));
  //  assign(P.R(0), BR(0) * upper(0));
  //  assign(P.R(1), BR(1) * upper(1));
#ifdef __MWERKS__
  // ??? Really need to find if there's a way to make CW4 accept as op=()
  assign(P.R(0), myMesh.get_origin()(0) + meshspacing * BR(0) * upper(0));
  assign(P.R(1), myMesh.get_origin()(1) + meshspacing * BR(1) * upper(1));
#else
  P.R(0) = myMesh.get_origin()(0) + meshspacing * BR(0) * upper(0);
  P.R(1) = myMesh.get_origin()(1) + meshspacing * BR(1) * upper(1);
#endif // __MWERKS__
  //  P.R = myMesh.get_origin() + meshspacing * BR * upper;
  P.update();

  // print out the contents of each node's data
  testmsg << "------------------" << endl;
  testinfo << "Local Particles = " << P.getLocalNum() << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testinfo << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testinfo << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testinfo << P.R[i] << endl;
  }

  // test out random coordinate initialization
  testmsg << "Testing random coordinate initialization ..." << endl;
  Vektor<RNGSimpleSequence,Dim> RND(RNGSimpleSequence(Ippl::myNode()),
                                    RNGSimpleSequence(Ippl::myNode()+
                                                      Ippl::getNodes()));
  //  assign(P.R(0), RND(0) * upper(0));
  //  assign(P.R(1), RND(1) * upper(1));
#ifdef __MWERKS__
  // ??? Really need to find if there's a way to make CW4 accept as op=()
  assign(P.R(0), myMesh.get_origin()(0) + meshspacing * RND(0) * upper(0));
  assign(P.R(1), myMesh.get_origin()(1) + meshspacing * RND(1) * upper(1));
#else
  P.R(0) = myMesh.get_origin()(0) + meshspacing * RND(0) * upper(0);
  P.R(1) = myMesh.get_origin()(1) + meshspacing * RND(1) * upper(1);
#endif // __MWERKS__
  //  P.R = myMesh.get_origin() + meshspacing * RND * upper;
  P.update();

  // print out the contents of each node's data
  testmsg << "------------------" << endl;
  testinfo << "Local Particles = " << P.getLocalNum() << endl;
  for (i=0; i < P.getLocalNum(); i++) {
    testinfo << "LID=" << i << ", GID=" << P.ID[i] << ", data[0]=";
    testinfo << P.data[0][i] << ", data[1]=" << P.data[1][i] << ", R=";
    testinfo << P.R[i] << endl;
  }

  testmsg << "Particle test: Spatial layout: End." << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: spatial.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:03:46 $
 * IPPL_VERSION_ID: $Id: spatial.cpp,v 1.1.1.1 2000/11/30 21:03:46 adelmann Exp $ 
 ***************************************************************************/
