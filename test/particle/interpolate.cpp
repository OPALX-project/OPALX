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
 ***************************************************************************/

#include "Ippl.h"
#include "Utility/RNGBitReverse.h"

// SpecialParticle class definition
template <class PLayout, class T, unsigned Dim>
class SpecialParticle : public GenArrayParticle<PLayout,T,Dim> {

    typedef typename PLayout::Position_t      Position_t;
    typedef typename InterpolatorTraits<Position_t,Dim,IntNGP>::Cache_t NGP_t;
    typedef typename InterpolatorTraits<Position_t,Dim,IntCIC>::Cache_t CIC_t;
    typedef typename InterpolatorTraits<Position_t,Dim,IntSUDS>::Cache_t SUDS_t;
    
public:
    // additional attributes for this class
    ParticleAttrib< NGP_t > CacheNGP;
    ParticleAttrib< CIC_t > CacheCIC;
    ParticleAttrib< SUDS_t > CacheSUDS;
    
    SpecialParticle(PLayout* L) : GenArrayParticle<PLayout,T,Dim>(L) {
      addAttribute(CacheNGP);
      addAttribute(CacheCIC);
      addAttribute(CacheSUDS);
    }

private:
    SpecialParticle();
};

// dimension of our positions
const unsigned Dim = 2;

// typedef our particle layout type
typedef ParticleSpatialLayout<double,Dim> playout_t;

// sizes of domain
const int nx=8, ny=8;


int main(int argc, char *argv[])
{
// Macro to invoke profiling of MC++.
  
// This profiling macro and all others are a no-op if profiling is off.

// Create some tau profile timers
  
  
  
  Ippl ippl(argc, argv);
  Inform testinfo(argv[0]);
  Inform testmsg(argv[0], INFORM_ALL_NODES);

  testinfo << "Particle test: Spatial layout: Begin." << endl;

  unsigned int vnodes;
  int tag=9999, parent=0;
  if (Ippl::myNode() == parent) {
    cout << "Enter number of vnodes: ";
    cin >> vnodes;
    Message* mess = new Message;
    putMessage(*mess,vnodes);
    Ippl::Comm->broadcast_others(mess,tag);
  }
  else {
    Message* rmsg = Ippl::Comm->receive_block(parent,tag);
    getMessage(*rmsg,vnodes);
    delete rmsg;
  }
  // create layout objects
  Index I(nx), J(ny);
  Index I1(nx+1), J1(ny+1);
  UniformCartesian<Dim> mymesh(I1,J1);
  CenteredFieldLayout<Dim,UniformCartesian<Dim>,Cell>
    FL(mymesh,PARALLEL,PARALLEL,vnodes);

  // create an empty particle object, with 2D position and 2
  // double attributes, using a spatial layout
  playout_t* PL = new playout_t(FL);
  SpecialParticle< playout_t, double, 2> P(PL);

  // initialize the particle object:
  // The basic method is this: on just one node, call create(N)
  // and then set the values for the N new particles.  Then, on
  // every node, call update to distribute the particles where
  // they are supposed to go.  This can be done several times
  // if the number of particles is very large an cannot all fit
  // on one processor.
  testinfo << "Start of initialization ..." << endl;
  int i, len;
  if (Ippl::myNode() == parent) { // only on Node 0
    P.create(10);	// makes new particles at end of current array
    Vektor<double,Dim> upper(nx,ny);
#ifdef __MWERKS__
    assign(P.R(0), IpplRandom * upper(0));
    assign(P.R(1), IpplRandom * upper(1));
#else
    P.R(0) = IpplRandom * upper(0);
    P.R(1) = IpplRandom * upper(1);
#endif // __MWERKS__
    // initialize the particle data fields
    P.data[0] = 0;
    RNGBitReverseSequence BR(2,1);
    P.data[1] = BR;
  }
  testinfo << "Performing initial update ..." << endl;
  P.update();         // do update to move particles to proper loc
  len = P.getLocalNum();
  testinfo << "End of particle initialization: len = " << len << endl;

  // initialize some Fields for interpolation
  GuardCellSizes<Dim> GC(1);
  BConds<double,Dim,UniformCartesian<Dim>,Cell> BC;
  BC[0] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Cell>(0);
  BC[1] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Cell>(1);
  BC[2] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Cell>(2);
  BC[3] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Cell>(3);
  Field<double,Dim,UniformCartesian<Dim>,Cell>
    A(mymesh,FL,GC,BC), B(mymesh,FL,GC,BC);
  A[I][J] = I*J + 1.0;
  B = 0.0;

  testinfo << "End of data initialization " << endl;
  // Just print out the contents of particle data
  testinfo << "------------------" << endl;
  for (i=0; i < len; i++) {
    testmsg << "LID = " << i << ", GID = " << P.ID[i] << ", R = ";
    testmsg << P.R[i] << ", data[0] = " << P.data[0][i] << ", data[1] = ";
    testmsg << P.data[1][i] << endl;
  }
  // print out the fields
  testinfo << "------------------" << endl;
  testinfo << "Field A:" << endl;
  A.write(cout);
  cout << endl;
  testinfo << "------------------" << endl;
  testinfo << "Field B:" << endl;
  B.write(cout);
  cout << endl;

  // NGP test: gather Field A into P.data[0] and scatter P.data[1] into B.
  testinfo << "Performing NGP gather of A into P.data[0] ... " << endl;
  
  // perform gather and cache mesh info in index attribute
  P.data[0].gather(A, P.R, IntNGP(), P.CacheNGP);
  
  testinfo << "Performing NGP scatter of P.data[1] onto B ... " << endl;
  
  // perform scatter, reusing cached mesh info
  P.data[1].scatter(B, IntNGP(), P.CacheNGP);
  

  // Just print out the contents of particle data
  testinfo << "------------------" << endl;
  for (i=0; i < len; i++) {
    testmsg << "R = " << P.R[i] << ", data[0] = " << P.data[0][i] << endl;
  }

  // print out the fields
  testinfo << "------------------" << endl;
  testinfo << "Field B:" << endl;
  B.write(cout);
  cout << endl;

  // Perform 2nd NGP scatter, to test handling of guard cells
  testinfo << "Performing 2nd NGP scatter of P.data[1] onto B ... " << endl;
  // B = 0.0;
  
  // perform scatter, reusing cached mesh info
  P.data[1].scatter(B, IntNGP(), P.CacheNGP);
  

  // Print out the fields
  testinfo << "------------------" << endl;
  testinfo << "Field B:" << endl;
  B.write(cout);
  cout << endl;

  // CIC test: gather Field A into P.data[0] and scatter P.data[1] into B.
  testinfo << "Performing CIC gather of A into P.data[0] ... " << endl;
  
  // perform gather and cache mesh info in index and delta attributes
  P.data[0].gather(A, P.R, IntCIC(), P.CacheCIC);
  
  testinfo << "Performing CIC scatter of P.data[1] onto B ... " << endl;
  B = 0.0;
  
  // perform scatter, reusing cached mesh info
  P.data[1].scatter(B, IntCIC(), P.CacheCIC);
  

  // Just print out the contents of node 0's data
  testinfo << "------------------" << endl;
  for (i=0; i < len; i++) {
    testmsg << "R = " << P.R[i] << ", data[0] = " << P.data[0][i] << endl;
  }

  // print out the fields
  testinfo << "------------------" << endl;
  testinfo << "Field B:" << endl;
  B.write(cout);
  cout << endl;

  // do a second CIC scatter, to test handling of guard cells
  testinfo << "Performing 2nd CIC scatter of P.data[1] onto B ... " << endl;
  // B = 0.0;
  
  // perform scatter, reusing cached mesh info
  P.data[1].scatter(B, IntCIC(), P.CacheCIC);
  

  // print out the fields
  testinfo << "------------------" << endl;
  testinfo << "Field B:" << endl;
  B.write(cout);
  cout << endl;

  // SUDS test: gather Field A into P.data[0] and scatter P.data[1] into B.
  testinfo << "Performing SUDS gather of A into P.data[0] ... " << endl;
  
  // perform gather and cache mesh info in index and delta attributes
  P.data[0].gather(A, P.R, IntSUDS(), P.CacheSUDS);
  
  testinfo << "Performing SUDS scatter of P.data[1] onto B ... " << endl;
  B = 0.0;
  
  // perform scatter, reusing cached mesh info
  P.data[1].scatter(B, IntSUDS(), P.CacheSUDS);
  

  // Just print out the contents of node 0's data
  testinfo << "------------------" << endl;
  for (i=0; i < len; i++) {
    testmsg << "R = " << P.R[i] << ", data[0] = " << P.data[0][i] << endl;
  }

  // print out the fields
  testinfo << "------------------" << endl;
  testinfo << "Field B:" << endl;
  B.write(cout);
  cout << endl;

  // Perform 2nd SUDS scatter, to test handling of guard cells
  testinfo << "Performing 2nd SUDS scatter of P.data[1] onto B ... " << endl;
  // B = 0.0;
  
  // perform scatter, reusing cached mesh info
  P.data[1].scatter(B, IntSUDS(), P.CacheSUDS);
  

  // print out the fields
  testinfo << "------------------" << endl;
  testinfo << "Field B:" << endl;
  B.write(cout);
  cout << endl;

  Ippl::Comm->barrier();
  testinfo << "Particle test: Interpolation: End." << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: interpolate.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: interpolate.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/
