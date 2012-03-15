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

/* This test program sets up a simple sine-wave electric field in 2D, 
   creates a population of particles with random q/m values (charge-to-mass
   ratio) and velocities, and then tracks their motions in the static
   electric field using nearest-grid-point interpolation. */

#include "Ippl.h"
#include "Particle/IntNGP.h"
#include "Utility/RNGBitReverse.h"


template<class PL>
class ChargedParticles : public ParticleBase<PL> {
public:
  ParticleAttrib<double>     qm;       // charge-to-mass ratio
  typename PL::ParticlePos_t V;        // particle velocity
  typename PL::ParticlePos_t E;        // electric field at particle position
  ChargedParticles(PL* pl) : ParticleBase<PL>(pl) {
    // register the particle attributes
    addAttribute(qm);
    addAttribute(V);
    addAttribute(E);
  }
};

// dimension of our positions
const unsigned Dim = 2;

// typedef our particle layout type
typedef ParticleSpatialLayout<double,Dim> playout_t;


const int nx=200, ny=200;              // size of domain is nx X ny
const unsigned int totalP = 400;      // number of particles to create
const int nt = 20;                    // total number of timesteps

const double pi = acos(-1.0);
const double qmmax = 1.0;       // maximum value for particle q/m
const double dt = 1.0;          // size of timestep

int main(int argc, char *argv[]){
  TAU_PROFILE("main()", "int (int, char **)", TAU_DEFAULT);
  unsigned int len;
  int i, it;

  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);
  testmsg << "Particle test PIC2d: Begin." << endl;

  // potential phi = phi0 * sin(2*pi*x/Lx) * cos(4*pi*y/Ly)
  double phi0 = 0.01*nx;             // electric potential amplitude

  // create interpolater object (nearest-grid-point method)
  IntNGP myinterp;

  // create layout objects
  Index I(nx), J(ny);
  Index I1(nx+1), J1(ny+1);
  UniformCartesian<Dim> mymesh(I1,J1);
  CenteredFieldLayout<Dim,UniformCartesian<Dim>,Cell> FL(mymesh);

  // initialize static electric field
  Field<Vektor<double,Dim>,Dim> EFD(mymesh,FL);
  Field<double,Dim> EFDMag(mymesh,FL);
  assign(EFD[I][J](0), -2.0*pi*phi0/nx * cos(2.0*pi*(I+0.5)/nx) * cos(4.0*pi*(J+0.5)/ny));
  assign(EFD[I][J](1),  4.0*pi*phi0/ny * sin(2.0*pi*(I+0.5)/nx) * sin(4.0*pi*(J+0.5)/ny));
  EFDMag = dot(EFD, EFD);

  // create an empty ChargedParticles object, setting it to use periodic BC's
  playout_t* PL = new playout_t(FL);
  ChargedParticles<playout_t> P(PL);
  P.getBConds()[0] = ParticlePeriodicBCond;
  P.getBConds()[1] = ParticlePeriodicBCond;
  P.getBConds()[2] = ParticlePeriodicBCond;
  P.getBConds()[3] = ParticlePeriodicBCond;
  DataConnect *dc = DataConnectCreator::create(argv[0]);
  dc->connect("P", P);
  dc->connect("P.qm", P.qm);
  dc->connect("P.velocity", P.V);
  dc->connect("P.E-Field", P.E);
  dc->connect("S E-Field", EFDMag);
  dc->connect("V E-Field", EFD);

  // initialize the particle object: initialize equal number of particles 
  // on each node
  int localP = totalP / Ippl::getNodes();
  P.create(localP);

  // quiet start for particle positions

  // Vektor<unsigned int,Dim> base(2,3); // bases for bit-reversal
  // Vektor<double,Dim> upper(nx,ny); // range for particle positions
  // Assign_bit_reverse(P.R, P.ID, base, upper);
  //assign(P.R(1), ny * RNGBitReverse(3, localP*Ippl::myNode()+1));
  P.R(0) = nx * IpplRandom;
  P.R(1) = ny * IpplRandom;

  // random initialization for charge-to-mass ratio
  P.qm = (2 * qmmax * IpplRandom) - qmmax;

  // redistribute particles based on spatial layout
  P.update(); 

  testmsg << "PIC2d setup complete." << endl;
  testmsg << "------------------" << endl;

  testmsg << "Initial particle data:" << endl;
  Ippl::Comm->barrier();  // synch before I/O
  len = P.getLocalNum();
  for (i=0; i < len; i++) {
    testmsg << "LID=" << setw(3) << i << ", GID=" << setw(3) << P.ID[i];
    testmsg << ", R=" << P.R[i] << ", V=" << P.V[i];
    testmsg << ", q/m=" << P.qm[i] << endl;
  }
  Ippl::Comm->barrier();  // wait for I/O to complete

  dc->updateConnections();
  EFDMag.interact();

  // begin main timestep loop
  for (it=0; it<nt; it++) {
    testmsg << "Time = " << it*dt << endl;
    //    testmsg << "Current particles = " << P << endl;

    /*
    if (it % pfreq == 0) {
      // time to do a data printout
      Ippl::Comm->barrier();  // synch before I/O
      len = P.getLocalNum();
      for (int i=0; i < len; i++) {
        testmsg << "GID=" << P.ID[i] << ", R=" << P.R[i] << ", V=";
        testmsg << P.V[i] << endl;
      }
      Ippl::Comm->barrier();  // wait for I/O to complete
    }
    */

    // advance the particle positions
    // basic leapfrogging timestep scheme.  velocities are offset
    // by half a timestep from the positions.
    testmsg << "Advance particle positions..." << endl;
    P.R = P.R + dt * P.V;

    // update particle distribution across processors
    testmsg << "Particle update..." << endl;
    // BinaryRepartition(P);
    P.update();
    //    testmsg << "New particles = " << P << endl;

    testmsg << "Gather E field..." << endl;
    // gather the local value of the E field
    P.E.gather(EFD, P.R, myinterp);

    // advance the particle velocities
    testmsg << "Advance particle velocities..." << endl;
    P.V = P.V + dt * P.qm * P.E;

    dc->updateConnections();
    EFDMag.interact();
  }

  // final data output
  testmsg << "PIC2d timestep loop complete!" << endl;
  testmsg << "------------------" << endl;
  testmsg << "Final particle data:" << endl;
  Ippl::Comm->barrier();  // synchronize before final I/O
  len = P.getLocalNum();
  for (i=0; i < len; i++) {
    testmsg << "LID=" << setw(3) << i << ", GID=" << setw(3) << P.ID[i];
    testmsg << ", R=" << P.R[i] << ", V=" << P.V[i] << endl;
  }
  Ippl::Comm->barrier();  // wait for I/O to complete

  testmsg << "Particle test PIC2d: End." << endl;

  delete dc;

  return 0;
}

/***************************************************************************
 * $RCSfile: PIC2d.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:02:35 $
 * IPPL_VERSION_ID: $Id: PIC2d.cpp,v 1.1.1.1 2000/11/30 21:02:35 adelmann Exp $ 
 ***************************************************************************/
