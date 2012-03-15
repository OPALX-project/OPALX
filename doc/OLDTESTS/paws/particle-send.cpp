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

#include "Ippl.h"


template<class PL>
class ChargedParticles : public ParticleBase<PL> {
public:
  ParticleAttrib<double>     qm;       // charge-to-mass ratio
  typename PL::ParticlePos_t V;        // particle velocity

  ChargedParticles(PL* pl) : ParticleBase<PL>(pl) {
    // register the particle attributes
    addAttribute(qm);
    addAttribute(V);
  }
};


// dimension of our positions
const unsigned Dim = 2;

// typedef our particle layout type
typedef ParticleSpatialLayout<double,Dim> playout_t;

// other constants
const int nx=20, ny=20;              // size of domain is nx X ny
const unsigned int totalP = 40;      // number of particles to create


int main(int argc, char *argv[])
{
  // initialize IPPL
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);

  msg << "Paws Particle send test code." << endl;
  msg << "-----------------------------" << endl;

  // create layout objects
  Index I(nx), J(ny);
  Index I1(nx+1), J1(ny+1);
  UniformCartesian<Dim> mymesh(I1,J1);
  CenteredFieldLayout<Dim,UniformCartesian<Dim>,Cell> FL(mymesh);

  // create a test Field object
  msg << "Creating field object ..." << endl;
  Field<double,Dim> EFDMag(mymesh, FL);
  EFDMag[I][J] = 0.001 * I * J;

  // create an empty ChargedParticles object, setting it to use periodic BC's
  msg << "Creating particle object ..." << endl;
  playout_t* PL = new playout_t(FL);
  ChargedParticles<playout_t> P(PL);
  P.getBConds()[0] = ParticlePeriodicBCond;
  P.getBConds()[1] = ParticlePeriodicBCond;
  P.getBConds()[2] = ParticlePeriodicBCond;
  P.getBConds()[3] = ParticlePeriodicBCond;

  // connect up the Field and Particles via Paws
  msg << "Registering Paws data objects ..." << endl;
  DataConnect *dc = new PawsDataConnect("particle-send");
  dc->connect("EField", EFDMag, DataSource::OUTPUT);
  dc->connect("R", P.R, DataSource::OUTPUT);
  dc->connect("ID", P.ID, DataSource::OUTPUT);
  dc->connect("QM", P.qm, DataSource::OUTPUT);

  // wait for ready
  msg << "Waiting for ready on sending end ..." << endl;
  dc->ready();

  // create some initial particles, and assign them random positions and
  // random charge-to-mass ratios
  msg << "Initializing particles ..." << endl;
  P.globalCreate(totalP);
  P.qm = (10.0 * IpplRandom);
  assign(P.R(0), nx * IpplRandom);
  assign(P.R(1), ny * IpplRandom);
  P.update();

  // send the initial values of the paws objects
  msg << "Updating PAWS connection for the first time ..." << endl;
  dc->updateConnections();

  // go through a loop, deleting some particles, and updating the
  // PAWS connection
  msg << "\nStarting iteration loop:" << endl;
  msg <<   "------------------------" << endl;
  for (int i=0; i < totalP/2; ++i) {
    msg << "On iteration " << i << ":" << endl;

    // Delete a particle from one node's list, and update the particles
    msg << "  Deleting particle 0 from node " << (i % Ippl::getNodes());
    msg << " ..." << endl;
    if (P.R.size() > 0 && Ippl::myNode() == (i % Ippl::getNodes()))
      P.destroy(1, 0);
    P.update();

    // Update the Field data as well
    msg << "  Multipplying EFDMag by 10 ..." << endl;
    EFDMag *= 10.0;

    // Print out first ten particle positions
    int maxn = 10;
    if (P.R.size() < 10)
      maxn = P.R.size();
    msg << "  Positions and ID's for first " << maxn << " particles:" << endl;
    for (int j=0; j < maxn; ++j) {
      msg << "  R[" << j << "] = " << P.R[j];
      msg << ", ID[" << j << "] = " << P.ID[j] << endl;
    }

    // Update the PAWS connection
    msg << "  Updating PAWS connection ..." << endl;
    dc->updateConnections();
  }

  msg << "\nFinishing up:" << endl;
  msg <<   "-------------" << endl;

  // delete everything and shut down
  msg << "Deleting DataConnect ..." << endl;
  delete dc;

  msg << "Returning from main()." << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

