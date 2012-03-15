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


// dimension of our positions
const unsigned Dim = 2;

// typedef our particle layout type
typedef ParticleSpatialLayout<double,Dim> playout_t;

// other constants
const int nx=20, ny=20;              // size of domain is nx X ny
const unsigned int totalP = 40;      // number of particles to create


class ChargedParticles {
public:
  ParticleAttrib<Vektor<double,Dim> > R;
  ParticleAttrib<unsigned>            ID;
  ParticleAttrib<double>              qm;
};


int main(int argc, char *argv[])
{
  // initialize IPPL
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);

  msg << "Paws Particle recv test code." << endl;
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
  msg << "Creating particle attributes ..." << endl;
  ChargedParticles P;

  // connect up the Field and Particles via Paws
  msg << "Registering Paws data objects ..." << endl;
  DataConnect *dc = new PawsDataConnect("particle-recv");
  dc->connect("EField", EFDMag, DataSource::INPUT);
  dc->connect("R", P.R, DataSource::INPUT);
  dc->connect("ID", P.ID, DataSource::INPUT);
  dc->connect("QM", P.qm, DataSource::INPUT);

  // wait for ready
  msg << "Waiting for ready on receive end ..." << endl;
  dc->ready();

  // put in a couple particles first
  P.R.create(3);
  P.R = IpplRandom + 10.0;

  // receive the initial values of the paws objects
  msg << "Updating PAWS connection for the first time ..." << endl;
  dc->updateConnections();

  // go through a loop, receiving new particles, and printing them out
  msg << "\nStarting iteration loop:" << endl;
  msg <<   "------------------------" << endl;
  for (int i=0; i < totalP/2; ++i) {
    msg << "On iteration " << i << ":" << endl;

    P.R.create(3);
    P.R = IpplRandom + 10.0;
    msg << "  Before update, positions are:" << endl;
    msg << "   R[0] = " << P.R[0] << endl;
    msg << "   R[1] = " << P.R[1] << endl;
    msg << "   R[2] = " << P.R[2] << endl;

    // Update the PAWS connection
    msg << "  Updating PAWS connection ..." << endl;
    dc->updateConnections();

    // Print out first ten particle positions
    int maxn = 10;
    if (P.R.size() < 10)
      maxn = P.R.size();
    msg << "  Positions and ID's for first " << maxn << " particles:" << endl;
    for (int j=0; j < maxn; ++j) {
      msg << "   R[" << j << "] = " << P.R[j];
      msg << ", ID[" << j << "] = " << P.ID[j] << endl;
    }
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

