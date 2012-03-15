//////////////////////////////////////////////////////////////////////////////
//
// Test program for initial exchange of particles to simulate what happens
// when the initializer supplies arrays of location, velocity and tags.
// Because of particle movement, a processor might have particles that
// really should be alive on another processor and we want to place them
// there before calling ParticleExchange.
//
// This simulates by calling ParticleDistribute to read files which have
// been written for ONE_TO_ONE.  Next we call a test class TestExchange which
// takes a veneer of the good alive particles on a processor and sends them
// to neighbors, and deletes them on this processor.  Now when 
// InitialExchange is called, the goal is to have the same particles that
// were present immediately after ParticleDistribute was called.  So we
// send away some particles and then get them back in the test.
//
// Added a further test which writes all particles in ASCII after the
// reading from files, dist.[0-N], and again after the InitialExchange,
// init.[0-N].  Each processor writes its own files so if you do
//
//    % cat dist* | sort -n > x.dist
//    % cat init* | sort -n > x.init
//    % diff x.dist x.init
//
// there should be no difference if all data transferred together.
//
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>

#include "Partition.h"
#include "ParticleDistribute.h"
#include "TestExchange.h"
#include "InitialExchange.h"

#include <mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc != 6) {
    cout << "Usage: mpirun -np # InitialTest in rL d ";
    cout << "[RECORD|BLOCK] ONE_TO_ONE" << endl;
  }

  // Base file name (Actual file names are basename.proc#
  int i = 1;
  string inFile = argv[i++];

  // Physical coordinate box size
  POSVEL_T rL = atof(argv[i++]);

  // Physical coordinate dead zone area
  POSVEL_T deadSize = atof(argv[i++]);

  // Input is BLOCK or RECORD structured data
  string dataType = argv[i++];

  // Distribution mechanism is ROUND_ROBIN which selects alive and dead
  // as the data files are read and passed around
  // or ONE_TO_ONE where the alive data is immediately read on a processor
  // and the dead particles must be bundled and shared
  string distributeType = argv[i++];

  // Initialize the partitioner which uses MPI Cartesian Topology
  //Partition::initialize(argc, argv);
  MPI_Init(&argc, &argv);
  Partition::initialize();

  // Construct the particle distributor, exchanger and halo finder
  ParticleDistribute distribute;
  TestExchange test;
  InitialExchange exchange;

  // Initialize classes for reading, exchanging and calculating
  // TestExchange is initialized with a small deadSize to send away
  distribute.setParameters(inFile, rL, dataType);
  test.setParameters(rL, deadSize / 100.0);
  exchange.setParameters(rL, deadSize);

  distribute.initialize();
  test.initialize();
  exchange.initialize();

  vector<POSVEL_T>* xx = new vector<POSVEL_T>;
  vector<POSVEL_T>* yy = new vector<POSVEL_T>;
  vector<POSVEL_T>* zz = new vector<POSVEL_T>;
  vector<POSVEL_T>* vx = new vector<POSVEL_T>;
  vector<POSVEL_T>* vy = new vector<POSVEL_T>;
  vector<POSVEL_T>* vz = new vector<POSVEL_T>;
  vector<ID_T>* tag = new vector<ID_T>;
  vector<STATUS_T>* status = new vector<STATUS_T>;

  distribute.setParticles(xx,yy,zz,vx,vy,vz,tag,status);

  // Read alive particles only from files
  // In ONE_TO_ONE every processor reads its own processor in the topology
  // which has already been populated with the correct alive particles
  distribute.readParticlesOneToOne();

  // Create the mask and potential vectors which will be filled in elsewhere
  int numberOfParticles = (*xx).size();
  vector<POTENTIAL_T>* potential = new vector<POTENTIAL_T>(numberOfParticles);
  vector<MASK_T>* mask = new vector<MASK_T>(numberOfParticles);

  // Output files for testing
  ostringstream outFile1, outFile2, outFile3;
  outFile1 << "dist." << Partition::getMyProc();
  outFile2 << "test." << Partition::getMyProc();
  outFile3 << "init." << Partition::getMyProc();
  ofstream distStr(outFile1.str().c_str());
  ofstream testStr(outFile2.str().c_str());
  ofstream initStr(outFile3.str().c_str());

  // Write particles after the ParticleDistribute
  for (int p = 0; p < xx->size(); p++) {
    distStr << (*tag)[p] << " " 
            << (*xx)[p] << " " << (*yy)[p] << " " << (*zz)[p] << " "
            << (*vx)[p] << " " << (*vy)[p] << " " << (*vz)[p] << endl;
  }

  // TestExchange will take some alive particles on a processor and
  // send them to the neighbor as ALIVE, deleting them from this processor
  test.setParticles(xx, yy, zz, vx, vy, vz, tag, status);
  test.exchangeParticles();

  // Write into arrays because that is what the initializer will do
  int particleCount = xx->size();
  POSVEL_T* xxInit = new POSVEL_T[particleCount];
  POSVEL_T* yyInit = new POSVEL_T[particleCount];
  POSVEL_T* zzInit = new POSVEL_T[particleCount];
  POSVEL_T* vxInit = new POSVEL_T[particleCount];
  POSVEL_T* vyInit = new POSVEL_T[particleCount];
  POSVEL_T* vzInit = new POSVEL_T[particleCount];
  ID_T* tagInit = new ID_T[particleCount];
  POTENTIAL_T* potInit = new POTENTIAL_T[particleCount];
  MASK_T* maskInit = new MASK_T[particleCount];
  
  for (int p = 0; p < particleCount; p++) {
    xxInit[p] = (*xx)[p];
    yyInit[p] = (*yy)[p];
    zzInit[p] = (*zz)[p];
    vxInit[p] = (*vx)[p];
    vyInit[p] = (*vy)[p];
    vzInit[p] = (*vz)[p];
    tagInit[p] = (*tag)[p];
    potInit[p] = p * 1.0;
    maskInit[p] = 011;
  }

  // Write particles after the TestExchange
  for (int p = 0; p < particleCount; p++) {
    testStr << tagInit[p] << " " 
            << xxInit[p] << " " << yyInit[p] << " " << zzInit[p] << " "
            << vxInit[p] << " " << vyInit[p] << " " << vzInit[p] << endl;
  }

  // Clear the original vectors for the final answer
  xx->clear();
  yy->clear();
  zz->clear();
  vx->clear();
  vy->clear();
  vz->clear();
  tag->clear();
  status->clear();
  potential->clear();
  mask->clear();

  // InitialExchange is called and should send every particle back
  exchange.setParticleArrays(particleCount,
                             xxInit, yyInit, zzInit, 
                             vxInit, vyInit, vzInit, potInit,
                             tagInit, maskInit);
  exchange.setParticleVectors(xx, yy, zz, vx, vy, vz, potential, 
                              tag, mask, status);
  exchange.exchangeParticles();

  // Write particles after the InitialExchange
  for (int p = 0; p < xx->size(); p++) {
    initStr << (*tag)[p] << " " 
            << (*xx)[p] << " " << (*yy)[p] << " " << (*zz)[p] << " "
            << (*vx)[p] << " " << (*vy)[p] << " " << (*vz)[p] << endl;
  }

  // Shut down MPI
  Partition::finalize();
  MPI_Finalize();

  delete [] xxInit;
  delete [] yyInit;
  delete [] zzInit;
  delete [] vxInit;
  delete [] vyInit;
  delete [] vzInit;
  delete [] tagInit;
  delete [] potInit;
  delete [] maskInit;

  return 0;
}
