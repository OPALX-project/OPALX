/*=========================================================================

  Program:   HaloFinder driver program
  Module:    $RCSfile: HaloTestP.cxx,v $

  Copyright (c) Chung-Hsing Hsu
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This example shows the driver for the HaloFinder parallel test. 
// It takes either a .cosmo file or Gadget file arranged for round robin
// passing of particles or dedicated single file per processor and outputs
// two files, one with all particles and the other with the halo catalog.
//

#include <iostream>
#include "Partition.h"
#include "ParticleDistribute.h"
#include "ParticleExchange.h"
#include "CosmoHaloFinderP.h"

using namespace std;

int main(int argc, char* argv[])
{
  if (argc != 10) {
    cout << "Usage: mpirun -np # HaloTestP in out rL d np bb pmin ";
    cout << "[RECORD|BLOCK] [ROUND_ROBIN|ONE_TO_ONE]" << endl;
  }

  // Base file name (Actual file names are basename.proc#
  int i = 1;
  string inFile = argv[i++];
  string outFile = argv[i++];

  // Physical coordinate box size
  float rL = atof(argv[i++]);

  // Physical coordinate dead zone area
  float deadSize = atof(argv[i++]);

  // Superimposed grid on physical box used to determine wraparound
  int np = atoi(argv[i++]);

  // BB parameter for distance between particles comprising a halo
  float bb = atof(argv[i++]);

  // Minimum number of particles to make a halo
  int pmin = atoi(argv[i++]);

  // Input is BLOCK or RECORD structured data
  string dataType = argv[i++];

  // Distribution mechanism is ROUND_ROBIN which selects alive and dead
  // as the data files are read and passed around
  // or EXCHANGE where the alive data is immediately read on a processor
  // and the dead particles must be bundled and shared
  string distributeType = argv[i++];

  // Initialize the partitioner which uses MPI Cartesian Topology
  Partition::initialize(argc, argv);

  // Construct the particle distributor, exchanger and halo finder
  ParticleDistribute distribute;
  ParticleExchange exchange;
  CosmoHaloFinderP haloFinder;

  // Initialize classes for reading, exchanging and calculating
  distribute.setParameters(inFile, rL, dataType);
  exchange.setParameters(rL, deadSize);
  haloFinder.setParameters(outFile, rL, np, pmin, bb);

  distribute.initialize();
  exchange.initialize();

  // Read alive particles only from files
  // In ROUND_ROBIN all files are read and particles are passed round robin
  // to every other processor so that every processor chooses its own
  // In ONE_TO_ONE every processor reads its own processor in the topology
  // which has already been populated with the correct alive particles

  if (distributeType == "ROUND_ROBIN")
    distribute.readParticlesRoundRobin();
  else if (distributeType == "ONE_TO_ONE")
    distribute.readParticlesOneToOne();

  // Collect alive particle information per processor
  vector<float>* xx = distribute.getXLocation();
  vector<float>* yy = distribute.getYLocation();
  vector<float>* zz = distribute.getZLocation();
  vector<float>* vx = distribute.getXVelocity();
  vector<float>* vy = distribute.getYVelocity();
  vector<float>* vz = distribute.getZVelocity();
  vector<int>* tag = distribute.getTag();
  vector<int>* status = distribute.getStatus();

  // Exchange particles adds dead particles to all the vectors
  exchange.setParticles(xx, yy, zz, vx, vy, vz, tag, status);
  exchange.exchangeParticles();

  // Set alive and dead particle information per processor
  haloFinder.setParticles(xx, yy, zz, vx, vy, vz, tag, status);

  // Run serial halo finder on alive and dead particles for this processor
  haloFinder.executeHaloFinder();

  // Collect the serial halo finder results
  haloFinder.collectHalos();

  // Merge the halos so that only one copy of each is written
  // Parallel halo finder must consult with each of the 26 possible neighbor
  // halo finders to see who will report a particular halo
  haloFinder.mergeHalos();

  // Write halo results on each processor
  haloFinder.writeHalos();

  // Shut down MPI
  Partition::finalize();

  return 0;
}
