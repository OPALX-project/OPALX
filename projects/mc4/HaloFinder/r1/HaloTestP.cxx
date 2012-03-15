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
#include "HaloProperties.h"

#include <mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc != 11) {
    cout << "Usage: mpirun -np # HaloTestP in out boxsz Hubble d np bb pmin ";
    cout << "[RECORD|BLOCK] [ROUND_ROBIN|ONE_TO_ONE]" << endl;
  }

  // Base file name (Actual file names are basename.proc#
  int i = 1;
  string inFile = argv[i++];
  string outFile = argv[i++];

  // Physical coordinate box size
  POSVEL_T boxSize = (POSVEL_T) (atof(argv[i++]));
  POSVEL_T hubbleConstant = (POSVEL_T) (atof(argv[i++]));
  POSVEL_T rL = boxSize / hubbleConstant;

  // Physical coordinate dead zone area
  POSVEL_T deadSize = (POSVEL_T) (atof(argv[i++]));

  // Superimposed grid on physical box used to determine wraparound
  int np = atoi(argv[i++]);

  // BB parameter for distance between particles comprising a halo
  POSVEL_T bb = (POSVEL_T) (atof(argv[i++]));

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
  //Partition::initialize(argc, argv);
  MPI_Init(&argc, &argv);
  Partition::initialize();

  // Construct the particle distributor, exchanger and halo finder
  ParticleDistribute distribute;
  ParticleExchange exchange;
  CosmoHaloFinderP haloFinder;

  // Initialize classes for reading, exchanging and calculating
  distribute.setParameters(inFile, rL, dataType);
  exchange.setParameters(rL, deadSize);
  haloFinder.setParameters(outFile, rL, deadSize, np, pmin, bb);

  distribute.initialize();
  exchange.initialize();

  // Read alive particles only from files
  // In ROUND_ROBIN all files are read and particles are passed round robin
  // to every other processor so that every processor chooses its own
  // In ONE_TO_ONE every processor reads its own processor in the topology
  // which has already been populated with the correct alive particles

/*
  double time1, time2;
  MPI_Barrier(MPI_COMM_WORLD);
  if (Partition::getMyProc() == MASTER)
    time1 = MPI_Wtime();
*/

  vector<POSVEL_T>* xx = new vector<POSVEL_T>;
  vector<POSVEL_T>* yy = new vector<POSVEL_T>;
  vector<POSVEL_T>* zz = new vector<POSVEL_T>;
  vector<POSVEL_T>* vx = new vector<POSVEL_T>;
  vector<POSVEL_T>* vy = new vector<POSVEL_T>;
  vector<POSVEL_T>* vz = new vector<POSVEL_T>;
  vector<ID_T>* tag = new vector<ID_T>;
  vector<STATUS_T>* status = new vector<STATUS_T>;

  distribute.setParticles(xx,yy,zz,vx,vy,vz,tag);

  if (distributeType == "ROUND_ROBIN")
    distribute.readParticlesRoundRobin();
  else if (distributeType == "ONE_TO_ONE")
    distribute.readParticlesOneToOne();

/*
  MPI_Barrier(MPI_COMM_WORLD);
  if (Partition::getMyProc() == MASTER) {
    time2 = MPI_Wtime();
    cout << "Time to distribute alive particles: " << (time2 - time1) << endl;
  }
*/

  // Create the mask and potential vectors which will be filled in elsewhere
  int numberOfParticles = (*xx).size();
  vector<POTENTIAL_T>* potential = new vector<POTENTIAL_T>(numberOfParticles);
  vector<MASK_T>* mask = new vector<MASK_T>(numberOfParticles);

  // Exchange particles adds dead particles to all the vectors
/*
  MPI_Barrier(MPI_COMM_WORLD);
  if (Partition::getMyProc() == MASTER)
    time1 = MPI_Wtime();
*/

  exchange.setParticles(xx, yy, zz, vx, vy, vz, potential, tag, mask, status);
  exchange.exchangeParticles();

/*
  MPI_Barrier(MPI_COMM_WORLD);
  if (Partition::getMyProc() == MASTER) {
    time2 = MPI_Wtime();
    cout << "Time to exchange dead particles: " << (time2 - time1) << endl;
  }
*/

  // Set alive and dead particle information per processor
  haloFinder.setParticles(xx, yy, zz, vx, vy, vz, potential, tag, mask, status);

  // Run serial halo finder on alive and dead particles for this processor
/*
  MPI_Barrier(MPI_COMM_WORLD);
  if (Partition::getMyProc() == MASTER)
    time1 = MPI_Wtime();
*/

  haloFinder.executeHaloFinder();

  // Collect the serial halo finder results
  haloFinder.collectHalos();

  // Merge the halos so that only one copy of each is written
  // Parallel halo finder must consult with each of the 26 possible neighbor
  // halo finders to see who will report a particular halo
  haloFinder.mergeHalos();

/*
  MPI_Barrier(MPI_COMM_WORLD);
  if (Partition::getMyProc() == MASTER) {
    time2 = MPI_Wtime();
    cout << "Time to find halos: " << (time2 - time1) << endl;
  }
*/

  // Collect information from the halo finder needed for halo properties
  // Vector halos is the index of the first particle for halo in the haloList
  // Following the chain of indices in the haloList retrieves all particles
  int numberOfHalos = haloFinder.getNumberOfHalos();
  int* halos = haloFinder.getHalos();
  int* haloCount = haloFinder.getHaloCount();
  int* haloList = haloFinder.getHaloList();

  // Construct the halo center finder giving it references to halo data
  POSVEL_T chainSize = 2.0;
  HaloProperties property;
  property.setHalos(numberOfHalos, halos, haloCount, haloList);
  property.setParameters(outFile, rL, deadSize);
  property.setParticles(xx, yy, zz, vx, vy, vz, potential, tag, mask, status);
  property.createChainingMesh(rL, deadSize, chainSize);

  // Following method is just to show how to iterate over the chaining mesh
  property.chainingMeshCentroids();

  // Not needed but here to regroup before doing the SOD calculation
  //MPI_Barrier(MPI_COMM_WORLD);

  // Find the index of the particle at the center of each halo
  vector<int>* haloCenter = new vector<int>;
  property.findHaloCenterN2_2(haloCenter);

  // Find the mass of each halo
  POSVEL_T particleMass = 1.0;
  vector<POSVEL_T>* haloMass = new vector<POSVEL_T>;
  property.findHaloMass(particleMass, haloMass);

  // Find the average position of each halo
  vector<POSVEL_T>* xAvgPos = new vector<POSVEL_T>;
  vector<POSVEL_T>* yAvgPos = new vector<POSVEL_T>;
  vector<POSVEL_T>* zAvgPos = new vector<POSVEL_T>;
  property.findAveragePosition(xAvgPos, yAvgPos, zAvgPos);

  // Find the average velocity of each halo
  vector<POSVEL_T>* xAvgVel = new vector<POSVEL_T>;
  vector<POSVEL_T>* yAvgVel = new vector<POSVEL_T>;
  vector<POSVEL_T>* zAvgVel = new vector<POSVEL_T>;
  property.findAverageVelocity(xAvgVel, yAvgVel, zAvgVel);

  //MPI_Barrier(MPI_COMM_WORLD);

  // Spherical over-density calculation
  int minHaloSize = 1000;
  property.calculateSOD(minHaloSize, particleMass, haloCenter);

  // Write halo results on each processor
  haloFinder.writeHalos();

  // Shut down MPI
  Partition::finalize();
  MPI_Finalize();

  return 0;
}
