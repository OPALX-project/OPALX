/*=========================================================================
                                                                                
Copyright (c) 2007, Los Alamos National Security, LLC

All rights reserved.

Copyright 2007. Los Alamos National Security, LLC. 
This software was produced under U.S. Government contract DE-AC52-06NA25396 
for Los Alamos National Laboratory (LANL), which is operated by 
Los Alamos National Security, LLC for the U.S. Department of Energy. 
The U.S. Government has rights to use, reproduce, and distribute this software. 
NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,
EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
If software is modified to produce derivative works, such modified software 
should be clearly marked, so as not to confuse it with the version available 
from LANL.
 
Additionally, redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following conditions 
are met:
-   Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer. 
-   Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution. 
-   Neither the name of Los Alamos National Security, LLC, Los Alamos National
    Laboratory, LANL, the U.S. Government, nor the names of its contributors
    may be used to endorse or promote products derived from this software 
    without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
                                                                                
=========================================================================*/

// .NAME HaloProperties - calculate properties of all halos
//
// .SECTION Description
// HaloProperties takes data from CosmoHaloFinderP about individual halos
// and data from all particles and calculates properties.
//

#ifndef HaloProperties_h
#define HaloProperties_h

#include "Definition.h"
#include <string>
#include <vector>

using namespace std;

class HaloProperties {
public:
  HaloProperties();
  ~HaloProperties();

  // Set parameters for sizes of the dead/alive space
  void setParameters(
	const string& outName,	// Base name of output halo files
	POSVEL_T rL,		// Box size of the physical problem
	POSVEL_T deadSize);	// Dead size used to normalize for non periodic

  // Set alive particle vectors which were created elsewhere
#ifdef MC4HALOFINDER
    void setParticles(
                        POSVEL_T* xLoc,
			POSVEL_T* yLoc,
                        POSVEL_T* zLoc,
                        POSVEL_T* xVel,
                        POSVEL_T* yVel,
                        POSVEL_T* zVel,
                        POTENTIAL_T* potential,
                        ID_T* id,
                        MASK_T* maskData,
                        STATUS_T* state,
			int np);
#else
  void setParticles(
        vector<POSVEL_T>* xLoc,
        vector<POSVEL_T>* yLoc,
        vector<POSVEL_T>* zLod,
        vector<POSVEL_T>* xVel,
        vector<POSVEL_T>* yVel,
        vector<POSVEL_T>* zVel,
        vector<POTENTIAL_T>* potential,
        vector<ID_T>* id,
        vector<MASK_T>* mask,
        vector<STATUS_T>* state);
#endif
  // Set the halo information from the FOF halo finder
  void setHalos(
	int  numberOfHalos,	// Number of halos found
	int* halos,		// Index into haloList of first particle
	int* haloCount,		// Number of particles in the matching halo
	int* haloList);		// Chain of indices of all particles in halo

  /////////////////////////////////////////////////////////////////////
  //
  // FOF (Friend of friends) halo analysis
  //
  /////////////////////////////////////////////////////////////////////

  // Find the halo centers (minimum potential) finding minimum of array
  void findHaloCenter(vector<int>* haloCenter);

  // Find the halo centers (minimum potential) using N^2 algorithm
  void findHaloCenterN2_2(vector<int>* haloCenter);
  int  minimumPotentialN2_2(int h, POTENTIAL_T* minPotential);

  // Find the mass of each halo
  void findHaloMass(
	POSVEL_T particleMass,
	vector<POSVEL_T>* haloMass);

  // Find the average position of a halo's particles
  void findAveragePosition(
	vector<POSVEL_T>* xPos,
	vector<POSVEL_T>* yPos,
	vector<POSVEL_T>* zPos);

  // Find the average velocity of a halo's particles
  void findAverageVelocity(
	vector<POSVEL_T>* xVel,
	vector<POSVEL_T>* yVel,
	vector<POSVEL_T>* zVel);

  // Kahan summation of floating point numbers to reduce roundoff error
  double KahanSummation(int halo, POSVEL_T* data);

  // Incremental mean, possibly needed for very large halos
  POSVEL_T incrementalMean(int halo, POSVEL_T* data);

  /////////////////////////////////////////////////////////////////////
  //
  // SOD (Spherical over density) halo analysis
  //
  /////////////////////////////////////////////////////////////////////

  // Create the chaining mesh to organize particles by location
  void createChainingMesh(
	POSVEL_T rL,
	POSVEL_T deadSize,
	POSVEL_T chainSize);

  // Demonstration method to show how to iterate over the chaining mesh
  void chainingMeshCentroids();

  // Spherical over-density (SOD) mass profile, velocity dispersion
  void calculateSOD(
	int minHaloSize,
	double particleMass,
	vector<int>* haloCenter);

  void collectSODParticles(
	int particleCenter,
	double radius);

private:
  int    myProc;		// My processor number
  int    numProc;		// Total number of processors

  int    layoutSize[DIMENSION]; // Decomposition of processors
  int    layoutPos[DIMENSION];  // Position of this processor in decomposition

  string outFile;		// File of particles written by this processor

  POSVEL_T boxSize;		// Physical box size of the data set
  POSVEL_T deadSize;		// Border size for dead particles

  long   particleCount;		// Total particles on this processor

  POSVEL_T* xx;			// X location for particles on this processor
  POSVEL_T* yy;			// Y location for particles on this processor
  POSVEL_T* zz;			// Z location for particles on this processor
  POSVEL_T* vx;			// X velocity for particles on this processor
  POSVEL_T* vy;			// Y velocity for particles on this processor
  POSVEL_T* vz;			// Z velocity for particles on this processor
  POTENTIAL_T* pot;		// Particle potential
  ID_T* tag;			// Id tag for particles on this processor
  MASK_T* mask;			// Particle information
  STATUS_T* status;		// Particle is ALIVE or labeled with neighbor
				// processor index where it is ALIVE

  // Information about halos from FOF halo finder
  int  numberOfHalos;		// Number of halos found
  int* halos;			// First particle index into haloList
  int* haloCount;		// Size of each halo 
  int* haloList;		// Indices of next particle in halo

  // Information about particle location in chaining mesh buckets
  POSVEL_T minMine[DIMENSION];	// Lowest value of particle location
  POSVEL_T maxMine[DIMENSION];	// Highest value of particle location

  POSVEL_T chainSize;		// Grid size in chaining mesh
  int  meshSize[DIMENSION];	// Chaining mesh grid dimension

  int*** buckets;		// First particle index into bucketList
  int*** bucketCount;		// Size of each bucket 
  int* bucketList;		// Indices of next particle in halo

  vector<int> sphereParticles;	// Particles in one SOD halo
};

#endif
