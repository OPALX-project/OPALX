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

// .NAME ParticleExchange - read or get pointer to alive particles on this
//                          process and exchange dead particles with neighbors
//
// .SECTION Description
// ParticleExchange reads one file per processor containing RECORD style
// .cosmo data or Gadget style BLOCK data (or receives a pointer to memory)
// along with parameters defining the box size for the data and for
// determining halos within the particle data.  It exchanges the particles
// between neighbors so that every processor has alive and dead particles.
// By definition all halos can be determined
// completely for any processor because of this dead zone.  The serial
// halo finder is called on each processor.
//

#ifndef ParticleExchange_h
#define ParticleExchange_h

// This is needed when mpi is used in a program
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END

#include "Definition.h"
#include <string>
#include <vector>
#include <mpi.h>

using namespace std;

class ParticleExchange {
public:
  ParticleExchange();
  ~ParticleExchange();

  // Set parameters particle distribution
  void setParameters(
	float rL,		// Box size of the physical problem
	float deadSize);	// Dead delta border for each processor

  // Set neighbor processor numbers and calculate dead regions
  void initialize();

  // Calculate physical range of alive particles which must be shared
  void calculateExchangeRegions();

  // Set alive particle vectors which were created elsewhere
  void setParticles(
	vector<float>* xx,
	vector<float>* yy,
	vector<float>* zz,
	vector<float>* vx,
	vector<float>* vy,
	vector<float>* vz,
	vector<int>* tag,
	vector<int>* status);

  // Identify and exchange alive particles which must be shared with neighbors
  void exchangeParticles();
  void identifyExchangeParticles();
  void exchangeNeighborParticles();
  void exchange(
	int sendTo,		// Neighbor to send particles to
	int recvFrom,		// Neighbor to receive particles from
	float* sendBuffer,
	float* recvBuffer,
	int bufferSize);

  // Return data needed by other software
  int getParticleCount()		{ return this->particleCount; }

private:
  int    myProc;		// My processor number
  int    numProc;		// Total number of processors

  int    totalParticles;	// Number of particles on all files
  int    headerSize;		// For BLOCK files

  int    layoutSize[DIMENSION];	// Decomposition of processors
  int    layoutPos[DIMENSION];	// Position of this processor in decomposition

  int    np;			// Number of particles in the problem
  float  boxSize;		// Physical box size (rL)
  float  deadSize;		// Border size for dead particles

  int    numberOfAliveParticles;
  int    numberOfDeadParticles;

  int    particleCount;		// Running index used to store data
				// Ends up as the number of alive plus dead

  float  minMine[DIMENSION];	// Minimum alive particle not exchanged
  float  maxMine[DIMENSION];	// Maximum alive particle not exchanged
  float  minShare[DIMENSION];	// Minimum alive particle shared
  float  maxShare[DIMENSION];	// Maximum alive particle shared

  int    neighbor[NUM_OF_NEIGHBORS];		// Neighbor processor ids
  float  minRange[NUM_OF_NEIGHBORS][DIMENSION];	// Range of dead particles
  float  maxRange[NUM_OF_NEIGHBORS][DIMENSION];	// Range of dead particles

  vector<int> neighborParticles[NUM_OF_NEIGHBORS];
				// Particle ids sent to each neighbor as DEAD

  vector<float>* xx;		// X location for particles on this processor
  vector<float>* yy;		// Y location for particles on this processor
  vector<float>* zz;		// Z location for particles on this processor
  vector<float>* vx;		// X velocity for particles on this processor
  vector<float>* vy;		// Y velocity for particles on this processor
  vector<float>* vz;		// Z velocity for particles on this processor
  vector<int>* tag;		// Id tag for particles on this processor

  vector<int>* status;		// Particle is ALIVE or labeled with neighbor
				// processor index where it is ALIVE
};

#endif
