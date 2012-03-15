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

// .NAME ParticleDistribute - distribute particles to processors
//
// .SECTION Description
// ParticleDistribute takes a series of data files containing RECORD style
// .cosmo data or Gadget style BLOCK data
// along with parameters defining the box size for the data and for
// determining halos within the particle data.  It distributes the data
// across processors including a healthy dead zone of particles belonging
// to neighbor processors.  By definition all halos can be determined
// completely for any processor because of this dead zone.  The serial
// halo finder is called on each processor.
//

#ifndef ParticleDistribute_h
#define ParticleDistribute_h

// This is needed when mpi is used in a program
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END

#include "Definition.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <mpi.h>

using namespace std;

class ParticleDistribute {
public:
  ParticleDistribute();
  ~ParticleDistribute();

  // Set parameters particle distribution
  void setParameters(
	const string& inName,	// Base file name to read from
	float rL,		// Box size of the physical problem
	string dataType);	// BLOCK or RECORD structured input data

  // Set neighbor processor numbers and calculate dead regions
  void initialize();

  // Read particle files per processor and share round robin with others
  // extracting only the alive particles
  void readParticlesRoundRobin();
  void partitionInputFiles();

  // Read one particle file per processor with alive particles 
  // and correct topology
  void readParticlesOneToOne();

  // Get particle counts for allocating buffers
  void findFileParticleCount();

  // Round robin version must buffer for MPI sends to other processors
  void readFromRecordFile(
	ifstream* inStream,	// Stream to read from
	int firstParticle,	// First particle index to read in this chunk
	int numberOfParticles,	// Number of particles to read in this chunk
	float* fblock,          // Buffer for read in data
	int* iblock,            // Buffer for read in data
	float* buffer1);        // Reordered data

  void readFromBlockFile(
	ifstream* inStream,	// Stream to read from
	int firstParticle,	// First particle index to read in this chunk
	int numberOfParticles,	// Number of particles to read in this chunk
	int totParticles,	// Total particles (used to get offset)
	float* fblock,          // Buffer for read in data
	int* iblock,            // Buffer for read in data
	float* buffer1);        // Reordered data

  // One to one version of read is simpler with no MPI buffering
  void readFromRecordFile();
  void readFromBlockFile();

  // Collect local alive particles from the input buffers
  void distributeParticles(
	float* fBlock,		// Data read from file into this block
	int* iBlock,		// Data read from file into this block
	float* buffer1,		// Double buffering for reads
	float* buffer2, 	// Double buffering for reads
	int bufferSize);	// Size of double buffers
  void collectLocalParticles(float* buffer);

  // Return data needed by other software
  int     getParticleCount()	{ return this->particleCount; }

  vector<float>* getXLocation()	{ return this->xx; }
  vector<float>* getYLocation()	{ return this->yy; }
  vector<float>* getZLocation()	{ return this->zz; }
  vector<float>* getXVelocity()	{ return this->vx; }
  vector<float>* getYVelocity()	{ return this->vx; }
  vector<float>* getZVelocity()	{ return this->vx; }
  vector<int>* getTag()		{ return this->tag; }
  vector<int>* getStatus()	{ return this->status; }

private:
  int    myProc;		// My processor number
  int    numProc;		// Total number of processors

  string baseFile;		// Base name of input particle files
  int    inputType;		// BLOCK or RECORD structure
  int    maxFiles;		// Maximum number of files per processor
  vector<string> inFiles;	// Files read by this processor
  vector<int> fileParticles;	// Number of particles in files on processor

  struct CosmoHeader cosmoHeader;// Gadget file header

  int    maxParticles;		// Largest number of particles in any file
  int    maxRead;		// Largest number of particles read at one time
  int    maxReadsPerFile;	// Max number of reads per file

  int    totalParticles;	// Number of particles on all files
  int    headerSize;		// For BLOCK files

  int    nextProc;		// Where to send buffers to be shared
  int    prevProc;		// Where to receive buffers from be shared
  int    numberOfFiles;		// Number of input files total
  int    processorsPerFile;	// Multiple processors read same file
  int    numberOfFileSends;	// Number of round robin sends to share buffers

  int    layoutSize[DIMENSION];	// Decomposition of processors
  int    layoutPos[DIMENSION];	// Position of this processor in decomposition

  int    np;			// Number of particles in the problem
  float  boxSize;		// Physical box size (rL)

  int    numberOfAliveParticles;

  int    particleCount;		// Running index used to store data
				// Ends up as the number of alive plus dead

  float  minAlive[DIMENSION];	// Minimum alive particle location on processor
  float  maxAlive[DIMENSION];	// Maximum alive particle location on processor

  int    neighbor[NUM_OF_NEIGHBORS];		// Neighbor processor ids

  vector<float>* xx;		// X location for particles on this processor
  vector<float>* yy;		// Y location for particles on this processor
  vector<float>* zz;		// Z location for particles on this processor
  vector<float>* vx;		// X velocity for particles on this processor
  vector<float>* vy;		// Y velocity for particles on this processor
  vector<float>* vz;		// Z velocity for particles on this processor
  vector<float>* ms;		// Mass for particles on this processor
  vector<int>* tag;		// Id tag for particles on this processor

  vector<int>* status;		// Particle is ALIVE or labeled with neighbor
				// processor index where it is ALIVE
};

#endif
