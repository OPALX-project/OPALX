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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <sys/types.h>
#include <dirent.h>

#include "Partition.h"
#include "ParticleExchange.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////
//
// Particle data space is partitioned for the number of processors
// which currently is a factor of two but is easily extended.  Particles
// are read in from files where each processor reads one file into a buffer,
// extracts the particles which really belong on the processor (ALIVE) and
// those in a buffer region around the edge (DEAD).  The buffer is then
// passed round robin to every other processor so that all particles are
// examined by all processors.  All dead particles are tagged with the
// neighbor zone (26 neighbors in 3D) so that later halos can be associated
// with zones.
//
/////////////////////////////////////////////////////////////////////////

ParticleExchange::ParticleExchange()
{
  // Get the number of processors running this problem and rank
  this->numProc = Partition::getNumProc();
  this->myProc = Partition::getMyProc();

  // Get the number of processors in each dimension
  Partition::getDecompSize(this->layoutSize);

  // Get my position within the Cartesian topology
  Partition::getMyPosition(this->layoutPos);

  // Get neighbors of this processor including the wraparound
  Partition::getNeighbors(this->neighbor);

  this->numberOfAliveParticles = 0;
  this->numberOfDeadParticles = 0;
}

ParticleExchange::~ParticleExchange()
{
}

/////////////////////////////////////////////////////////////////////////
//
// Set parameters for particle distribution
//
/////////////////////////////////////////////////////////////////////////

void ParticleExchange::setParameters(float rL, float deadSz)
{
  // Physical total space and amount of physical space to use for dead particles
  this->boxSize = rL;
  this->deadSize = deadSz;

  if (this->myProc == MASTER) {
    cout << endl << "------------------------------------" << endl;
    cout << "boxSize:  " << this->boxSize << endl;
    cout << "deltaBox: " << this->deadSize << endl;
  }
}

/////////////////////////////////////////////////////////////////////////
//
// All particles on this processor initially are alive, but some of those
// alive must be exchanged with neighbors.  Determine the physical range
// on this processor where an ALIVE particle will never be exchanged and
// the ranges for each neighbor's future DEAD particles.  Then when
// reading each particle it can quickly be assigned.
//
/////////////////////////////////////////////////////////////////////////

void ParticleExchange::initialize()
{
#ifdef DEBUG
  if (this->myProc == MASTER)
    cout << "Decomposition: [" << this->layoutSize[0] << ":"
         << this->layoutSize[1] << ":" << this->layoutSize[2] << "]" << endl;
#endif

  // Set subextents on particle locations for this processor
  float boxStep[DIMENSION];
  for (int dim = 0; dim < DIMENSION; dim++) {
    boxStep[dim] = this->boxSize / this->layoutSize[dim];

    // All particles are alive and available for sharing
    this->minShare[dim] = this->layoutPos[dim] * boxStep[dim];
    this->maxShare[dim] = this->minShare[dim] + boxStep[dim];
    if (this->maxShare[dim] > this->boxSize)
      this->maxShare[dim] = this->boxSize;

    // Particles in the middle of the shared region will not be shared
    // and will only belong on this processor
    // If the layout is of size 1 in any dimension no need to exchange
    if (this->layoutSize[dim] > 1) {
      this->minMine[dim] = this->minShare[dim] + this->deadSize;
      this->maxMine[dim] = this->maxShare[dim] - this->deadSize;
    } else {
      this->minMine[dim] = this->minShare[dim];
      this->maxMine[dim] = this->maxShare[dim];
    }
  }

  // Set the ranges on the dead particles for each neighbor direction
  calculateExchangeRegions();
}

/////////////////////////////////////////////////////////////////////////
//
// Each of the 26 neighbors will be sent a rectangular region of my particles
// Calculate the range in each dimension of the ghost area
//
/////////////////////////////////////////////////////////////////////////

void ParticleExchange::calculateExchangeRegions()
{
  // Initialize all neighbors to the entire available exchange range
  for (int i = 0; i < NUM_OF_NEIGHBORS; i++) {
    for (int dim = 0; dim < DIMENSION; dim++) {
      this->minRange[i][dim] = this->minShare[dim];
      this->maxRange[i][dim] = this->maxShare[dim];
    }
  }

  // Left face
  this->minRange[X0][0] = this->minShare[0];
  this->maxRange[X0][0] = this->minMine[0];

  // Right face
  this->minRange[X1][0] = this->maxMine[0];
  this->maxRange[X1][0] = this->maxShare[0];

  // Bottom face
  this->minRange[Y0][1] = this->minShare[1];
  this->maxRange[Y0][1] = this->minMine[1];

  // Top face
  this->minRange[Y1][1] = this->maxMine[1];
  this->maxRange[Y1][1] = this->maxShare[1];

  // Front face
  this->minRange[Z0][2] = this->minShare[2];
  this->maxRange[Z0][2] = this->minMine[2];

  // Back face
  this->minRange[Z1][2] = this->maxMine[2];
  this->maxRange[Z1][2] = this->maxShare[2];

  // Left bottom and top bars
  this->minRange[X0_Y0][0] = this->minShare[0];
  this->maxRange[X0_Y0][0] = this->minMine[0];
  this->minRange[X0_Y0][1] = this->minShare[1];
  this->maxRange[X0_Y0][1] = this->minMine[1];

  this->minRange[X0_Y1][0] = this->minShare[0];
  this->maxRange[X0_Y1][0] = this->minMine[0];
  this->minRange[X0_Y1][1] = this->maxMine[1];
  this->maxRange[X0_Y1][1] = this->maxShare[1];

  // Right bottom and top bars
  this->minRange[X1_Y0][0] = this->maxMine[0];
  this->maxRange[X1_Y0][0] = this->maxShare[0];
  this->minRange[X1_Y0][1] = this->minShare[1];
  this->maxRange[X1_Y0][1] = this->minMine[1];

  this->minRange[X1_Y1][0] = this->maxMine[0];
  this->maxRange[X1_Y1][0] = this->maxShare[0];
  this->minRange[X1_Y1][1] = this->maxMine[1];
  this->maxRange[X1_Y1][1] = this->maxShare[1];

  // Bottom front and back bars
  this->minRange[Y0_Z0][1] = this->minShare[1];
  this->maxRange[Y0_Z0][1] = this->minMine[1];
  this->minRange[Y0_Z0][2] = this->minShare[2];
  this->maxRange[Y0_Z0][2] = this->minMine[2];

  this->minRange[Y0_Z1][1] = this->minShare[1];
  this->maxRange[Y0_Z1][1] = this->minMine[1];
  this->minRange[Y0_Z1][2] = this->maxMine[2];
  this->maxRange[Y0_Z1][2] = this->maxShare[2];

  // Top front and back bars 
  this->minRange[Y1_Z0][1] = this->maxMine[1];
  this->maxRange[Y1_Z0][1] = this->maxShare[1];
  this->minRange[Y1_Z0][2] = this->minShare[2];
  this->maxRange[Y1_Z0][2] = this->minMine[2];

  this->minRange[Y1_Z1][1] = this->maxMine[1];
  this->maxRange[Y1_Z1][1] = this->maxShare[1];
  this->minRange[Y1_Z1][2] = this->maxMine[2];
  this->maxRange[Y1_Z1][2] = this->maxShare[2];

  // Left front and back bars (vertical)
  this->minRange[Z0_X0][0] = this->minShare[0];
  this->maxRange[Z0_X0][0] = this->minMine[0];
  this->minRange[Z0_X0][2] = this->minShare[2];
  this->maxRange[Z0_X0][2] = this->minMine[2];

  this->minRange[Z1_X0][0] = this->minShare[0];
  this->maxRange[Z1_X0][0] = this->minMine[0];
  this->minRange[Z1_X0][2] = this->maxMine[2];
  this->maxRange[Z1_X0][2] = this->maxShare[2];

  // Right front and back bars (vertical)
  this->minRange[Z0_X1][0] = this->maxMine[0];
  this->maxRange[Z0_X1][0] = this->maxShare[0];
  this->minRange[Z0_X1][2] = this->minShare[2];
  this->maxRange[Z0_X1][2] = this->minMine[2];

  this->minRange[Z1_X1][0] = this->maxMine[0];
  this->maxRange[Z1_X1][0] = this->maxShare[0];
  this->minRange[Z1_X1][2] = this->maxMine[2];
  this->maxRange[Z1_X1][2] = this->maxShare[2];

  // Left bottom front corner
  this->minRange[X0_Y0_Z0][0] = this->minShare[0];
  this->maxRange[X0_Y0_Z0][0] = this->minMine[0];
  this->minRange[X0_Y0_Z0][1] = this->minShare[1];
  this->maxRange[X0_Y0_Z0][1] = this->minMine[1];
  this->minRange[X0_Y0_Z0][2] = this->minShare[2];
  this->maxRange[X0_Y0_Z0][2] = this->minMine[2];

  // Left bottom back corner
  this->minRange[X0_Y0_Z1][0] = this->minShare[0];
  this->maxRange[X0_Y0_Z1][0] = this->minMine[0];
  this->minRange[X0_Y0_Z1][1] = this->minShare[1];
  this->maxRange[X0_Y0_Z1][1] = this->minMine[1];
  this->minRange[X0_Y0_Z1][2] = this->maxMine[2];
  this->maxRange[X0_Y0_Z1][2] = this->maxShare[2];

  // Left top front corner
  this->minRange[X0_Y1_Z0][0] = this->minShare[0];
  this->maxRange[X0_Y1_Z0][0] = this->minMine[0];
  this->minRange[X0_Y1_Z0][1] = this->maxMine[1];
  this->maxRange[X0_Y1_Z0][1] = this->maxShare[1];
  this->minRange[X0_Y1_Z0][2] = this->minShare[2];
  this->maxRange[X0_Y1_Z0][2] = this->minMine[2];

  // Left top back corner
  this->minRange[X0_Y1_Z1][0] = this->minShare[0];
  this->maxRange[X0_Y1_Z1][0] = this->minMine[0];
  this->minRange[X0_Y1_Z1][1] = this->maxMine[1];
  this->maxRange[X0_Y1_Z1][1] = this->maxShare[1];
  this->minRange[X0_Y1_Z1][2] = this->maxMine[2];
  this->maxRange[X0_Y1_Z1][2] = this->maxShare[2];

  // Right bottom front corner
  this->minRange[X1_Y0_Z0][0] = this->maxMine[0];
  this->maxRange[X1_Y0_Z0][0] = this->maxShare[0];
  this->minRange[X1_Y0_Z0][1] = this->minShare[1];
  this->maxRange[X1_Y0_Z0][1] = this->minMine[1];
  this->minRange[X1_Y0_Z0][2] = this->minShare[2];
  this->maxRange[X1_Y0_Z0][2] = this->minMine[2];

  // Right bottom back corner
  this->minRange[X1_Y0_Z1][0] = this->maxMine[0];
  this->maxRange[X1_Y0_Z1][0] = this->maxShare[0];
  this->minRange[X1_Y0_Z1][1] = this->minShare[1];
  this->maxRange[X1_Y0_Z1][1] = this->minMine[1];
  this->minRange[X1_Y0_Z1][2] = this->maxMine[2];
  this->maxRange[X1_Y0_Z1][2] = this->maxShare[2];

  // Right top front corner
  this->minRange[X1_Y1_Z0][0] = this->maxMine[0];
  this->maxRange[X1_Y1_Z0][0] = this->maxShare[0];
  this->minRange[X1_Y1_Z0][1] = this->maxMine[1];
  this->maxRange[X1_Y1_Z0][1] = this->maxShare[1];
  this->minRange[X1_Y1_Z0][2] = this->minShare[2];
  this->maxRange[X1_Y1_Z0][2] = this->minMine[2];

  // Right top back corner
  this->minRange[X1_Y1_Z1][0] = this->maxMine[0];
  this->maxRange[X1_Y1_Z1][0] = this->maxShare[0];
  this->minRange[X1_Y1_Z1][1] = this->maxMine[1];
  this->maxRange[X1_Y1_Z1][1] = this->maxShare[1];
  this->minRange[X1_Y1_Z1][2] = this->maxMine[2];
  this->maxRange[X1_Y1_Z1][2] = this->maxShare[2];
}

/////////////////////////////////////////////////////////////////////////
//
// Set the particle vectors that have already been read and which
// contain only the alive particles for this processor
//
/////////////////////////////////////////////////////////////////////////

void ParticleExchange::setParticles(
			vector<float>* xLoc,
			vector<float>* yLoc,
			vector<float>* zLoc,
			vector<float>* xVel,
			vector<float>* yVel,
			vector<float>* zVel,
			vector<int>* id,
			vector<int>* type)
{
  this->particleCount = xLoc->size();
  this->numberOfAliveParticles = this->particleCount;
  this->xx = xLoc;
  this->yy = yLoc;
  this->zz = zLoc;
  this->vx = xVel;
  this->vy = yVel;
  this->vz = zVel;
  this->tag = id;
  this->status = type;
}
	
/////////////////////////////////////////////////////////////////////////////
//
// Alive particles are contained on each processor.  Identify the border
// particles which will be dead on other processors and exchange them
//
/////////////////////////////////////////////////////////////////////////////

void ParticleExchange::exchangeParticles()
{
  // Identify alive particles on this processor which must be shared
  // because they are dead particles on neighbor processors
  identifyExchangeParticles();

  // Exchange those particles with appropriate neighbors
  exchangeNeighborParticles();

  // Count the particles across processors
  int totalAliveParticles = 0;
  int totalDeadParticles = 0;
  MPI_Allreduce((void*) &this->numberOfAliveParticles, 
                (void*) &totalAliveParticles, 
                1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce((void*) &this->numberOfDeadParticles,
                (void*) &totalDeadParticles, 
                1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#ifdef DEBUG
  cout << "Exchange Particles Rank " << setw(3) << this->myProc 
       << " #alive = " << this->numberOfAliveParticles
       << " #dead = " << this->numberOfDeadParticles << endl;
#endif
 
  if (this->myProc == MASTER) {
    cout << "TotalAliveParticles " << totalAliveParticles << endl;
    cout << "TotalDeadParticles  " << totalDeadParticles << endl << endl;
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Iterate over all the alive particles on this processor and determine
// which must be shared and add them to the vector for that neighbor
//
/////////////////////////////////////////////////////////////////////////////

void ParticleExchange::identifyExchangeParticles()
{
  int notSharedCount = 0;
  int sharedCount = 0;

  for (int i = 0; i < this->particleCount; i++) {
    if (((*this->xx)[i] > this->minMine[0] && 
         (*this->xx)[i] < this->maxMine[0]) &&
        ((*this->yy)[i] > this->minMine[1] && 
         (*this->yy)[i] < this->maxMine[1]) &&
        ((*this->zz)[i] > this->minMine[2] && 
         (*this->zz)[i] < this->maxMine[2])) {
          notSharedCount++;
    } else {
      // Particle is alive here but which processors need it as dead
      for (int n = 0; n < NUM_OF_NEIGHBORS; n++) {
        if ((*this->xx)[i] >= minRange[n][0] && 
            (*this->xx)[i] <= maxRange[n][0] &&
            (*this->yy)[i] >= minRange[n][1] && 
            (*this->yy)[i] <= maxRange[n][1] &&
            (*this->zz)[i] >= minRange[n][2] && 
            (*this->zz)[i] <= maxRange[n][2]) {
                this->neighborParticles[n].push_back(i);
                sharedCount++;
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Exchange the appropriate particles with neighbors
// Only the index of the particle to be exchanged is stored so fill out
// the message with location, velocity, tag.  Status information doesn't
// have to be sent because when the message is received, the neighbor
// containing the new dead particle will be known
//
// Use the Cartesian communicator for neighbor exchange
//
/////////////////////////////////////////////////////////////////////////////

void ParticleExchange::exchangeNeighborParticles()
{
  // Calculate the maximum number of particles to share for calculating buffer
  int myShareSize = 0;
  for (int n = 0; n < NUM_OF_NEIGHBORS; n++)
    if (myShareSize < this->neighborParticles[n].size())
      myShareSize = this->neighborParticles[n].size();


  int maxShareSize;
  MPI_Allreduce((void*) &myShareSize,
                (void*) &maxShareSize,
                1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // Allocate send and receive MPI buffers
  // Use the last float which was mass to store the tag
  // First float holds the number of particles being sent

  int bufferSize = 1 + (maxShareSize * COSMO_FLOAT);
  float* sendBuffer = new float[bufferSize];
  float* recvBuffer = new float[bufferSize];

  // Exchange with each neighbor, with everyone sending in one direction and
  // receiving from the other.  Data corresponding to the particle index
  // must be packed in the buffer.  When the data is received it is unpacked
  // into the location, velocity and tag vectors and the status is set
  // to the neighbor who sent it

  for (int n = 0; n < NUM_OF_NEIGHBORS; n=n+2) {
    // Neighbor pairs in Definition.h must match so that every processor
    // sends and every processor receives on each exchange
    exchange(n, n+1, sendBuffer, recvBuffer, bufferSize);
    exchange(n+1, n, sendBuffer, recvBuffer, bufferSize);
  }

  delete [] sendBuffer;
  delete [] recvBuffer;
}

/////////////////////////////////////////////////////////////////////////////
//
// Pack particle data for the indicated neighbor into MPI message
// Send that message and receive from opposite neighbor
// Unpack the received particle data and add to particle buffers with
// an indication of dead and the neighbor on which particle is alive
//
/////////////////////////////////////////////////////////////////////////////

void ParticleExchange::exchange(
			int sendTo, 
			int recvFrom, 
			float* sendBuffer, 
			float* recvBuffer, 
			int bufferSize)
{
  // Pack the send buffer
  int sendParticleCount = this->neighborParticles[sendTo].size();
  sendBuffer[0] = *(reinterpret_cast<float*>(&sendParticleCount));

  int index = 1;
  for (int i = 0; i < sendParticleCount; i++) {
    int deadIndex = this->neighborParticles[sendTo][i];
    sendBuffer[index++] = (*this->xx)[deadIndex];
    sendBuffer[index++] = (*this->yy)[deadIndex];
    sendBuffer[index++] = (*this->zz)[deadIndex];
    sendBuffer[index++] = (*this->vx)[deadIndex];
    sendBuffer[index++] = (*this->vy)[deadIndex];
    sendBuffer[index++] = (*this->vz)[deadIndex];
    sendBuffer[index++] = *(reinterpret_cast<float*>(&(*this->tag)[deadIndex]));
  }

  // Send the buffer
  MPI_Request mpiRequest;
  MPI_Isend(sendBuffer, bufferSize, MPI_FLOAT, this->neighbor[sendTo],
            0, Partition::getComm(), &mpiRequest);

  // Receive the buffer from neighbor on other side
  MPI_Status mpiStatus;
  MPI_Recv(recvBuffer, bufferSize, MPI_FLOAT, this->neighbor[recvFrom],
           0, Partition::getComm(), &mpiStatus);

  MPI_Barrier(Partition::getComm());

  // Process the received buffer
  int recvParticleCount = *(reinterpret_cast<int*>(&recvBuffer[0]));

  index = 1;
  for (int i = 0; i < recvParticleCount; i++) {
    this->xx->push_back(recvBuffer[index++]);
    this->yy->push_back(recvBuffer[index++]);
    this->zz->push_back(recvBuffer[index++]);
    this->vx->push_back(recvBuffer[index++]);
    this->vy->push_back(recvBuffer[index++]);
    this->vz->push_back(recvBuffer[index++]);
    this->tag->push_back(*(reinterpret_cast<int*>(&recvBuffer[index++])));
    this->status->push_back(recvFrom);

    this->numberOfDeadParticles++;
    this->particleCount++;
  }
}
