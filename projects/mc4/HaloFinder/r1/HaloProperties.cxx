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
#include <set>
#include <math.h>
#ifdef MC4HALOFINDER

#else
#include "Partition.h"
#endif

#include "HaloProperties.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////
//
// HaloProperties uses the results of the CosmoHaloFinder to locate the
// particle within every halo in order to calculate properties on halos
//
/////////////////////////////////////////////////////////////////////////

HaloProperties::HaloProperties()
{
#ifdef MC4HALOFINDER

#else
  // Get the number of processors and rank of this processor
  this->numProc = Partition::getNumProc();
  this->myProc = Partition::getMyProc();

  // Get the number of processors in each dimension
  Partition::getDecompSize(this->layoutSize);

  // Get my position within the Cartesian topology
  Partition::getMyPosition(this->layoutPos);
#endif
}

HaloProperties::~HaloProperties()
{
}

/////////////////////////////////////////////////////////////////////////
//
// Set chaining structure which will locate all particles in a halo
//
/////////////////////////////////////////////////////////////////////////

void HaloProperties::setHalos(
			int numberHalos,
			int* haloStartIndex,
			int* haloParticleCount,
			int* nextParticleIndex)
{
  this->numberOfHalos = numberHalos;
  this->halos = haloStartIndex;
  this->haloCount = haloParticleCount;
  this->haloList = nextParticleIndex;
}

/////////////////////////////////////////////////////////////////////////
//
// Set parameters for the halo center finder
//
/////////////////////////////////////////////////////////////////////////

void HaloProperties::setParameters(
			const string& outName,
			POSVEL_T rL,
			POSVEL_T deadSz)
{
  // Particles for this processor output to file
  ostringstream oname, hname;
  if (this->numProc == 1) {
    oname << outName;
    hname << outName;
  } else {
    oname << outName << "." << myProc;
    hname << outName << ".halo." << myProc;
  }
  this->outFile = oname.str();

  // Halo finder parameters
  this->boxSize = rL;
  this->deadSize = deadSz;
}

/////////////////////////////////////////////////////////////////////////
//
// Set the particle vectors that have already been read and which
// contain only the alive particles for this processor
//
/////////////////////////////////////////////////////////////////////////
#ifdef MC4HALOFINDER
void HaloProperties::setParticles(
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
			int np)
{
    this->particleCount = np;
    this->xx = xLoc;
    this->yy = yLoc;
    this->zz = zLoc;
    this->vx = xVel;
    this->vy = yVel;
    this->vz = zVel;
    this->pot = potential;
    this->tag = id;
    this->mask = maskData;
    this->status = state;
}
#else
void HaloProperties::setParticles(
                        vector<POSVEL_T>* xLoc,
                        vector<POSVEL_T>* yLoc,
                        vector<POSVEL_T>* zLoc,
                        vector<POSVEL_T>* xVel,
                        vector<POSVEL_T>* yVel,
                        vector<POSVEL_T>* zVel,
                        vector<POTENTIAL_T>* potential,
                        vector<ID_T>* id,
                        vector<MASK_T>* maskData,
                        vector<STATUS_T>* state)
{
  this->particleCount = xLoc->size();

  // Extract the contiguous data block from a vector pointer
  this->xx = &(*xLoc)[0];
  this->yy = &(*yLoc)[0];
  this->zz = &(*zLoc)[0];
  this->vx = &(*xVel)[0];
  this->vy = &(*yVel)[0];
  this->vz = &(*zVel)[0];
  this->pot = &(*potential)[0];
  this->tag = &(*id)[0];
  this->mask = &(*maskData)[0];
  this->status = &(*state)[0];
}
#endif
/////////////////////////////////////////////////////////////////////////
//
// Find the index of the particle at the center of each halo and return
// in the given vector.  It is the minimum value in the potential array.
//
/////////////////////////////////////////////////////////////////////////

void HaloProperties::findHaloCenter(vector<int>* haloCenter)
{
  for (int halo = 0; halo < this->numberOfHalos; halo++) {
                        
    // First particle in halo
    int p = this->halos[halo];
    POTENTIAL_T minPotential = this->pot[p];
    int centerIndex = p;
  
    // Next particle
    p = this->haloList[p];
  
    // Search for minimum
    while (p != -1) {
      if (minPotential > this->pot[p]) {
        minPotential = this->pot[p];
        centerIndex = p;
      }
      p = this->haloList[p];
    }

    // Save the minimum potential index for this halo
    (*haloCenter).push_back(centerIndex);
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Find the index of the particle at the center of each halo and return
// in the given vector.  Use the minimumPotential() algorithm
//
/////////////////////////////////////////////////////////////////////////

void HaloProperties::findHaloCenterN2_2(vector<int>* haloCenter)
{
  int smallHalo = 0;
  int largeHalo = 0;
  POTENTIAL_T minPotential;

  for (int halo = 0; halo < this->numberOfHalos; halo++) {
    if (this->haloCount[halo] <= 100) {
      smallHalo++;
      int centerIndex = minimumPotentialN2_2(halo, &minPotential);
      (*haloCenter).push_back(centerIndex);
    } else {
      largeHalo++;
      (*haloCenter).push_back(-1);
    }
  }
  cout << "Rank " << this->myProc 
       << " #small = " << smallHalo
       << " #large = " << largeHalo << endl;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate the mass of each halo
//
/////////////////////////////////////////////////////////////////////////

void HaloProperties::findHaloMass(
			POSVEL_T particleMass,
			vector<POSVEL_T>* haloMass)
{
  for (int halo = 0; halo < this->numberOfHalos; halo++) {
    int p = this->halos[halo];
    while (p != -1) {
      POSVEL_T mass = particleMass * this->haloCount[halo];
      (*haloMass).push_back(mass);
      p = this->haloList[p];
    }
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate the average position of particles in every halo
//
/////////////////////////////////////////////////////////////////////////

void HaloProperties::findAveragePosition(
			vector<POSVEL_T>* xMeanPos,
			vector<POSVEL_T>* yMeanPos,
			vector<POSVEL_T>* zMeanPos)
{
  POSVEL_T xMean, yMean, zMean;
  double xKahan, yKahan, zKahan;

  for (int halo = 0; halo < this->numberOfHalos; halo++) {
    xKahan = KahanSummation(halo, this->xx);
    yKahan = KahanSummation(halo, this->yy);
    zKahan = KahanSummation(halo, this->zz);

    xMean = (POSVEL_T) (xKahan / this->haloCount[halo]);
    yMean = (POSVEL_T) (yKahan / this->haloCount[halo]);
    zMean = (POSVEL_T) (zKahan / this->haloCount[halo]);

    (*xMeanPos).push_back(xMean);
    (*yMeanPos).push_back(yMean);
    (*zMeanPos).push_back(zMean);
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate the average velocity of particles in every halo
//
/////////////////////////////////////////////////////////////////////////

void HaloProperties::findAverageVelocity(
                        vector<POSVEL_T>* xMeanVel,
                        vector<POSVEL_T>* yMeanVel,
                        vector<POSVEL_T>* zMeanVel)
{
  POSVEL_T xMean, yMean, zMean;
  double xKahan, yKahan, zKahan;

  for (int halo = 0; halo < this->numberOfHalos; halo++) {
    xKahan = KahanSummation(halo, this->vx);
    yKahan = KahanSummation(halo, this->vy);
    zKahan = KahanSummation(halo, this->vz);

    xMean = (POSVEL_T) (xKahan / this->haloCount[halo]);
    yMean = (POSVEL_T) (yKahan / this->haloCount[halo]);
    zMean = (POSVEL_T) (zKahan / this->haloCount[halo]);

    (*xMeanVel).push_back(xMean);
    (*yMeanVel).push_back(yMean);
    (*zMeanVel).push_back(zMean);
  }
}

/////////////////////////////////////////////////////////////////////////
//                      
// Calculate the Kahan summation
// Reduces roundoff error in floating point arithmetic
//                      
/////////////////////////////////////////////////////////////////////////
                        
double HaloProperties::KahanSummation(int halo, POSVEL_T* data)
{                       
  double dataSum, dataRem, value, v, w;
                        
  // First particle in halo and first step in Kahan summation
  int p = this->halos[halo];
  dataSum = data[p];
  dataRem = 0.0;
  
  // Next particle
  p = this->haloList[p];
  
  // Remaining steps in Kahan summation
  while (p != -1) {
    v = data[p] - dataRem;
    w = dataSum + v;
    dataRem = (w - dataSum) - v;
    dataSum = w;

    p = this->haloList[p];
  }
  return dataSum;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate the incremental mean using Kahan summation
//
/////////////////////////////////////////////////////////////////////////

POSVEL_T HaloProperties::incrementalMean(int halo, POSVEL_T* data)
{
  double dataMean, dataRem, diff, value, v, w;

  // First particle in halo and first step in incremental mean
  int p = this->halos[halo];
  dataMean = data[p];
  dataRem = 0.0;
  int count = 1;

  // Next particle
  p = this->haloList[p];
  count++;

  // Remaining steps in incremental mean
  while (p != -1) {
    diff = data[p] - dataMean;
    value = diff / count;
    v = value - dataRem;
    w = dataMean + v;
    dataRem = (w - dataMean) - v;
    dataMean = w;

    p = this->haloList[p];
    count++;
  }
  return (POSVEL_T) dataMean;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate the mimimum potential or halo center using (N*(N-1)) / 2 algorithm
//
/////////////////////////////////////////////////////////////////////////

int HaloProperties::minimumPotentialN2_2(int halo, POTENTIAL_T* minPotential)
{
  int x = this->halos[halo];
  int cnt = 0;
  while (x != -1) {
    cnt++;
    x = this->haloList[x];
  }

  // Arrange in an upper triangular grid to save computation
  POTENTIAL_T* pot = new POTENTIAL_T[cnt];
  int* actualIndx = new int[cnt];
  x = this->halos[halo];
  for (int i = 0; i < cnt; i++) {
    pot[i] = 0.0;
    actualIndx[i] = x;
    x = this->haloList[x];
  }

  // First particle in halo to calculate minimum potential on
  int p = this->halos[halo];
  *minPotential = MAX_FLOAT;
  int minIndex = this->halos[halo];
  int count = 0;
  int indx1 = 0;

  while (p != -1) {

    // Next particle in halo in minimum potential loop
    int q = this->haloList[p];
    POTENTIAL_T potential = 0.0;
    int indx2 = indx1 + 1;
    count++;

    while (q != -1) {

      POSVEL_T xdist = fabs(this->xx[p] - this->xx[q]);
      POSVEL_T ydist = fabs(this->yy[p] - this->yy[q]);
      POSVEL_T zdist = fabs(this->zz[p] - this->zz[q]);

      POSVEL_T xval = min(xdist, this->boxSize - xdist);
      POSVEL_T yval = min(ydist, this->boxSize - ydist);
      POSVEL_T zval = min(zdist, this->boxSize - zdist);

      POSVEL_T r = sqrt((xval * xval) + (yval * yval) + (zval * zval));
      if (r != 0.0) {
        pot[indx1] = pot[indx1] - (1.0 / r);
        pot[indx2] = pot[indx2] - (1.0 / r);
      }

      // Next particle
      q = this->haloList[q];
      indx2++;
    }

    // Next particle
    p = this->haloList[p];
    indx1++;
  }

  for (int i = 0; i < cnt; i++) {
    if (pot[i] < *minPotential) {
      *minPotential = pot[i];
      minIndex = i;
    }
  } 
  delete [] pot;
  delete [] actualIndx;

  return actualIndx[minIndex];
}

/////////////////////////////////////////////////////////////////////////
//
// Create the chaining mesh which organizes particles into location grids
// by creating buckets of locations and chaining the indices of the
// particles so that all particles in a bucket can be located
//
/////////////////////////////////////////////////////////////////////////

void HaloProperties::createChainingMesh(
			POSVEL_T boxSz,
			POSVEL_T deadSz,
			POSVEL_T chainSz)

{
  this->boxSize = boxSz;
  this->deadSize = deadSz;
  this->chainSize = chainSz;

  // Calculate the physical boundary on this processor for alive particles
  POSVEL_T boxStep[DIMENSION];
  POSVEL_T minAlive[DIMENSION];
  POSVEL_T maxAlive[DIMENSION];

  for (int dim = 0; dim < DIMENSION; dim++) {
    boxStep[dim] = this->boxSize / this->layoutSize[dim];
   
    // Region of particles that are alive on this processor
    minAlive[dim] = this->layoutPos[dim] * boxStep[dim];
    maxAlive[dim] = minAlive[dim] + boxStep[dim];
    if (maxAlive[dim] > this->boxSize)
      maxAlive[dim] = this->boxSize;
      
    // Allow for the boundary of dead particles, normalized to 0
    // Overall boundary will be [0:(rL+2*deadSize)]
    this->minMine[dim] = minAlive[dim];
    this->maxMine[dim] = maxAlive[dim] + (2.0 * this->deadSize);
  }   

  // How many chain mesh grids will fit
  for (int dim = 0; dim < DIMENSION; dim++) {
    this->meshSize[dim] = (int) (((this->maxMine[dim] - this->minMine[dim]) / 
                            this->chainSize)) + 1;
  }

  // Create the bucket grid and initialize to -1
  this->buckets = new int**[this->meshSize[0]];
  this->bucketCount = new int**[this->meshSize[0]];

  for (int i = 0; i < this->meshSize[0]; i++) {
    this->buckets[i] = new int*[this->meshSize[1]];
    this->bucketCount[i] = new int*[this->meshSize[1]];

    for (int j = 0; j < this->meshSize[1]; j++) {
      this->buckets[i][j] = new int[this->meshSize[2]];
      this->bucketCount[i][j] = new int[this->meshSize[2]];

      for (int k = 0; k < this->meshSize[2]; k++) {
        this->buckets[i][j][k] = -1;
        this->bucketCount[i][j][k] = 0;
      }
    }
  }

  // Create the chaining list of particles and initialize to -1
  this->bucketList = new int[this->particleCount];
  for (int p = 0; p < this->particleCount; p++)
    this->bucketList[p] = -1;

  // Iterate on all particles on this processor and assign a bucket grid
  // First particle index is assigned to the actual bucket grid
  // Next particle found is assigned to the bucket grid, only after the
  // index which is already there is assigned to the new particles index
  // position in the bucketList.
  // Then to iterate through all particles in a bucket, start with the index
  // in the buckets grid and follow it throught the bucketList until -1 reached

  for (int p = 0; p < this->particleCount; p++) {

    // Normalize to zero
    POSVEL_T loc[DIMENSION];
    loc[0] = this->xx[p] + this->deadSize;
    loc[1] = this->yy[p] + this->deadSize;
    loc[2] = this->zz[p] + this->deadSize;

    int i = (int) ((loc[0] - this->minMine[0]) / this->chainSize);
    int j = (int) ((loc[1] - this->minMine[1]) / this->chainSize);
    int k = (int) ((loc[2] - this->minMine[2]) / this->chainSize);

    // First particle in bucket
    if (this->buckets[i][j][k] == -1) {
      this->buckets[i][j][k] = p;
      this->bucketCount[i][j][k]++;
    }

    // Other particles in same bucket
    else {
      this->bucketList[p] = this->buckets[i][j][k];
      this->buckets[i][j][k] = p;
      this->bucketCount[i][j][k]++;
    }
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Demonstration method to show how to iterate over the chaining mesh
// Calculate the centroid of each bucket
//
/////////////////////////////////////////////////////////////////////////

void HaloProperties::chainingMeshCentroids()
{
  // Test by calculating centroid of each bucket grid
  POSVEL_T centroid[DIMENSION];

  // Iterate on all particles in a bucket
  for (int i = 0; i < this->meshSize[0]; i++) {
    for (int j = 0; j < this->meshSize[1]; j++) {
      for (int k = 0; k < this->meshSize[2]; k++) {
        centroid[0] = 0.0;
        centroid[1] = 0.0;
        centroid[2] = 0.0;
    
        // First particle in the bucket
        int p = this->buckets[i][j][k];

        while (p != -1) {
          centroid[0] += this->xx[p];
          centroid[1] += this->yy[p];
          centroid[2] += this->zz[p];

          // Next particle in the bucket
          p = this->bucketList[p];
        }
        for (int dim = 0; dim < DIMENSION; dim++) {
          if (centroid[dim] != 0.0)
            centroid[dim] /= this->bucketCount[i][j][k];
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Calculations involving spherical over density (SOD)
// Estimate characteristic radius
// Mass profile
// Re-estimate characteristic radius
// Velocity dispersion
//
/////////////////////////////////////////////////////////////////////////

void HaloProperties::calculateSOD(
			int minHaloSize,
			double particleMass,
			vector<int>* haloCenter)
{
  // Iterate on all halos
  for (int halo = 0; halo < this->numberOfHalos; halo++) {

    if (this->haloCount[halo] >= minHaloSize) {

      this->sphereParticles.clear();

      double haloMass = this->haloCount[halo] * particleMass;
      double cubeRoot = 1.0 / 3.0;
      double estRadius = pow((haloMass / (1.0e14 * particleMass)), cubeRoot);

      double maxRadius = 6.0;

      // Index of the particle at the center of this halo
// PKF - Note that I don't know this yet because I need the minimum potential
// for all particles to find the center.  So for testing, just choose the
// first particle in the halo

      //int particleCenter = (*haloCenter)[halo];
      int particleCenter = this->halos[halo];

      // Radii go out logarithmically from the center particle
      double radius = maxRadius;
      collectSODParticles(particleCenter, radius);

    }
  }
}

/////////////////////////////////////////////////////////////////////////
//
// For the given radius and particle center collect the indices of all
// particles within the described sphere using the bucketList
//
/////////////////////////////////////////////////////////////////////////

void HaloProperties::collectSODParticles(
			int particleCenter,
			double radius)
{
  // Grid index containing the center particle, remembering to normalize 0
  POSVEL_T center[DIMENSION], location[DIMENSION];
  center[0] = this->xx[particleCenter] + this->deadSize;
  center[1] = this->yy[particleCenter] + this->deadSize;
  center[2] = this->zz[particleCenter] + this->deadSize;

  int centerIndex[DIMENSION];
  for (int dim = 0; dim < DIMENSION; dim++)
    centerIndex[dim] = (int) ((center[dim] - this->minMine[dim]) / 
                               this->chainSize);

  // How do I decide what grids to look at
  // Then take the location of each particle in those grids and find the
  // distance to the center.  If < radius, then collect it

  // Number of grids to look at in each direction
  int gridOffset = (int) ((radius / this->chainSize)) + 1;

  // Range of grid positions to examine for this particle center
  int first[DIMENSION], last[DIMENSION];
  for (int dim = 0; dim < DIMENSION; dim++) {
    first[dim] = centerIndex[dim] - gridOffset;
    last[dim] = centerIndex[dim] + gridOffset;
    if (first[dim] < 0)
      first[dim] = 0;
    if (last[dim] > this->meshSize[dim])
      last[dim] = this->meshSize[dim];
  }

  // Iterate over every possible grid and collect particles if any
  // Probably can look at the closest point in the grid to the center
  // to exclude an entire grid
  for (int i = first[0]; i < last[0]; i++) {
    for (int j = first[1]; j < last[1]; j++) {
      for (int k = first[2]; k < last[2]; k++) {

        // Iterate on all particles in this bucket
        // Index of first particle in bucket
        int p = this->buckets[i][j][k];
        while (p != -1) {

          location[0] = this->xx[p] + this->deadSize;
          location[1] = this->yy[p] + this->deadSize;
          location[2] = this->zz[p] + this->deadSize;

          // Calculate distance between this particle and the center
          double diff[DIMENSION];
          for (int dim = 0; dim < DIMENSION; dim++)
            diff[dim] = fabs(location[dim] - center[dim]);

          double distance = sqrt((diff[0] * diff[0]) +
                                 (diff[1] * diff[1]) +
                                 (diff[2] * diff[2]));
          // If the distance is less than the radius record this index
          if (distance < radius)
            sphereParticles.push_back(p);

          // Next particle in bucket
          p = this->bucketList[p];
        }
      }
    }
  }
}
