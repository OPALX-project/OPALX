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

#include "Partition.h"
#include "HaloCenterFinder.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////
//
// HaloCenterFinder uses the results of the CosmoHaloFinder to locate the
// particle within every halo with the minimum potential
//
/////////////////////////////////////////////////////////////////////////

HaloCenterFinder::HaloCenterFinder(
			vector<int>& haloStartIndex,
			vector<int>& haloParticleCount,
			int* nextParticleIndex) :
				halos(haloStartIndex),
				haloCount(haloParticleCount),
				haloList(nextParticleIndex)
{
  // Get the number of processors and rank of this processor
  this->numProc = Partition::getNumProc();
  this->myProc = Partition::getMyProc();

  // Get the number of processors in each dimension
  Partition::getDecompSize(this->layoutSize);

  // Get my position within the Cartesian topology
  Partition::getMyPosition(this->layoutPos);
}

HaloCenterFinder::~HaloCenterFinder()
{
}

/////////////////////////////////////////////////////////////////////////
//
// Set parameters for the halo center finder
//
/////////////////////////////////////////////////////////////////////////

void HaloCenterFinder::setParameters(
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

void HaloCenterFinder::setParticles(
                        vector<POSVEL_T>* xLoc,
                        vector<POSVEL_T>* yLoc,
                        vector<POSVEL_T>* zLoc,
                        vector<POSVEL_T>* xVel,
                        vector<POSVEL_T>* yVel,
                        vector<POSVEL_T>* zVel,
                        vector<ID_T>* id,
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
  this->tag = &(*id)[0];
  this->status = &(*state)[0];
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate the minimum potential or halo center using using a combination
// of the N^2 algorithm for small halos and the A* algorithm for large
//
/////////////////////////////////////////////////////////////////////////

void HaloCenterFinder::minimumPotential()
{
  int smallHalo = 0;
  int largeHalo = 0;
  for (int h = 0; h < this->halos.size(); h++) {
    if (this->haloCount[h] <= 100) {
      smallHalo++;
      minimumPotentialN2_2(h);
    } else {
      largeHalo++;
      //minimumPotentialAStar(h);
    }
  }
  cout << "Rank " << this->myProc 
       << " #small = " << smallHalo
       << " #large = " << largeHalo << endl;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate the mimimum potential or halo center using N^2 algorithm
//
/////////////////////////////////////////////////////////////////////////

void HaloCenterFinder::minimumPotentialN2(int h)
{
  // First particle in halo to calculate minimum potential on
  int p = this->halos[h];
  POSVEL_T minPotential = MAX_FLOAT;
  int minIndex = this->halos[h];
  int count = 0;

  while (p != -1) {

    // First particle in halo in minimum potential loop
    int q = this->halos[h];
    POSVEL_T potential = 0.0;
    count++;

    while (q != -1) {

      POSVEL_T xdist = fabs(this->xx[p] - this->xx[q]);
      POSVEL_T ydist = fabs(this->yy[p] - this->yy[q]);
      POSVEL_T zdist = fabs(this->zz[p] - this->zz[q]);

      POSVEL_T xval = min(xdist, this->boxSize - xdist);
      POSVEL_T yval = min(ydist, this->boxSize - ydist);
      POSVEL_T zval = min(zdist, this->boxSize - zdist);

      POSVEL_T r = sqrt((xval * xval) + (yval * yval) + (zval * zval));
      if (r != 0.0)
        potential = potential - (1.0 / r);

      // Next particle
      q = this->haloList[q];
    }

    // Finished calculating potential for this halo, is it the least so far
    if (potential < minPotential) {
      minPotential = potential;
      minIndex = p;
    }

    // Next particle
    p = this->haloList[p];
  }

  // What is the center for this halo
/*
  cout << "Halo " << h << " Center is " << this->xx[minIndex] << "   "
                      << this->yy[minIndex] << "    "
                      << this->zz[minIndex] << " Size " << count
                      << "  Potential " << minPotential << endl;
*/
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate the mimimum potential or halo center using (N*(N-1)) / 2 algorithm
//
/////////////////////////////////////////////////////////////////////////

void HaloCenterFinder::minimumPotentialN2_2(int h)
{
  int x = this->halos[h];
  int cnt = 0;
  while (x != -1) {
    cnt++;
    x = this->haloList[x];
  }
  POSVEL_T* pot = new POSVEL_T[cnt];
  int* actualIndx = new int[cnt];
  x = this->halos[h];
  for (int i = 0; i < cnt; i++) {
    pot[i] = 0.0;
    actualIndx[i] = x;
    x = this->haloList[x];
  }

  // First particle in halo to calculate minimum potential on
  int p = this->halos[h];
  POSVEL_T minPotential = MAX_FLOAT;
  int minIndex = this->halos[h];
  int count = 0;
  int indx1 = 0;

  while (p != -1) {

    // Next particle in halo in minimum potential loop
    int q = this->haloList[p];
    POSVEL_T potential = 0.0;
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
    if (pot[i] < minPotential) {
      minPotential = pot[i];
      minIndex = i;
    }
  } 

  // What is the center for this halo
/*
  cout << "Rank " << myProc << " Halo " << h << " Center is " 
                  << this->xx[actualIndx[minIndex]] << "   "
                  << this->yy[actualIndx[minIndex]] << "    "
                  << this->zz[actualIndx[minIndex]] << " Size " << count 
                  << "  Potential " << minPotential << endl;
*/
  delete [] pot;
  delete [] actualIndx;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate the minimum potential using the A* search algorithm
//
/////////////////////////////////////////////////////////////////////////

void HaloCenterFinder::minimumPotentialAStar(int h)
{
  // Save the minimum potential for particles in a halo
  // Allocate space for all particles for now
  POSVEL_T* minPot = new POSVEL_T[this->particleCount];

  // Processor must know its alive particle physical range in order to
  // handle periodic conditions
  POSVEL_T minDead[DIMENSION];
  POSVEL_T maxDead[DIMENSION];
  POSVEL_T boxStep[DIMENSION];

  for (int dim = 0; dim < DIMENSION; dim++) {
    boxStep[dim] = this->boxSize / this->layoutSize[dim];

    minDead[dim] = (this->layoutPos[dim] * boxStep[dim]) - this->deadSize;
    maxDead[dim] = minDead[dim] + boxStep[dim] + this->deadSize;
    if (maxDead[dim] > (this->boxSize + this->deadSize))
      maxDead[dim] = this->boxSize + this->deadSize;
  }

  // Iterate over all particles in this halo to find the minimum bounding box
  // PKF Must deal with the periodic condition
  // How do we know a particle meets the periodic condition?
  int p = this->halos[h];
  POSVEL_T minBound[DIMENSION], maxBound[DIMENSION];
  for (int dim = 0; dim < DIMENSION; dim++) {
    minBound[dim] = MAX_FLOAT;
    maxBound[dim] = MIN_FLOAT;
  }

  // Iterate on particles
  while (p != -1) {
    POSVEL_T xLoc = this->xx[p];
    POSVEL_T yLoc = this->yy[p];
    POSVEL_T zLoc = this->zz[p];

    // Make adjustments to physical location depending on periodic condition
    if (this->xx[p] < minDead[0]) xLoc += this->boxSize;
    if (this->xx[p] > maxDead[0]) xLoc -= this->boxSize;
    if (this->yy[p] < minDead[1]) yLoc += this->boxSize;
    if (this->yy[p] > maxDead[1]) yLoc -= this->boxSize;
    if (this->zz[p] < minDead[2]) zLoc += this->boxSize;
    if (this->zz[p] > maxDead[2]) zLoc -= this->boxSize;

    if (minBound[0] > xLoc) minBound[0] = xLoc;
    if (minBound[1] > yLoc) minBound[1] = yLoc;
    if (minBound[2] > zLoc) minBound[2] = zLoc;
    if (maxBound[0] < xLoc) maxBound[0] = xLoc;
    if (maxBound[1] < yLoc) maxBound[1] = yLoc;
    if (maxBound[2] < zLoc) maxBound[2] = zLoc;
    p = this->haloList[p];
  }

  // Impose a grid over the bounding box using the given physical step
  int gridSize[DIMENSION];
  POSVEL_T step = 1.0;
  int numberOfGrids = 1;
  for (int dim = 0; dim < DIMENSION; dim++) {
    minBound[dim] = floorf(minBound[dim]);
    maxBound[dim] = ceilf(maxBound[dim]);
    gridSize[dim] = (int) ((maxBound[dim] - minBound[dim]) / step) + 1;
    numberOfGrids * gridSize[dim];
  }

  // Create the grid of vectors to hold particles indices which have been
  // assigned to the physical grid
  vector<int>* grid = new vector<int>[numberOfGrids];

  // Iterate over all particles and assign to grid positions
  p = this->halos[h];
  while (p != -1) {
    POSVEL_T xLoc = this->xx[p];
    POSVEL_T yLoc = this->yy[p];
    POSVEL_T zLoc = this->zz[p];

    // Make adjustments to physical location depending on periodic condition
    if (this->xx[p] < minDead[0]) xLoc += this->boxSize;
    if (this->xx[p] > maxDead[0]) xLoc -= this->boxSize;
    if (this->yy[p] < minDead[1]) yLoc += this->boxSize;
    if (this->yy[p] > maxDead[1]) yLoc -= this->boxSize;
    if (this->zz[p] < minDead[2]) zLoc += this->boxSize;
    if (this->zz[p] > maxDead[2]) zLoc -= this->boxSize;

    // Place in the correct vector
    int iindx = (int) ((xLoc - minBound[0]) / step);
    int jindx = (int) ((yLoc - minBound[1]) / step);
    int kindx = (int) ((zLoc - minBound[2]) / step);

    int index = kindx * (gridSize[0] * gridSize[1]) +
                jindx * (gridSize[0]) +
                iindx;
cout << "Grid size " << gridSize[0] << "," << gridSize[1] << "," << gridSize[2] << "  index = " << index << endl;
    grid[index].push_back(p);
    p = this->haloList[p];
  }

  // For each particle in the halo calculate the g(x) which is the actual
  // minimum potential between that particle and every other particle in
  // the 27 neighbor block surrounding it.

  // This is being done on a single processor so can we use a sliding window
  // to save on calculations?  Also how do we handle periodic condition?

  // For every particle in an entire grid position we can calculate the
  // h(x) estimate by looking at all other grids in the system and taking
  // the count of particles in that grid times the nearest corner

  // Now every particle has a g(x) actual and h(x) estimate which gives
  // f(x) which is what A* will use

  // Loop 
  // A* selects a best minimum potential particle
  // For that particle we calculate a new g(x) by walking out one more
  // level and getting better actual values, and must also subtract from
  // h(x) the former estimate

  // Go through the loop again and let A* select another best particle and
  // refine for it  (every refinement would have to catch up with the current
  // level of refinement so if we have a new best point that still has
  // the level one refinement and we have done level 5, then we must look
  // 5 levels out

  cout << "Rank " << myProc << " Halo " << h 
                  << " X bound " << minBound[0] << " to " << maxBound[0]
                  << " Y bound " << minBound[1] << " to " << maxBound[1]
                  << " Z bound " << minBound[2] << " to " << maxBound[2]
                  << " COUNT " << this->haloCount[h] << endl;
  delete [] grid;
  delete [] minPot;
}
