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

#include "CosmoHaloFinderP.h"

//#ifdef MC4HALOFINDER
#include "../../ChargedParticles/ChargedParticles.hh"
//#else
//#include "Partition.h"
//#endif

using namespace std;

/////////////////////////////////////////////////////////////////////////
//
// Parallel manager for serial CosmoHaloFinder
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
// The serial halo finder is called on each processor and returns enough
// information so that it can be determined if a halo is completely ALIVE,
// completely DEAD, or mixed.  A mixed halo that is shared between on two
// processors is kept by the processor that contains it in one of its
// high plane neighbors, and is given up if contained in a low plane neighbor.
//
// Mixed halos that cross more than two processors are bundled up and sent
// to the MASTER processor which decides the processor that should own it.
//
/////////////////////////////////////////////////////////////////////////


#ifdef MC4HALOFINDER
CosmoHaloFinderP::CosmoHaloFinderP(ChargedParticles<TT,3> *myU ) 
{
 // Get the number of processors and rank of this processor
    this->numProc = Ippl::getNodes();
    this->myProc = Ippl::myNode();

    // For each neighbor zone, how many dead particles does it contain to start
    // and how many dead halos does it contain after the serial halo finder
    // For analysis but not necessary to run the code
    //
    for (int n = 0; n < NUM_OF_NEIGHBORS; n++) {
	this->deadParticle[n] = 0;
	this->deadHalo[n] = 0;
    }
}  
#else
CosmoHaloFinderP::CosmoHaloFinderP()
{
  // Get the number of processors and rank of this processor
    this->numProc = Partition::getNumProc();
    this->myProc = Partition::getMyProc();

    // Get the number of processors in each dimension
    Partition::getDecompSize(this->layoutSize);
  
    // Get my position within the Cartesian topology
    Partition::getMyPosition(this->layoutPos);

    // Get the neighbors of this processor
    Partition::getNeighbors(this->neighbor);
    // For each neighbor zone, how many dead particles does it contain to start
    // and how many dead halos does it contain after the serial halo finder
    // For analysis but not necessary to run the code
    //
    for (int n = 0; n < NUM_OF_NEIGHBORS; n++) {
	this->deadParticle[n] = 0;
	this->deadHalo[n] = 0;
    }
}
#endif

CosmoHaloFinderP::~CosmoHaloFinderP()
{
  for (int i = 0; i < this->myMixedHalos.size(); i++)
    delete this->myMixedHalos[i];

  delete [] this->haloList;
  delete [] this->haloStart;
  delete [] this->haloSize;

  for (int dim = 0; dim < DIMENSION; dim++)
    delete haloData[dim];
  delete [] haloData;
}

/////////////////////////////////////////////////////////////////////////
//
// Set parameters for the serial halo finder
//
/////////////////////////////////////////////////////////////////////////

void CosmoHaloFinderP::setParameters(
			const string& outName,
			POSVEL_T rL,
			POSVEL_T deadSz,
			long np,
			int pmin,
			POSVEL_T bb)
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
  this->outHaloFile = hname.str();

  // Halo finder parameters
  this->np = np;
  this->pmin = pmin;
  this->bb = bb;
  this->boxSize = rL;
  this->deadSize = deadSz;

  // First version of this code distributed the dead particles on a processor
  // by taking the x,y,z position and adding or subtracting boxSize in all
  // combinations.  This revised x,y,z was then normalized and added to the
  // haloData array which was passed to each serial halo finder.
  // This did not get the same answer as the standalone serial version which
  // read the x,y,z and normalized without adding or subtracting boxSize first.
  // Then when comparing distance the normalized "np" was used for subtraction.
  // By doing things in this order some particles were placed slightly off,
  // which was enough for particles to be included in halos where they should
  // not have been.  In this first version, since I had placed particles by
  // subtracting first, I made "periodic" false figuring all particles were
  // placed where they should go.
  //
  // In the second version I normalize the dead particles, even from wraparound,
  // using the actual x,y,z.  So when looking at a processor the alive particles
  // will appear all together and the wraparound will properly be on the other
  // side of the box.  Combined with doing this is setting "periodic" to true
  // so that the serial halo finder works as it does in the standalone case
  // and the normalization and subtraction from np happens in the same order.
  //
  // Third version went back to the first version because we need 
  // contiguous locations coming out of the halo finder for the center finder

  this->haloFinder.np = np;
  this->haloFinder.pmin = pmin;
  this->haloFinder.bb = bb;
  this->haloFinder.rL = rL;
  this->haloFinder.periodic = false;
  this->haloFinder.textmode = "ascii";

  // Serial halo finder wants normalized locations on a grid superimposed
  // on the physical rL grid.  Grid size is np and number of particles in
  // the problem is np^3
  this->normalizeFactor = rL / (1.0 * np);

  if (this->myProc == MASTER) {
    cout << endl << "------------------------------------" << endl;
    cout << "np:       " << this->np << endl;
    cout << "bb:       " << this->bb << endl;
    cout << "pmin:     " << this->pmin << endl << endl;
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Set the particle vectors that have already been read and which
// contain only the alive particles for this processor
//
/////////////////////////////////////////////////////////////////////////
#ifdef MC4HALOFINDER
void CosmoHaloFinderP::setParticles(ChargedParticles<TT,3> *myU)				
{
    this->particleCount = myU->getLocalNum() + myU->getGhostNum();
    cout << "I have " << myU->getGhostNum() << " ghost particles" << endl;
    
    vector<POSVEL_T>* xx = new vector<POSVEL_T>;
    vector<POSVEL_T>* yy = new vector<POSVEL_T>;
    vector<POSVEL_T>* zz = new vector<POSVEL_T>;
    vector<POSVEL_T>* vx = new vector<POSVEL_T>;
    vector<POSVEL_T>* vy = new vector<POSVEL_T>;
    vector<POSVEL_T>* vz = new vector<POSVEL_T>;
    vector<ID_T>* tag = new vector<ID_T>;
    vector<STATUS_T>* status = new vector<STATUS_T>;
   
    for (int i=0; i<this->particleCount; i++) {
        xx->push_back((float)myU->R[i](0));
        yy->push_back((float)myU->R[i](1));
        zz->push_back((float)myU->R[i](2));
        vx->push_back((float)myU->V[i](0));
        vy->push_back((float)myU->V[i](1));
        vz->push_back((float)myU->V[i](2));
        tag->push_back(myU->ID[i]);
        status->push_back(-1);
     }

    /**
       Now find out the particles needed
       from neighboring nodes and add them.
   
    */

    size_t nn = myU->getNeighborInformation(this->bb);
    
    this->xx = &(*xx)[0];
    this->yy = &(*yy)[0];
    this->zz = &(*zz)[0];
    this->vx = &(*vx)[0];
    this->vy = &(*vy)[0];
    this->vz = &(*vz)[0];
    this->tag = &(*tag)[0];

    this->status = &(*status)[0];

    // not used
    this->pot = NULL;
    this->mask = NULL;
}
#else
void CosmoHaloFinderP::setParticles(
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
// Execute the serial halo finder on the particles for this processor
// Both ALIVE and DEAD particles we collected and normalized into
// haloData which is in the form that the serial halo finder wants
//
/////////////////////////////////////////////////////////////////////////

void CosmoHaloFinderP::executeHaloFinder()
{
  // Allocate data pointer which is sent to the serial halo finder
  this->haloData = new POSVEL_T*[DIMENSION];
  for (int dim = 0; dim < DIMENSION; dim++)
    this->haloData[dim] = new POSVEL_T[this->particleCount];

  // Fill it with normalized x,y,z of all particles on this processor
  for (int p = 0; p < this->particleCount; p++) {

    this->haloData[0][p] = this->xx[p] / this->normalizeFactor;
    this->haloData[1][p] = this->yy[p] / this->normalizeFactor;
    this->haloData[2][p] = this->zz[p] / this->normalizeFactor;
  }

  this->haloFinder.setParticleLocations(haloData);
  this->haloFinder.setNumberOfParticles(this->particleCount);
  this->haloFinder.setMyProc(this->myProc);
  this->haloFinder.setOutFile(this->outFile);

  cout << "Rank " << setw(3) << this->myProc 
       << " RUNNING SERIAL HALO FINDER on " 
       << particleCount << " particles" << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  if (this->particleCount > 0)
    this->haloFinder.Finding();
  MPI_Barrier(MPI_COMM_WORLD);
}

/////////////////////////////////////////////////////////////////////////
//
// At this point each serial halo finder ran and
// the particles handed to it included alive and dead.  Get back the 
// halo tag array and figure out the indices of the particles in each halo
// and translate that into absolute particle tags and note alive or dead
//
// After the serial halo finder has run the halo tag is the INDEX of the
// lowest particle in the halo on this processor.  It is not the absolute
// particle tag id over the entire problem.
//
//    Serial partindex i = 0 haloTag = 0 haloSize = 1
//    Serial partindex i = 1 haloTag = 1 haloSize = 1
//    Serial partindex i = 2 haloTag = 2 haloSize = 1
//    Serial partindex i = 3 haloTag = 3 haloSize = 1
//    Serial partindex i = 4 haloTag = 4 haloSize = 2
//    Serial partindex i = 5 haloTag = 5 haloSize = 1
//    Serial partindex i = 6 haloTag = 6 haloSize = 1616
//    Serial partindex i = 7 haloTag = 7 haloSize = 1
//    Serial partindex i = 8 haloTag = 8 haloSize = 2
//    Serial partindex i = 9 haloTag = 9 haloSize = 1738
//    Serial partindex i = 10 haloTag = 10 haloSize = 4
//    Serial partindex i = 11 haloTag = 11 haloSize = 1
//    Serial partindex i = 12 haloTag = 12 haloSize = 78
//    Serial partindex i = 13 haloTag = 12 haloSize = 0
//    Serial partindex i = 14 haloTag = 12 haloSize = 0
//    Serial partindex i = 15 haloTag = 12 haloSize = 0
//    Serial partindex i = 16 haloTag = 16 haloSize = 2
//    Serial partindex i = 17 haloTag = 17 haloSize = 1
//    Serial partindex i = 18 haloTag = 6 haloSize = 0
//    Serial partindex i = 19 haloTag = 6 haloSize = 0
//    Serial partindex i = 20 haloTag = 6 haloSize = 0
//    Serial partindex i = 21 haloTag = 6 haloSize = 0
//
// Halo of size 1616 has the low particle tag of 6 and other members are
// 18,19,20,21 indicated by a tag of 6 and haloSize of 0
//
/////////////////////////////////////////////////////////////////////////

void CosmoHaloFinderP::collectHalos()
{
  // Halo tag returned from the serial halo finder is actually the index
  // of the particle on this processor.  Must map to get to actual tag
  // which is common information between all processors.
  this->haloTag = haloFinder.getHaloTag();

  // Record the halo size of each particle on this processor
  this->haloSize = new int[this->particleCount];
  this->haloAliveSize = new int[this->particleCount];
  this->haloDeadSize = new int[this->particleCount];

  // Create a list of particles in any halo by recording the index of the
  // first particle and having that index give the index of the next particle
  // Last particle index reports a -1
  // List is built by iterating on the tags and storing in opposite order so
  this->haloList = new int[this->particleCount];
  this->haloStart = new int[this->particleCount];

  for (int p = 0; p < this->particleCount; p++) {
    this->haloList[p] = -1;
    this->haloStart[p] = p;
    this->haloSize[p] = 0;
    this->haloAliveSize[p] = 0;
    this->haloDeadSize[p] = 0;
  }

  // Build the chaining mesh of particles in all the halos and count particles
  buildHaloStructure();

  // Mixed halos are saved separately so that they can be merged
  processMixedHalos();

  delete [] this->haloAliveSize;
  delete [] this->haloDeadSize;
}

/////////////////////////////////////////////////////////////////////////
//
// Examine every particle on this processor, both ALIVE and DEAD
// For that particle increment the count for the corresponding halo
// which is indicated by the lowest particle index in that halo
// Also build the haloList so that we can find all particles in any halo
//
/////////////////////////////////////////////////////////////////////////

void CosmoHaloFinderP::buildHaloStructure()
{
  // Build the chaining mesh so that all particles in a halo can be found
  // This will include even small halos which will be excluded later
  for (int p = 0; p < this->particleCount; p++) {

    // Chain backwards the halo particles
    // haloStart is the index of the last particle in a single halo in haloList
    // The value found in haloList is the index of the next particle
    if (this->haloTag[p] != p) {
      this->haloList[p] = haloStart[this->haloTag[p]];
      this->haloStart[this->haloTag[p]] = p;
    }

    // Count particles in the halos
    if (this->status[p] == ALIVE)
      this->haloAliveSize[this->haloTag[p]]++;
    else
      this->haloDeadSize[this->haloTag[p]]++;
    this->haloSize[this->haloTag[p]]++;
  }

  // Iterate over particles and create a CosmoHalo for halos with size > pmin
  // only for the mixed halos, not for those completely alive or dead
  this->numberOfAliveHalos = 0;
  this->numberOfDeadHalos = 0;
  this->numberOfMixedHalos = 0;

  // Only the first particle id for a halo records the size
  // Succeeding particles which are members of a halo have a size of 0
  // Record the start index of any legal halo which will allow the
  // following of the chaining mesh to identify all particles in a halo
  this->numberOfHaloParticles = 0;
  for (ID_T p = 0; p < this->particleCount; p++) {

    if (this->haloSize[p] >= this->pmin) {

      if (this->haloAliveSize[p] > 0 && this->haloDeadSize[p] == 0) {
        this->numberOfAliveHalos++;
        this->numberOfHaloParticles += this->haloAliveSize[p];

        // Save start of legal alive halo for center finding
        this->halos.push_back(this->haloStart[p]);
        this->haloCount.push_back(this->haloAliveSize[p]);
      }
      else if (this->haloDeadSize[p] > 0 && this->haloAliveSize[p] == 0) {
        this->numberOfDeadHalos++;
      }
      else {
        this->numberOfMixedHalos++;
        CosmoHalo* halo = new CosmoHalo(p, 
                                this->haloAliveSize[p], this->haloDeadSize[p]);
        this->myMixedHalos.push_back(halo);

        // Save start of legal mixed halo for center finding
        this->halos.push_back(this->haloStart[p]);
        this->haloCount.push_back(this->haloAliveSize[p]+this->haloDeadSize[p]);
      }
    }
  }
#ifdef DEBUG
  cout << "Rank " << this->myProc 
       << " #alive halos = " << this->numberOfAliveHalos
       << " #dead halos = " << this->numberOfDeadHalos
       << " #mixed halos = " << this->numberOfMixedHalos << endl;
#endif
}

/////////////////////////////////////////////////////////////////////////
//
// Mixed halos (which cross several processors) have been collected
// By applying a high/low rule most mixed halos are assigned immediately
// to one processor or another.  This requires extra processing so that
// it is known which neighbor processors share the halo.
//
/////////////////////////////////////////////////////////////////////////

void CosmoHaloFinderP::processMixedHalos()
{
  // Iterate over all particles and add tags to large mixed halos 
  for (ID_T p = 0; p < this->particleCount; p++) {

    // All particles in the same halo have the same haloTag
    if (this->haloSize[this->haloTag[p]] >= pmin && 
        this->haloAliveSize[this->haloTag[p]] > 0 && 
        this->haloDeadSize[this->haloTag[p]] > 0) {

          // Check all each mixed halo to see which this particle belongs to
          for (int h = 0; h < this->myMixedHalos.size(); h++) {

            // If the tag of the particle matches the halo ID it belongs
            if (this->haloTag[p] == this->myMixedHalos[h]->getHaloID()) {

              // Add the index to that mixed halo.  Also record which neighbor
              // the dead particle is associated with for merging
              this->myMixedHalos[h]->addParticle(
                                        p, this->tag[p], this->status[p]);

              // For debugging only
              if (this->status[p] > 0)
                this->deadHalo[this->status[p]]++;

              // Do some bookkeeping for the final output
              // This processor should output all ALIVE particles, unless they
              // are in a mixed halo that ends up being INVALID
              // This processor should output none of the DEAD particles,
              // unless they are in a mixed halo that ends up being VALID

              // So since this particle is in a mixed halo set it to MIXED
              // which is going to be one less than ALIVE.  Later when we
              // determine we have a VALID mixed, we'll add one to the status
              // for every particle turning all into ALIVE

              // Now when we output we only do the ALIVE particles
              this->status[p] = MIXED;
            }
          }
    }
  }

  // Iterate over the mixed halos that were just created checking to see if
  // the halo is on the "high" side of the 3D data space or not
  // If it is on the high side and is shared with one other processor, keep it
  // If it is on the low side and is shared with one other processor, delete it
  // Any remaining halos are shared with more than two processors and must
  // be merged by having the MASTER node decide
  //
  for (int h = 0; h < this->myMixedHalos.size(); h++) {
    int lowCount = 0;
    int highCount = 0;
    set<int>* neighbors = this->myMixedHalos[h]->getNeighbors();
    set<int>::iterator iter;
    set<int> haloNeighbor;

    for (iter = neighbors->begin(); iter != neighbors->end(); ++iter) {
      if ((*iter) == X1 || (*iter) == Y1 || (*iter) == Z1 || 
          (*iter) == X1_Y1 || (*iter) == Y1_Z1 || (*iter) == Z1_X1 ||
          (*iter) == X1_Y1_Z1) {
            highCount++;
      } else {
            lowCount++;
      }
      // Neighbor zones are on what actual processors
      haloNeighbor.insert(this->neighbor[(*iter)]);
    }

    // Halo is kept by this processor and is marked as VALID
    // May be in multiple neighbor zones, but all the same processor neighbor
    if (highCount > 0 && lowCount == 0 && haloNeighbor.size() == 1) {
      this->numberOfAliveHalos++;
      this->numberOfMixedHalos--;
      this->myMixedHalos[h]->setValid(VALID);
      this->numberOfHaloParticles += this->myMixedHalos[h]->getAliveCount();
      this->numberOfHaloParticles += this->myMixedHalos[h]->getDeadCount();

      // Output trick - since the status of this particle was marked MIXED
      // when it was added to the mixed CosmoHalo vector, and now it has
      // been declared VALID, change it to ALIVE even if it was dead before
      vector<ID_T>* particles = this->myMixedHalos[h]->getParticles();
      vector<ID_T>::iterator iter;
      for (iter = particles->begin(); iter != particles->end(); ++iter)
        this->status[(*iter)] = ALIVE;
    }

    // Halo will be kept by some other processor and is marked INVALID
    // May be in multiple neighbor zones, but all the same processor neighbor
    else if (highCount == 0 && lowCount > 0 && haloNeighbor.size() == 1) {
      this->numberOfDeadHalos++;
      this->numberOfMixedHalos--;
      this->myMixedHalos[h]->setValid(INVALID);
    }

    // Remaining mixed halos must be examined by MASTER and stay UNMARKED
    // Sort them on the tag field for easy comparison
    else {
      this->myMixedHalos[h]->setValid(UNMARKED);
      this->myMixedHalos[h]->sortParticleTags();
    }
  }

  // If only one processor is running there are no halos to merge
  if (this->numProc == 1)
    for (int h = 0; h < this->myMixedHalos.size(); h++)
       this->myMixedHalos[h]->setValid(INVALID);
}

/////////////////////////////////////////////////////////////////////////
//
// Using the MASTER node merge all mixed halos so that only one processor
// takes credit for them.
//
// Each processor containing mixed halos that are UNMARKED sends:
//    Rank
//    Number of mixed halos to merge
//    for each halo
//      id
//      number of alive (for debugging)
//      number of dead  (for debugging)
//      first MERGE_COUNT particle ids (for merging)
//
/////////////////////////////////////////////////////////////////////////

void CosmoHaloFinderP::mergeHalos()
{
  // What size integer buffer is needed to hold the largest halo data
  int maxNumberOfMixed;
  int numberOfMixed = this->myMixedHalos.size();
  MPI_Allreduce((void*) &numberOfMixed, (void*) &maxNumberOfMixed,
                1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // If there are no halos to merge, return
  if (maxNumberOfMixed == 0)
    return;

  // Everyone creates the buffer for maximum halos
  // MASTER will receive into it, others will send from it
  int haloBufSize = maxNumberOfMixed * MERGE_COUNT * 2;
  ID_T* haloBuffer = new ID_T[haloBufSize];

  // MASTER moves its own mixed halos to mixed halo vector (change index to tag)
  // then gets messages from others and creates those mixed halos
  collectMixedHalos(haloBuffer, haloBufSize);
  MPI_Barrier(MPI_COMM_WORLD);

  // MASTER has all data and runs algorithm to make decisions
  assignMixedHalos();
  MPI_Barrier(MPI_COMM_WORLD);

  // MASTER sends merge results to all processors
  sendMixedHaloResults(haloBuffer, haloBufSize);
  MPI_Barrier(MPI_COMM_WORLD);

  // Collect totals for result checking
  int totalAliveHalos;
  MPI_Allreduce((void*) &this->numberOfAliveHalos, (void*) &totalAliveHalos,
                1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int totalAliveHaloParticles;
  MPI_Allreduce((void*) &this->numberOfHaloParticles,
                (void*) &totalAliveHaloParticles,
                1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (this->myProc == MASTER) {
    cout << endl;
    cout << "Total halos found:    " << totalAliveHalos << endl;
    cout << "Total halo particles: " << totalAliveHaloParticles << endl;
  }

  for (int i = 0; i < this->allMixedHalos.size(); i++)
    delete this->allMixedHalos[i];
  delete [] haloBuffer;
}

/////////////////////////////////////////////////////////////////////////
//
// MASTER collects all mixed halos which are UNMARKED from all processors
// including its own mixed halos
//
/////////////////////////////////////////////////////////////////////////

void CosmoHaloFinderP::collectMixedHalos(ID_T* haloBuffer, int haloBufSize)
{
  // How many processors have mixed halos
  int haveMixedHalo = (this->numberOfMixedHalos > 0 ? 1 : 0);
  int processorsWithMixedHalos;
  MPI_Allreduce((void*) &haveMixedHalo, (void*) &processorsWithMixedHalos,
                1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // MASTER moves its own mixed halos to mixed halo vector (change index to tag)
  // then gets messages from others and creates those mixed halos
  if (this->myProc == MASTER) {

    // If MASTER has any mixed halos add them to the mixed halo vector
    if (this->numberOfMixedHalos > 0) {
      processorsWithMixedHalos--;

      for (int h = 0; h < this->myMixedHalos.size(); h++) {
        if (this->myMixedHalos[h]->getValid() == UNMARKED) {
          CosmoHalo* halo = new CosmoHalo(
                                     this->myMixedHalos[h]->getHaloID(),
                                     this->myMixedHalos[h]->getAliveCount(),
                                     this->myMixedHalos[h]->getDeadCount());
          halo->setRankID(this->myProc);
          this->allMixedHalos.push_back(halo);

          // Translate index of particle to tag of particle
          vector<ID_T>* tags = this->myMixedHalos[h]->getTags();
          for (int i = 0; i < MERGE_COUNT; i++)
            halo->addParticle((*tags)[i]);

        }
      }
    }

    // Wait on messages from other processors and process
    int notReceived = processorsWithMixedHalos;
    MPI_Status status;
    while (notReceived > 0) {

      // Get message containing mixed halo information
#ifdef ID_64
      MPI_Recv(haloBuffer, haloBufSize, MPI_LONG, MPI_ANY_SOURCE,
               0, MPI_COMM_WORLD, &status);
#else
      MPI_Recv(haloBuffer, haloBufSize, MPI_INT, MPI_ANY_SOURCE,
               0, MPI_COMM_WORLD, &status);
#endif

      // Gather halo information from the message
      int index = 0;
      int rank = haloBuffer[index++];
      int numMixed = haloBuffer[index++];

      for (int m = 0; m < numMixed; m++) {
        ID_T id = haloBuffer[index++];
        int aliveCount = haloBuffer[index++];
        int deadCount = haloBuffer[index++];

        // Create the CosmoHalo to hold the data and add to vector
        CosmoHalo* halo = new CosmoHalo(id, aliveCount, deadCount);

        halo->setRankID(rank);
        this->allMixedHalos.push_back(halo);

        for (int t = 0; t < MERGE_COUNT; t++)
          halo->addParticle(haloBuffer[index++]);
      }
      notReceived--;
    }
    cout << "Number of halos to merge: " << this->allMixedHalos.size() << endl;
  }

  // Other processors bundle up mixed and send to MASTER
  else {
    int index = 0;
    if (this->numberOfMixedHalos > 0) {

      haloBuffer[index++] = this->myProc;
      haloBuffer[index++] = this->numberOfMixedHalos;

      for (int h = 0; h < this->myMixedHalos.size(); h++) {
        if (this->myMixedHalos[h]->getValid() == UNMARKED) {

          haloBuffer[index++] = this->myMixedHalos[h]->getHaloID();
          haloBuffer[index++] = this->myMixedHalos[h]->getAliveCount();
          haloBuffer[index++] = this->myMixedHalos[h]->getDeadCount();

          vector<ID_T>* tags = this->myMixedHalos[h]->getTags();
          for (int i = 0; i < MERGE_COUNT; i++) {
            haloBuffer[index++] = (*tags)[i];
          }
        }
      }
      MPI_Request request;
#ifdef ID_64
      MPI_Isend(haloBuffer, haloBufSize, MPI_LONG, MASTER,
                0, MPI_COMM_WORLD, &request);
#else
      MPI_Isend(haloBuffer, haloBufSize, MPI_INT, MASTER,
                0, MPI_COMM_WORLD, &request);
#endif
    }
  }
}

/////////////////////////////////////////////////////////////////////////
//
// MASTER has collected all the mixed halos and decides which processors
// will get which by matching them up
//
/////////////////////////////////////////////////////////////////////////

void CosmoHaloFinderP::assignMixedHalos()
{
  // MASTER has all data and runs algorithm to make decisions
  if (this->myProc == MASTER) {

#ifdef DEBUG
    for (int m = 0; m < this->allMixedHalos.size(); m++) {
      vector<ID_T>* tags = this->allMixedHalos[m]->getTags();
      cout << "Mixed Halo " << m << ": "
           << " rank=" << this->allMixedHalos[m]->getRankID()
           << " index=" << this->allMixedHalos[m]->getHaloID()
           << " tag=" << (*tags)[0]
           << " alive=" << this->allMixedHalos[m]->getAliveCount()
           << " dead=" << this->allMixedHalos[m]->getDeadCount() << endl; 
    }
#endif

    // Iterate over mixed halo vector and match and mark
    // Remember that I can have 3 or 4 that match
    for (int m = 0; m < this->allMixedHalos.size(); m++) {

      // If this halo has not already been paired with another
      if (this->allMixedHalos[m]->getPartners()->empty() == true) {

        // Iterate on the rest of the mixed halos
        int n = m + 1;
        while (n < this->allMixedHalos.size()) {

          // Compare to see if there are a number of tags in common
          int match = compareHalos(this->allMixedHalos[m], 
                                   this->allMixedHalos[n]);
          if (match > 0) {
            this->allMixedHalos[m]->addPartner(n);
            this->allMixedHalos[n]->addPartner(-1);
            this->allMixedHalos[m]->setValid(VALID);
            this->allMixedHalos[n]->setValid(INVALID);
          }
          n++;
        }
      }
    }

#ifdef DEBUG
    for (int m = 0; m < this->allMixedHalos.size(); m++) {
      cout << "Mixed Halo " << m << " partners with ";
      set<int>* partner = this->allMixedHalos[m]->getPartners();
      set<int>::iterator iter;
      for (iter = partner->begin(); iter != partner->end(); ++iter)
        cout << (*iter) << " ";
      cout << endl;
    }
#endif
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Compare the tags of two halos to see if they are somewhat the same
// This needs to be made better
//
/////////////////////////////////////////////////////////////////////////

int CosmoHaloFinderP::compareHalos(CosmoHalo* halo1, CosmoHalo* halo2)
{
  vector<ID_T>* member1 = halo1->getTags();
  vector<ID_T>* member2 = halo2->getTags();

  int numFound = 0;
  for (int i = 0; i < member1->size(); i++) {
    bool done = false;
    int j = 0;
    while (!done && 
           (*member1)[i] >= (*member2)[j] && 
           j < member2->size()) {
      if ((*member1)[i] == (*member2)[j]) {
        done = true;
        numFound++;
      }
      j++;
    }
  }
  if (numFound == 0)
    return numFound;
  else
    return numFound;
}

/////////////////////////////////////////////////////////////////////////
//
// MASTER sends the result of the merge back to the processors which
// label their previously UNMARKED mixed halos as VALID or INVALID
// VALID halos have all their particles made ALIVE for output
// INVALID halos have all their particles made DEAD because other
// processors will report them
//
/////////////////////////////////////////////////////////////////////////

void CosmoHaloFinderP::sendMixedHaloResults(ID_T* haloBuffer, int haloBufSize)
{
  // MASTER sends merge results to all processors
  if (this->myProc == MASTER) {

    // Share the information
    // Send to each processor the rank, id, and valid status
    // Use the same haloBuffer
    int index = 0;
    haloBuffer[index++] = this->allMixedHalos.size();
    for (int m = 0; m < this->allMixedHalos.size(); m++) {
      haloBuffer[index++] = this->allMixedHalos[m]->getRankID();
      haloBuffer[index++] = this->allMixedHalos[m]->getHaloID();
      haloBuffer[index++] = this->allMixedHalos[m]->getValid();
    }

    MPI_Request request;
    for (int proc = 1; proc < this->numProc; proc++) {
#ifdef ID_64
      MPI_Isend(haloBuffer, haloBufSize, MPI_LONG, proc,
                0, MPI_COMM_WORLD, &request); 
#else
      MPI_Isend(haloBuffer, haloBufSize, MPI_INT, proc,
                0, MPI_COMM_WORLD, &request); 
#endif
    }

    // MASTER must claim the mixed halos assigned to him
    for (int m = 0; m < this->allMixedHalos.size(); m++) {
      if (this->allMixedHalos[m]->getRankID() == MASTER &&
          this->allMixedHalos[m]->getValid() == VALID) {

        // Locate the mixed halo in question
        for (int h = 0; h < this->myMixedHalos.size(); h++) {
          if (this->myMixedHalos[h]->getHaloID() == 
              this->allMixedHalos[m]->getHaloID()) {
            this->myMixedHalos[h]->setValid(VALID);
            this->numberOfHaloParticles += 
                           this->myMixedHalos[h]->getAliveCount();
            this->numberOfHaloParticles += 
                           this->myMixedHalos[h]->getDeadCount();
            this->numberOfAliveHalos++;

            // Output trick - since the status of this particle was marked MIXED
            // when it was added to the mixed CosmoHalo vector, and now it has
            // been declared VALID, change it to ALIVE even if it was dead
            vector<ID_T>* particles = this->myMixedHalos[h]->getParticles();
            vector<ID_T>::iterator iter;
            for (iter = particles->begin(); iter != particles->end(); ++iter)
              this->status[(*iter)] = ALIVE;
          }
        }
      }
    }
  }

  // Other processors wait for result and adjust their halo vector
  else {
    MPI_Status status;
#ifdef ID_64
    MPI_Recv(haloBuffer, haloBufSize, MPI_LONG, MASTER,
             0, MPI_COMM_WORLD, &status);
#else
    MPI_Recv(haloBuffer, haloBufSize, MPI_INT, MASTER,
             0, MPI_COMM_WORLD, &status);
#endif

    // Unpack information to see which of mixed halos are still valid
    int index = 0;
    int numMixed = haloBuffer[index++];
    for (int m = 0; m < numMixed; m++) {
      int rank = haloBuffer[index++];
      int id = haloBuffer[index++];
      int valid = haloBuffer[index++];

      // If this mixed halo is on my processor
      if (rank == this->myProc && valid == VALID) {

        // Locate the mixed halo in question
        for (int h = 0; h < this->myMixedHalos.size(); h++) {
          if (this->myMixedHalos[h]->getHaloID() == id) {
            this->myMixedHalos[h]->setValid(VALID);
            this->numberOfHaloParticles += 
                           this->myMixedHalos[h]->getAliveCount();
            this->numberOfHaloParticles += 
                           this->myMixedHalos[h]->getDeadCount();
            this->numberOfAliveHalos++;

            // Output trick - since the status of this particle was marked MIXED
            // when it was added to the mixed CosmoHalo vector, and now it has
            // been declared VALID, change it to ALIVE even if it was dead
            vector<ID_T>* particles = this->myMixedHalos[h]->getParticles();
            vector<ID_T>::iterator iter;
            for (iter = particles->begin(); iter != particles->end(); ++iter)
              this->status[(*iter)] = ALIVE;
          }
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Write the output of the halo finder in the form of the input .cosmo file
//
// Encoded mixed halo VALID or INVALID into the status array such that
// ALIVE particles that are part of an INVALID mixed array will not write
// but DEAD particles that are part of a VALID mixed array will be written
//
// In order to make the output consistent with the serial output where the
// lowest tagged particle in a halo owns the halo, work must be done to
// identify the lowest tag.  This is because as particles are read onto
// this processor using the round robin read of every particle, those
// particles are no longer in tag order.  When the serial halo finder is
// called it has to use the index of the particle on this processor which
// is no longer the tag.
//
//      p    haloTag     tag    haloSize
//      0          0     523           3
//      1          0     522           0
//      2          0     266           0
//
// In the above example the halo will be credited to 523 instead of 266
// because the index of 523 is 0 and the index of 266 is 2.  So we must
// make a pass to map the indexes.
//
/////////////////////////////////////////////////////////////////////////

void CosmoHaloFinderP::writeHalos()
{
  // Map the index of the particle on this process to the index of the
  // particle with the lowest tag value so that the written output refers
  // to the lowest tag as being the owner of the halo
  int* mapIndex = new int[this->particleCount];
  for (int p = 0; p < this->particleCount; p++)
    mapIndex[p] = p;

  // If the tag for the first particle of this halo is bigger than the tag
  // for this particle, change the map to identify this particle as the lowest
  for (int p = 0; p < this->particleCount; p++) {
    if (this->tag[mapIndex[this->haloTag[p]]] > this->tag[p])
      mapIndex[this->haloTag[p]] = p;
  }

  // Write the halo catalog
  ofstream* outStream = new ofstream(this->outFile.c_str(), ios::out);

  // Write halo information of halo tag followed by size for histgram
  ofstream* haloStream = new ofstream(this->outHaloFile.c_str(), ios::out);

  string textMode = "ascii";
  char str[1024];
  string ascii = "ascii";
  if (textMode == "ascii") {

    // Output all ALIVE particles that were not part of a mixed halo 
    // unless that halo is VALID.  Output only the DEAD particles that are
    // part of a VALID halo. This was encoded when mixed halos were found
    // so any ALIVE particle is VALID

    for (int p = 0; p < this->particleCount; p++) {

      if (this->status[p] == ALIVE) {
        // Every alive particle appears in the particle output
        sprintf(str, "%12.4E %12.4E ", this->xx[p], this->vx[p]);
        *outStream << str;
        sprintf(str, "%12.4E %12.4E ", this->yy[p], this->vy[p]);
        *outStream << str;
        sprintf(str, "%12.4E %12.4E ", this->zz[p], this->vz[p]);
        *outStream << str;
        int result = (this->haloSize[this->haloTag[p]] < this->pmin) 
                      ? -1: this->tag[mapIndex[this->haloTag[p]]];
        sprintf(str, "%12d %12d\n", result, this->tag[p]);
        *outStream << str;

        // Only ALIVE particles which are the first particle in a halo
        // get written to the halo file with the halo size
        if (this->haloSize[p] >= this->pmin) {
          *haloStream << this->tag[mapIndex[p]] << "\t"
                      << this->haloSize[p] << endl;
        }
      }
    }
  }

  else {

    // output in COSMO form
    for (int p = 0; p < this->particleCount; p++) {
      float fBlock[COSMO_FLOAT];
      fBlock[0] = this->xx[p];
      fBlock[1] = this->vx[p];
      fBlock[2] = this->yy[p];
      fBlock[3] = this->vy[p];
      fBlock[4] = this->zz[p];
      fBlock[5] = this->vz[p];
      fBlock[6] = (this->haloSize[this->haloTag[p]] < this->pmin) 
                   ? -1.0: 1.0 * this->tag[this->haloTag[p]];
      outStream->write((char *)fBlock, COSMO_FLOAT * sizeof(float));

      int iBlock[COSMO_INT];
      iBlock[0] = this->tag[p];
      outStream->write((char *)iBlock, COSMO_INT * sizeof(int));
    }
  }
  outStream->close();
  haloStream->close();

  delete outStream;
  delete haloStream;
  delete [] mapIndex;
}
