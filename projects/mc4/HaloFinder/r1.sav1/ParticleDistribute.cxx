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
#include "ParticleDistribute.h"

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

ParticleDistribute::ParticleDistribute()
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

  this->xx = new vector<float>;
  this->yy = new vector<float>;
  this->zz = new vector<float>;
  this->vx = new vector<float>;
  this->vy = new vector<float>;
  this->vz = new vector<float>;
  this->tag = new vector<int>;
  this->status = new vector<int>;
}

ParticleDistribute::~ParticleDistribute()
{
  delete this->xx;
  delete this->yy;
  delete this->zz;
  delete this->vx;
  delete this->vy;
  delete this->vz;
  delete this->tag;
  delete this->status;
}

/////////////////////////////////////////////////////////////////////////
//
// Set parameters for particle distribution
//
/////////////////////////////////////////////////////////////////////////

void ParticleDistribute::setParameters(
			const string& baseName,
			float rL,
			string dataType)
{
  // Base file name which will have processor id appended for actual files
  this->baseFile = baseName;

  // Physical total space and amount of physical space to use for dead particles
  this->boxSize = rL;

  // RECORD format is the binary .cosmo of one particle with all information
  if (dataType == "RECORD")
    this->inputType = RECORD;

  // BLOCK format is Gadget format with a header and x,y,z locations for
  // all particles, then x,y,z velocities for all particles, and all tags
  else if (dataType == "BLOCK")
    this->inputType = BLOCK;

  if (this->myProc == MASTER) {
    cout << endl << "------------------------------------" << endl;
    cout << "boxSize:  " << this->boxSize << endl;
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Set box sizes for determining if a particle is in the alive or dead
// region of this processor.  Data space is a DIMENSION torus.
//
/////////////////////////////////////////////////////////////////////////

void ParticleDistribute::initialize()
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

    // Alive particles
    this->minAlive[dim] = this->layoutPos[dim] * boxStep[dim];
    this->maxAlive[dim] = this->minAlive[dim] + boxStep[dim];
    if (this->maxAlive[dim] > this->boxSize)
      this->maxAlive[dim] = this->boxSize;
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Each processor reads 0 or more files, a buffer at a time, and shares
// the particles by passing the buffer round robin to every other processor
//
/////////////////////////////////////////////////////////////////////////

void ParticleDistribute::readParticlesRoundRobin()
{
  // Find how many input files there are and deal them between the processors
  // Calculates the max number of files per processor and max number of
  // particles per file so that buffering can be done
  // For round robin sharing determine where to send and receive buffers from
  partitionInputFiles();

  // Compute the total number of particles in the problem
  // Compute the maximum number of particles in any one file to set buffer size
  findFileParticleCount();

  // MPI buffer size might limit the number of particles read from a file
  // and passed round robin
  // Largest file will have a number of buffer chunks to send if it is too large
  // Every processor must send that number of chunks even if its own file
  // does not have that much information

  if (this->maxParticles > MAX_READ) {
    this->maxRead = MAX_READ;
    this->maxReadsPerFile = (this->maxParticles / this->maxRead) + 1;
  } else {
    this->maxRead = this->maxParticles;
    this->maxReadsPerFile = 1;
  }

  // Allocate space to hold buffer information for reading of files
  // Mass is constant use that float to store the tag
  // Number of particles is the first integer in the buffer
  int bufferSize = 1 + (this->maxRead * COSMO_FLOAT);
  float* buffer1 = new float[bufferSize];
  float* buffer2 = new float[bufferSize];

  // Allocate space for the data read from the file
  float* fBlock;
  int* iBlock;

  // RECORD format reads one particle at a time
  if (this->inputType == RECORD) {
    fBlock = new float[COSMO_FLOAT];
    iBlock = new int[COSMO_INT];
  }

  // BLOCK format reads all particles at one time for triples
  else if (this->inputType == BLOCK) {
    fBlock = new float[this->maxRead * 3];
    iBlock = new int[this->maxRead];
  }

  // Reserve particle storage to minimize reallocation
  int reserveSize = (int) (this->maxFiles * this->maxParticles * DEAD_FACTOR);

  // If multiple processors are reading the same file we can reduce size
  reserveSize /= this->processorsPerFile;

  this->xx->reserve(reserveSize);
  this->yy->reserve(reserveSize);
  this->zz->reserve(reserveSize);
  this->vx->reserve(reserveSize);
  this->vy->reserve(reserveSize);
  this->vz->reserve(reserveSize);
  this->tag->reserve(reserveSize);
  this->status->reserve(reserveSize);

  // Running total and index into particle data on this processor
  this->particleCount = 0;

  // Using the input files assigned to this processor, read the input
  // and push round robin to every other processor
  // this->maxFiles is the maximum number to read on any processor
  // Some processors may have no files to read but must still participate 
  // in the round robin distribution

  for (int file = 0; file < this->maxFiles; file++) {

    // Open file to read the data if any for this processor
    ifstream *inStream;
    int firstParticle = 0;
    int numberOfParticles = 0;
    int remainingParticles = 0;

    if (this->inFiles.size() > file) {
      inStream = new ifstream(this->inFiles[file].c_str(), ios::in);

      cout << "Rank " << this->myProc << " open file " << inFiles[file]
           << " with " << this->fileParticles[file] << " particles" << endl;

      // First particle to read in BLOCK format is 1
      firstParticle = 0;
      if (this->inputType == BLOCK)
        firstParticle = 1;

      // Number of particles read at one time depends on MPI buffer size
      numberOfParticles = this->fileParticles[file];
      if (numberOfParticles > this->maxRead)
        numberOfParticles = this->maxRead;

      // If a file is too large to be passed as an MPI message divide it up
      remainingParticles = this->fileParticles[file];

    } else {
      cout << "Rank " << this->myProc << " no file to open " << endl;
    }

    for (int piece = 0; piece < this->maxReadsPerFile; piece++) {

      // Processor has a file to read and share via round robin with others
      if (file < this->inFiles.size()) {
        if (this->inputType == RECORD) {
          readFromRecordFile(inStream, firstParticle, numberOfParticles,
                             fBlock, iBlock, buffer1);
        } else {
          readFromBlockFile(inStream, firstParticle, numberOfParticles,
                            this->fileParticles[file], fBlock, iBlock, buffer1);
        }
        firstParticle += numberOfParticles;
        remainingParticles -= numberOfParticles;
        if (remainingParticles <= 0)
          numberOfParticles = 0;
        else if (remainingParticles < numberOfParticles)
          numberOfParticles = remainingParticles;
      }

      // Processor does not have a file to open but must participate in the
      // round robin with an empty buffer
      else {
        // Store number of particles used in first position
        int numberOfParticles = 0;
        buffer1[0] = *(reinterpret_cast<float*>(&numberOfParticles));
      }

      // Particles belonging to this processor are put in vectors
      distributeParticles(fBlock, iBlock, buffer1, buffer2, bufferSize);
    }

    // Can delete the read buffers as soon as last file is read because
    // information has been transferred into the double buffer1
    if (file == (this->maxFiles - 1)) {
      delete [] fBlock;
      delete [] iBlock;
    }

#ifdef DEBUG
    cout << "Rank " << setw(3) << this->myProc << " FINISHED READ " << endl;
#endif

    if (this->inFiles.size() > file)
      inStream->close();
  }

  // After all particles have been distributed to vectors the double
  // buffers can be deleted
  delete [] buffer1;
  delete [] buffer2;

  // Count the particles across processors
  int totalAliveParticles = 0;
  MPI_Allreduce((void*) &this->numberOfAliveParticles, 
                (void*) &totalAliveParticles, 
                1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#ifdef DEBUG
  cout << "Rank " << setw(3) << this->myProc 
       << " #alive = " << this->numberOfAliveParticles << endl;
#endif
 
  if (this->myProc == MASTER) {
    cout << "TotalAliveParticles " << totalAliveParticles << endl;
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Using the base name of the data, go to the subdirectory and determine
// how many input files there are.  Parcel those files between all the
// processors which will be responsible for actually reading 0 or more.
//
/////////////////////////////////////////////////////////////////////////

void ParticleDistribute::partitionInputFiles()
{
  // Find number of input files for this problem given the base input name
  // Get the subdirectory name containing the input files
  string::size_type dirPos = this->baseFile.rfind("/");
  string subdirectory;
  string baseName;

  // If the directory is not given use the current directory
  // Get the base name of input files of the form filename."number" 
  if (dirPos == string::npos) {
    subdirectory = "./";
    baseName = this->baseFile + '.';
  } else {
    subdirectory = this->baseFile.substr(0, dirPos);
    baseName = this->baseFile.substr(dirPos+1);
    baseName = baseName + '.';
  }

  // Open the subdirectory and make a list of input files
  DIR* directory = opendir(subdirectory.c_str());
  struct dirent* directoryEntry;
  vector<string> files;

  if (directory != NULL) {
    while (directoryEntry = readdir(directory)) {

      // Only files containing the base name with a '.##' are recognized
      string fileName = directoryEntry->d_name;
      string::size_type pos = fileName.find(baseName.c_str());
      string::size_type pos1 = pos + baseName.size();

      if (pos != string::npos &&
          (fileName[pos1] >= '0' && fileName[pos1] <= '9')) {
        string inputFile = subdirectory + '/' + fileName;
        files.push_back(inputFile);
      }
    }
  }
  this->numberOfFiles = files.size();
  if (this->numberOfFiles == 0) {
    cout << "Rank " << this->myProc << " found no input files" << endl;
    exit(1);
  }

#ifdef DEBUG
  if (this->myProc == MASTER) {
    for (int i = 0; i < this->numberOfFiles; i++)
      cout << "   File " << i << ": " << files[i] << endl;
  }
#endif

  // Divide the files between all the processors
  // If there is at least one file per processor or if there are fewer but
  // the number of files does not divide the number of processors set the
  // buffering up with a full round robin between all processors

  if (this->numberOfFiles >= this->numProc || 
      ((this->numProc % this->numberOfFiles) != 0)) {

    // Number of round robin sends to share all the files
    this->processorsPerFile = 1;
    this->numberOfFileSends = this->numProc - 1;

    // Which files does this processor read
    for (int i = 0; i < this->numberOfFiles; i++)
      if ((i % this->numProc) == this->myProc)
        this->inFiles.push_back(files[i]);

    // Where is the file sent, and where is it received
    if (this->myProc == this->numProc - 1)
      this->nextProc = 0;
    else
      this->nextProc = this->myProc + 1;
    if (this->myProc == 0)
      this->prevProc = this->numProc - 1;
    else
      this->prevProc = this->myProc - 1;
  }

  // If fewer files than processors and files evenly divide into processors
  // set up smaller round robin circles
  else {

    // Number of round robin sends to share all the files
    this->processorsPerFile = this->numProc / this->numberOfFiles;
    this->numberOfFileSends = this->numberOfFiles - 1;

    // What file does this processor read (several procs will read same file)
    for (int i = 0; i < this->numberOfFiles; i++)
      if ((this->myProc % this->numberOfFiles) == i)
        this->inFiles.push_back(files[i]);

    // Where is the file sent, and where is it received
    if (((this->myProc + this->numberOfFiles) % this->numberOfFiles) == 
                                                  this->numberOfFiles - 1)
      this->nextProc = this->myProc - this->numberOfFiles + 1;
    else
      this->nextProc = this->myProc + 1;
    if (((this->myProc + this->numberOfFiles) % this->numberOfFiles) == 0)
      this->prevProc = this->myProc + this->numberOfFiles - 1;
    else
      this->prevProc = this->myProc - 1;
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Open each input file belonging to this processor and find the number
// of particles for setting buffer sizes
//
/////////////////////////////////////////////////////////////////////////

void ParticleDistribute::findFileParticleCount()
{
  // Compute the total number of particles in the problem
  // Compute the maximum number of particles in any one file to set buffer size
  int numberOfParticles = 0;
  int maxNumberOfParticles = 0;
  int numberOfMyFiles = this->inFiles.size();

  // Each processor counts the particles in its own files
  for (int i = 0; i < numberOfMyFiles; i++) {

    // Open my file
    ifstream *inStream = new ifstream(this->inFiles[i].c_str(), ios::in);
    if (inStream->fail()) {
      delete inStream;
      cout << "File: " << this->inFiles[i] << " cannot be opened" << endl;
      exit (-1);
    }

    if (this->inputType == RECORD) {

      // Compute the number of particles from file size
      inStream->seekg(0L, ios::end);
      int numberOfRecords = inStream->tellg() / RECORD_SIZE;
      this->fileParticles.push_back(numberOfRecords);

      numberOfParticles += numberOfRecords;
      if (maxNumberOfParticles < numberOfRecords)
        maxNumberOfParticles = numberOfRecords;
    }

    else if (this->inputType == BLOCK) {

      // Find the number of particles in the header
      inStream->read(reinterpret_cast<char*>(&this->cosmoHeader),
                     sizeof(this->cosmoHeader));

      this->headerSize = this->cosmoHeader.npart[0];
      if (sizeof(this->cosmoHeader) != this->headerSize)
        cout << "Mismatch of header size and header structure" << endl;

      int numberOfRecords = this->cosmoHeader.npart[2];
      this->fileParticles.push_back(numberOfRecords);

      numberOfParticles += numberOfRecords;
      if (maxNumberOfParticles < numberOfRecords)
        maxNumberOfParticles = numberOfRecords;
    }

    inStream->close();
    delete inStream;
  }

  // If multiple processors read the same file, just do the reduce on one set
  if (this->processorsPerFile > 1) {
    if (this->myProc >= this->numberOfFiles) {
      numberOfParticles = 0;
      maxNumberOfParticles = 0;
    }
  }

  // Share the information about total particles
  MPI_Allreduce((void*) &numberOfParticles,
                (void*) &this->totalParticles,
                1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // Share the information about max particles in a file for setting buffer size
  MPI_Allreduce((void*) &maxNumberOfParticles,
                (void*) &this->maxParticles,
                1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // Share the maximum number of files on a processor for setting the loop
  MPI_Allreduce((void*) &numberOfMyFiles,
                (void*) &this->maxFiles,
                1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

#ifdef DEBUG
  if (this->myProc == MASTER) {
    cout << "Total particle count: " << this->totalParticles << endl;
    cout << "Max particle count:   " << this->maxParticles << endl;
  }
#endif
}

/////////////////////////////////////////////////////////////////////////
//
// Each processor reads 0 or more files, a buffer at a time.
// The particles are processed by seeing if they are in the subextent of
// this processor and are tagged either ALIVE or if dead, by the index of
// the neighbor zone which contains that particle.  That buffer is sent
// round robin to (myProc + 1) % numProc where it is processed and sent on.
// After each processor reads one buffer and sends and receives numProc - 1
// times the next buffer from the file is read.  Must use a double buffering
// scheme so that on each send/recv we switch buffers.
//
// Input files may be BLOCK or RECORD structured
//
/////////////////////////////////////////////////////////////////////////

void ParticleDistribute::distributeParticles(
		float* fBlock,		// Read data from file into block
		int* iBlock,		// Read data from file into block
		float* buffer1,		// Double send/receive buffers
		float* buffer2,		// Double send/receive buffers
		int bufferSize)		// Size of the send/receive buffers
{
  // Each processor has filled a buffer with particles read from a file
  // or had no particles to read but set the count in the buffer to 0
  // Process the buffer to keep only those within range
  float* sendBuffer = buffer1;
  float* recvBuffer = buffer2;

  // Process the original send buffer of particles from the file
  collectLocalParticles(sendBuffer);

  // Distribute buffer round robin so that all processors see it
  for (int step = 0; step < this->numberOfFileSends; step++) {

    // Send buffer to the next processor
    MPI_Request request;
    MPI_Isend(sendBuffer, bufferSize, MPI_FLOAT, this->nextProc,
              0, MPI_COMM_WORLD, &request);

    // Receive buffer from the previous processor
    MPI_Status mpiStatus;
    MPI_Recv(recvBuffer, bufferSize, MPI_FLOAT, this->prevProc,
             0, MPI_COMM_WORLD, &mpiStatus);

    MPI_Barrier(MPI_COMM_WORLD);

    // Process the send buffer for alive and dead before sending on
    collectLocalParticles(recvBuffer);

#ifdef DEBUG
    if (this->myProc == MASTER)
      cout << "Round robin pass " << step+1 << endl;
#endif

    // Swap buffers for next round robin send
    float* tmp = sendBuffer;
    sendBuffer = recvBuffer;
    recvBuffer = tmp;
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Input file is RECORD structured so read each particle record and populate
// the double buffer in particle order for the rest of the processing
//
/////////////////////////////////////////////////////////////////////////////

void ParticleDistribute::readFromRecordFile(
			ifstream* inStream,	// Stream to read from
			int firstParticle,	// First particle index
			int numberOfParticles,	// Number to read this time
			float* fBlock,		// Buffer for read in data
			int* iBlock,		// Buffer for read in data
			float* buffer)		// Reordered data
{
  // Store number of particles used in first position
  buffer[0] = *(reinterpret_cast<float*>(&numberOfParticles));
  if (numberOfParticles == 0)
    return;

  // Seek to the first particle locations and read
  int floatSkip = COSMO_FLOAT * sizeof(float);
  int intSkip = COSMO_INT * sizeof(int);
  int skip = (floatSkip + intSkip) * firstParticle;
  inStream->seekg(skip, ios::beg);

  // Store each particle location, velocity in buffer and replace
  // constant mass by float tag id
  int index = 1;
  for (int i = 0; i < numberOfParticles; i++) {

    // Set file pointer to the requested particle
    inStream->read(reinterpret_cast<char*>(fBlock),
                   COSMO_FLOAT * sizeof(float));

    if (inStream->gcount() != COSMO_FLOAT * sizeof(float)) {
      cout << "Premature end-of-file" << endl;
      exit (-1);
    }

    inStream->read(reinterpret_cast<char*>(iBlock),
                   COSMO_INT * sizeof(int));
    if (inStream->gcount() != COSMO_INT * sizeof(int)) {
      cout << "Premature end-of-file" << endl;
      exit (-1);
    }

    // Store location and velocity in buffer but not the constant mass
    for (int i = 0; i < (COSMO_FLOAT - 1); i++) {
      buffer[index] = fBlock[i];
      index++;
    }

    // Store the integer tag in the float mass position reinterpreted
    buffer[index] = *(reinterpret_cast<float*>(&iBlock[0]));
    index++;
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Input file is BLOCK structured so read head and each block of data.
//    Header
//    Block of location data where each particle's x,y,z is stored together
//    Block of velocity data where each particle's xv,yv,zv is stored together
//    Block of tag data
// Reorder the data after it is read into the same structure as the
// RECORD data so that the rest of the code does not have to be changed
//
/////////////////////////////////////////////////////////////////////////////

void ParticleDistribute::readFromBlockFile(
			ifstream* inStream,	// Stream to read from
			int firstParticle,	// First particle index
			int numberOfParticles,	// Number to read this time
			int totParticles,	// Total particles in file
			float* fBlock,		// Buffer for read in data
			int* iBlock,		// Buffer for read in data
			float* buffer)		// Reordered data
{
  // Store number of particles used in first position
  buffer[0] = *(reinterpret_cast<float*>(&numberOfParticles));
  if (numberOfParticles == 0)
    return;

  // Seek to the first particle locations and read triples
  int skip = this->headerSize + (3 * sizeof(float) * firstParticle);
  inStream->seekg(skip, ios::beg);
  inStream->read(reinterpret_cast<char*>(fBlock),
                 3 * numberOfParticles * sizeof(float));

  // Store the locations in the double buffer in record order
  int findx = 0;
  int bindx = 1;
  for (int p = 0; p < numberOfParticles; p++) {
    buffer[bindx] = fBlock[findx];		// X value
    buffer[bindx+2] = fBlock[findx+1];		// Y value
    buffer[bindx+4] = fBlock[findx+2];		// Z value
    bindx += 7;
    findx += 3;
  }

  // Seek to first particle velocities and read triples
  skip = this->headerSize +                     // skip header
         (3 * sizeof(float) * totParticles) +   // skip all locations
         (3 * sizeof(float) * firstParticle);   // skip to first
  inStream->seekg(skip, ios::beg);
  inStream->read(reinterpret_cast<char*>(fBlock),
                 3 * numberOfParticles * sizeof(float));

  // Store the velocities in the double buffer in record order
  findx = 0;
  bindx = 2;
  for (int p = 0; p < numberOfParticles; p++) {
    buffer[bindx] = fBlock[findx];		// X value
    buffer[bindx+2] = fBlock[findx+1];		// Y value
    buffer[bindx+4] = fBlock[findx+2];		// Z value
    bindx += 7;
    findx += 3;
  }

  // Seek to first particle tags and read
  skip = this->headerSize +                     // skip header
         (3 * sizeof(float) * totParticles) +   // skip all locations
         (3 * sizeof(float) * totParticles) +   // skip all velocities
         (1 * sizeof(int) * firstParticle);     // skip to first tag
  inStream->seekg(skip, ios::beg);
  inStream->read(reinterpret_cast<char*>(iBlock),
                 1 * numberOfParticles * sizeof(int));

  // Store the tags in the double buffer for sharing in record order
  findx = 0;
  bindx = 7;
  for (int p = 0; p < numberOfParticles; p++) {

    // Store the integer tag reinterpreted as a float of same number of bytes
    // because tag can be large
    buffer[bindx] = *(reinterpret_cast<float*>(&iBlock[findx]));
    bindx += 7;
    findx += 1;
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Process the data buffer of particles to choose those which are ALIVE
// or DEAD on this processor.  Do wraparound tests to populate as for a
// 3D torus.  Dead particle status is the zone id of the neighbor processor
// which contains it as an ALIVE particle.
//
/////////////////////////////////////////////////////////////////////////////

void ParticleDistribute::collectLocalParticles(float* buf)
{
  int numParticles = *(reinterpret_cast<int*>(&buf[0]));
  float loc[DIMENSION];

  // Test each particle in the buffer to see if it is ALIVE or DEAD
  // If it is DEAD assign it to the neighbor zone that it is in
  // Check all combinations of wraparound

  int index = 1;
  for (int i = 0; i < numParticles; i++) {

    for (int dim = 0; dim < DIMENSION; dim++)
      loc[dim] = buf[index+(dim*2)];

    // Is the particle ALIVE on this processor
    if ((loc[0] >= minAlive[0] && loc[0] < maxAlive[0]) &&
        (loc[1] >= minAlive[1] && loc[1] < maxAlive[1]) &&
        (loc[2] >= minAlive[2] && loc[2] < maxAlive[2])) {

          this->xx->push_back(buf[index+0]);
          this->vx->push_back(buf[index+1]);
          this->yy->push_back(buf[index+2]);
          this->vy->push_back(buf[index+3]);
          this->zz->push_back(buf[index+4]);
          this->vz->push_back(buf[index+5]);
          this->tag->push_back(*(reinterpret_cast<int*>(&buf[index+6])));

          this->status->push_back(ALIVE);
          this->numberOfAliveParticles++;
          this->particleCount++;
    }
    index = index + SAVED_DATA;
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Each processor reads 1 file or gets a pointer to data eventually
// As the particle is read it will be stored as an alive particle on this
// processor and will be checked about neighbor ranges to see if it must
// be exchanged
//
/////////////////////////////////////////////////////////////////////////
  
void ParticleDistribute::readParticlesOneToOne()
{
  // File name is the base file name with processor id appended
  // Because MPI Cartesian topology is used the arrangement of files in
  // physical space must follow the rule of last dimension varies fastest
  ostringstream fileName;
  fileName << this->baseFile << "." << this->myProc;
  this->inFiles.push_back(fileName.str());

  // Compute the total number of particles in the problem
  // Compute the maximum number of particles in any one file to set buffer size
  findFileParticleCount();
  
  // Reserve particle storage to minimize reallocation
  int reserveSize = (int) (this->maxParticles * DEAD_FACTOR);
  
  this->xx->reserve(reserveSize);
  this->yy->reserve(reserveSize);
  this->zz->reserve(reserveSize);
  this->vx->reserve(reserveSize);
  this->vy->reserve(reserveSize);
  this->vz->reserve(reserveSize);
  this->tag->reserve(reserveSize);
  this->status->reserve(reserveSize);
  
  // Running total and index into particle data on this processor
  this->particleCount = 0;
  
  // Read the input file storing particles immediately because all are alive
  if (this->inputType == RECORD) {
    readFromRecordFile();
  } else {
    readFromBlockFile();
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Input file is RECORD structured so read each particle record and populate
// the vectors of particles marking all as ALIVE
//
/////////////////////////////////////////////////////////////////////////////

void ParticleDistribute::readFromRecordFile()
{
  // Only one file per processor named in index 0
  ifstream inStream(this->inFiles[0].c_str(), ios::in);
  int numberOfParticles = this->fileParticles[0];

  cout << "Rank " << this->myProc << " open file " << this->inFiles[0]
       << " with " << numberOfParticles << " particles" << endl;

  float* fBlock = new float[COSMO_FLOAT];
  int* iBlock = new int[COSMO_INT];

  // Store each particle location, velocity and tag
  int index = 1;
  for (int i = 0; i < numberOfParticles; i++) {

    // Set file pointer to the requested particle
    inStream.read(reinterpret_cast<char*>(fBlock),
                   COSMO_FLOAT * sizeof(float));

    if (inStream.gcount() != COSMO_FLOAT * sizeof(float)) {
      cout << "Premature end-of-file" << endl;
      exit (-1);
    }

    inStream.read(reinterpret_cast<char*>(iBlock),
                   COSMO_INT * sizeof(int));
    if (inStream.gcount() != COSMO_INT * sizeof(int)) {
      cout << "Premature end-of-file" << endl;
      exit (-1);
    }

    // Store location and velocity in buffer but not the constant mass
    this->xx->push_back(fBlock[0]);
    this->vx->push_back(fBlock[1]);
    this->yy->push_back(fBlock[2]);
    this->vy->push_back(fBlock[3]);
    this->zz->push_back(fBlock[4]);
    this->vz->push_back(fBlock[5]);
    this->tag->push_back(iBlock[0]);

    this->status->push_back(ALIVE);
    this->numberOfAliveParticles++;
    this->particleCount++;
  }
#ifdef DEBUG
    cout << "Rank " << setw(3) << this->myProc << " FINISHED READ " << endl;
#endif

  inStream.close();
  delete [] fBlock;
  delete [] iBlock;
}

/////////////////////////////////////////////////////////////////////////////
//
// Input file is BLOCK structured so read head and each block of data.
//    Header
//    Block of location data where each particle's x,y,z is stored together
//    Block of velocity data where each particle's xv,yv,zv is stored together
//    Block of tag data
//
/////////////////////////////////////////////////////////////////////////////

void ParticleDistribute::readFromBlockFile()
{
  // Only one file per processor named in index 0
  ifstream inStream(this->inFiles[0].c_str(), ios::in);
  int numberOfParticles = this->fileParticles[0];

  cout << "Rank " << this->myProc << " open file " << this->inFiles[0]
       << " with " << numberOfParticles << " particles" << endl;

  // Particle 0 in block data is not valid, tags start with 1 so
  // it must be skipped at the beginning of each block read
  int firstParticle = 1;

  // Allocate buffers for block reads
  float* fBlock = new float[numberOfParticles * 3];
  int* iBlock = new int[numberOfParticles];

  // Seek to the first particle locations and read triples
  int skip = this->headerSize + (3 * sizeof(float) * firstParticle);
  inStream.seekg(skip, ios::beg);
  inStream.read(reinterpret_cast<char*>(fBlock),
                3 * numberOfParticles * sizeof(float));

  // Store the locations in the double buffer in record order
  for (int p = 0; p < numberOfParticles; p++) {
    this->xx->push_back(fBlock[0]);
    this->yy->push_back(fBlock[1]);
    this->zz->push_back(fBlock[2]);
  }

  // Seek to first particle velocities and read triples
  skip = (3 * sizeof(float) * firstParticle);   // Particle 0 has no data
  inStream.seekg(skip, ios::cur);
  inStream.read(reinterpret_cast<char*>(fBlock),
                3 * numberOfParticles * sizeof(float));

  // Store the velocities in the double buffer in record order
  for (int p = 0; p < numberOfParticles; p++) {
    this->vx->push_back(fBlock[0]);
    this->vy->push_back(fBlock[1]);
    this->vz->push_back(fBlock[2]);
  }

  // Seek to first particle tags and read
  skip = (1 * sizeof(int) * firstParticle);     // skip to first tag
  inStream.seekg(skip, ios::cur);
  inStream.read(reinterpret_cast<char*>(iBlock),
                1 * numberOfParticles * sizeof(int));

  // Store the tags in the double buffer for sharing in record order
  for (int p = 0; p < numberOfParticles; p++) {
    this->tag->push_back(iBlock[0]);
    this->status->push_back(ALIVE);
    this->numberOfAliveParticles++;
    this->particleCount++;
  }
#ifdef DEBUG
    cout << "Rank " << setw(3) << this->myProc << " FINISHED READ " << endl;
#endif

  inStream.close();
  delete [] fBlock;
  delete [] iBlock;
}
