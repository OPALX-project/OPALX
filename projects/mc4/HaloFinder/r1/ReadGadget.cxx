#include <fstream>
#include <iostream>

using namespace std;

#define WORDSIZE 8

/////////////////////////////////////////////////////////////////////////
//
// Read in the number of items from the file pointer and
// byte swap if necessary
//
/////////////////////////////////////////////////////////////////////////

void readData(
        bool swapBytes,
        long int* data,
        unsigned long dataSize,
        unsigned long dataCount,
        FILE* fp)
{
   // Read all the data from the file
   fread(data, dataSize, dataCount, fp);

   if (swapBytes == true) {

      // Byte swap each INTEGER
      char* dataPtr = (char*) data;
      char temp;
      for (int item = 0; item < dataCount; item++) {

         // Do a byte-by-byte swap, reversing the order.
         for (int i = 0; i < dataSize / 2; i++) {
            temp = dataPtr[i];
            dataPtr[i] = dataPtr[dataSize - 1 - i];
            dataPtr[dataSize - 1 - i] = temp;
         }
         dataPtr += WORDSIZE;
      }
   }
}

void readData(
        bool swapBytes,
        int* data,
        unsigned long dataSize,
        unsigned long dataCount,
        FILE* fp)
{
   // Read all the data from the file
   fread(data, dataSize, dataCount, fp);

   if (swapBytes == true) {

      // Byte swap each INTEGER
      char* dataPtr = (char*) data;
      char temp;
      for (int item = 0; item < dataCount; item++) {

         // Do a byte-by-byte swap, reversing the order.
         for (int i = 0; i < dataSize / 2; i++) {
            temp = dataPtr[i];
            dataPtr[i] = dataPtr[dataSize - 1 - i];
            dataPtr[dataSize - 1 - i] = temp;
         }
         dataPtr += WORDSIZE;
      }
   }
}

void readData(
        bool swapBytes,
        short int* data,
        unsigned long dataSize,
        unsigned long dataCount,
        FILE* fp)
{
   // Read all the data from the file
   fread(data, dataSize, dataCount, fp);

   if (swapBytes == true) {

      // Byte swap each INTEGER
      char* dataPtr = (char*) data;
      char temp;
      for (int item = 0; item < dataCount; item++) {

         // Do a byte-by-byte swap, reversing the order.
         for (int i = 0; i < dataSize / 2; i++) {
            temp = dataPtr[i];
            dataPtr[i] = dataPtr[dataSize - 1 - i];
            dataPtr[dataSize - 1 - i] = temp;
         }
         dataPtr += WORDSIZE;
      }
   }
}

void readData(
        bool swapBytes,
        float* data,
        unsigned long dataSize,
        unsigned long dataCount,
        FILE* fp)
{
   // Read all the data from the file
   fread(data, dataSize, dataCount, fp);

   if (swapBytes == true) {

      // Byte swap each FLOAT
      char* dataPtr = (char*) data;
      char temp;
      for (int item = 0; item < dataCount; item++) {

         // Do a byte-by-byte swap, reversing the order.
         for (int i = 0; i < dataSize / 2; i++) {
            temp = dataPtr[i];
            dataPtr[i] = dataPtr[dataSize - 1 - i];
            dataPtr[dataSize - 1 - i] = temp;
         }
         dataPtr += WORDSIZE;
      }
   }
}

/////////////////////////////////////////////////////////////////////////
//
//
/////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  bool swapBytes = false;                // Bytes must be swapped

  string fileName = argv[1];

  struct io_header_1 {
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      flag_what;
    int      num_files;
    float   BoxSize;
    float   BoxSize1;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    int      flag_multiphase;
    int      flag_stellarage;
    int      flag_sfrhistogram;
    char     fill[76];  /* fills to 256 Bytes */
  } header;
cout << "Size of header " << sizeof(header) << endl;

  FILE *fp = fopen(fileName.c_str(), "r");
  fread(&header, sizeof(header), 1, fp);

  for (int i = 0; i < 6; i++)
    cout << "Header: npart " << header.npart[i] << endl;
  for (int i = 0; i < 6; i++)
    cout << "Header: mass " << header.mass[i] << endl;
  cout << "Header: time " << header.time << endl;
  cout << "Header: redshift " << header.redshift << endl;
  cout << "Header: flag_sfr " << header.flag_sfr << endl;
  cout << "Header: flag_feedback " << header.flag_feedback << endl;
  for (int i = 0; i < 6; i++)
    cout << "Header: npartTotal " << header.npartTotal[i] << endl;
  cout << "Header: flag_cooling " << header.flag_cooling << endl;
  cout << "Header: flag_what " << header.flag_what << endl;
  cout << "Header: num_files " << header.num_files << endl;
  cout << "Header: BoxSize " << header.BoxSize << endl;
  cout << "Header: Omega0 " << header.Omega0 << endl;
  cout << "Header: OmegaLambda " << header.OmegaLambda << endl;
  cout << "Header: HubbleParam " << header.HubbleParam << endl;
  cout << "Header: flag_multiphase " << header.flag_multiphase << endl;
  cout << "Header: flag_stellarage " << header.flag_stellarage << endl;
  cout << "Header: flag_sfrhistogram " << header.flag_sfrhistogram << endl;

  // Particle 0 is written but does not contain valid data
  int numberOfParticles = header.npart[2];
  int totalParticles = header.npartTotal[2];
  int numberOfFiles = header.num_files;

  cout << "Header: Particles this file " << numberOfParticles << endl;
  cout << "Header: Particles total " << totalParticles << endl;
  cout << "Header: Number of files " << numberOfFiles << endl;

  //
  // Read locations
  //

  // Skip array position 0 because particle tags start with 1
  float* loc = new float[3];
  int p = 0;
  fread(loc, sizeof(float), 3, fp);

  // Set initial min and max locations
  float minLoc[3], maxLoc[3];
  p = 1;
  fread(loc, sizeof(float), 3, fp);
  cout << p << " " << loc[0] << "," << loc[1] << "," << loc[2] << endl;
  for (int dim = 0; dim < 3; dim++)
    minLoc[dim] = maxLoc[dim] = loc[dim];

  p = 2;
  for (; p <= numberOfParticles; p++) {

    fread(loc, sizeof(float), 3, fp);

    for (int dim = 0; dim < 3; dim++) {
      if (minLoc[dim] > loc[dim])
        minLoc[dim] = loc[dim];
      if (maxLoc[dim] < loc[dim])
        maxLoc[dim] = loc[dim];
    }
    if (p < 5 || p > (numberOfParticles - 5))
      cout << p << " " << loc[0] << "," << loc[1] << "," << loc[2] << endl;
  }

  // Range of locations in file
  for (int dim = 0; dim < 3; dim++) {
    cout << "LOCATION RANGE Dimension " << dim << " [" 
         << minLoc[dim] << ":" << maxLoc[dim] << "]" << endl;
  }

  //
  // Read velocities
  //

  // Skip array position 0 because particle tags start with 1
  float* vel = new float[3];
  p = 0;
  fread(vel, sizeof(float), 3, fp);

  // Set initial min and max velocities
  float minVel[3], maxVel[3];
  p = 1;
  fread(vel, sizeof(float), 3, fp);
  cout << p << " " << vel[0] << "," << vel[1] << "," << vel[2] << endl;
  for (int dim = 0; dim < 3; dim++)
    minVel[dim] = maxVel[dim] = vel[dim];

  p = 2;
  for (; p <= numberOfParticles; p++) {

    fread(vel, sizeof(float), 3, fp);

    for (int dim = 0; dim < 3; dim++) {
      if (minVel[dim] > vel[dim])
        minVel[dim] = vel[dim];
      if (maxVel[dim] < vel[dim])
        maxVel[dim] = vel[dim];
    }
    if (p < 5 || p > (numberOfParticles - 5))
      cout << p << " " << vel[0] << "," << vel[1] << "," << vel[2] << endl;
  }

  // Range of velocities in file
  for (int dim = 0; dim < 3; dim++) {
    cout << "VELOCITY RANGE Dimension " << dim 
         << " [" << minVel[dim] << ":" << maxVel[dim] << "]" << endl;
  }

  //
  // Read tags
  //

  // Skip array position 0 because tags start with 1
  int xtag;
  fread(&xtag, sizeof(int), 1, fp);

  int* tag = new int[numberOfParticles];
  fread(tag, sizeof(int), numberOfParticles, fp);
  p = 0;

  int minTag = totalParticles + 1;
  int maxTag = -1;
  for (; p < numberOfParticles; p++) {
  
    if (minTag > tag[p])
      minTag = tag[p];
    if (maxTag < tag[p])
      maxTag = tag[p];

    if (p < 30 || p > (numberOfParticles - 10))
      cout << p << " " << tag[p] << endl;
  }

  if (maxTag >= totalParticles)
    cout << "Tag exceeded " << maxTag << endl;
  if (minTag < 1)
    cout << "Tag exceeded " << minTag << endl;

int x;
for (int i = 0; i < 20; i++) {
  fread(&x, sizeof(int), 1, fp);
  cout << i << " X " << x << endl;
}

  fclose(fp);
}

