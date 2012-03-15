/*=========================================================================

  Program:   Grid exchange driver program
  Module:    $RCSfile: GridTestP.cxx,v $

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// Test the parallel exchange of grids with neighbors which is not quite
// the standard ghost cell exchange.  In this case each processor has
// a contiguous array with an alive portion in the middle surrounded by
// a dead portion where the size on the front part of a dimension is not
// necessarily the size on the back.  When data is sent to a neighbor it
// comes from the alive part of the array, but when it is received it is
// unpacked into the dead part.
//

#include <iostream>
#include <iomanip>
#include "Partition.h"
#include "GridExchange.h"

#include <stdlib.h>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc != 4) {
    cout << "Usage: mpirun -np # GridTestP totalSize dead0 dead1" << endl;
  }

  // Size of the data array on a side
  int i = 1;
  int size = atoi(argv[i++]);

  // Size of dead regions on the front and back of each dimension
  int dead0 = atoi(argv[i++]);
  int dead1 = atoi(argv[i++]);

  // Initialize the partitioner which uses MPI Cartesian Topology
  //Partition::initialize(argc, argv);
  MPI_Init(&argc, &argv);
  Partition::initialize();
  int rank = Partition::getMyProc();

  // Allocate data to share
  int totalSize[DIMENSION];
  int dataSize = 1;

  for (int dim = 0; dim < DIMENSION; dim++) {
    totalSize[dim] = size;
    dataSize *= totalSize[dim];
  }

  GRID_T* data = new GRID_T[dataSize];
  for (int i = 0; i < dataSize; i++)
    data[i] = -1;

  for (int k = 0; k < totalSize[2]; k++) {
    for (int j = 0; j < totalSize[1]; j++) {
      for (int i = 0; i < totalSize[0]; i++) {
        if (k >= dead0 && k < (size - dead1) &&
            j >= dead0 && j < (size - dead1) &&
            i >= dead0 && i < (size - dead1)) {
              int index = (k * totalSize[0] * totalSize[1]) +
                          (j * totalSize[0]) + i;
              data[index] = rank;
        }
      }
    }
  }

//#ifdef DEBUG
  if (rank == 0) {
    cout << "BEFORE (Planes front to back)" << endl;
    for (int k = 0; k < totalSize[2]; k++) {
      cout << endl << " Plane " << k << endl << endl;
      for (int j = totalSize[1] - 1; j >= 0; j--) {
        for (int i = 0; i < totalSize[0]; i++) {
          int index = (k * totalSize[0] * totalSize[1]) +
                      (j * totalSize[0]) + i;
          cout << " " << setw(4) << data[index];
        }
        cout << endl;
      }
      cout << endl << endl << endl;
    }
  }
//#endif

  // Construct the grid exchanger
  GridExchange exchange(totalSize, dead0, dead1);

  // Exchange grid
  exchange.exchangeGrid(data);

//#ifdef DEBUG
  if (rank == 0) {
    cout << "AFTER (Planes front to back)" << endl;
    for (int k = 0; k < totalSize[2]; k++) {
      cout << endl << " Plane " << k << endl << endl;
      for (int j = totalSize[1] - 1; j >= 0; j--) {
        for (int i = 0; i < totalSize[0]; i++) {
          int index = (k * totalSize[0] * totalSize[1]) +
                      (j * totalSize[0]) + i;
          cout << " " << setw(4) << data[index];
        }
        cout << endl;
      }
      cout << endl << endl << endl;
    }
  }
//#endif

  // Shut down MPI
  Partition::finalize();
  MPI_Finalize();

  return 0;
}
