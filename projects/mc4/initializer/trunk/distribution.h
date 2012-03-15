#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#ifdef __cplusplus
extern "C" {

#ifndef NOFFTW

#ifdef FFTW2
    #include <fftw.h>
#endif
#ifdef FFTW3
    #include <fftw3-mpi.h>
#endif
typedef my_fftw_complex complex_t;
#else
    #include "TypesAndDefs.h"
    typedef my_fftw_complex complex_t;
#endif

#else

#include <complex.h>
typedef double complex complex_t;
#endif

#include <mpi.h>


///
// descriptor for a process grid
//   cart     Cartesian MPI communicator
//   nproc[]  dimensions of process grid
//   period[] periods of process grid
//   self[]   coordinate of this process in the process grid
//   n[]      local grid dimensions
///
typedef struct {
    MPI_Comm cart;
    int nproc[3];
    int period[3];
    int self[3];
    int n[3];
} process_topology_t;


///
// descriptor for data distribution
//   debug               toggle debug output
//   n[3]                (global) grid dimensions
//   padding[3]          padding applied to (local) arrays
//   process_topology_1  1-d process topology
//   process_topology_2  2-d process topology
//   process_topology_3  3-d process topology
///
typedef struct {
    bool debug;
    int n[3];
    int padding[3];
    process_topology_t process_topology_1;
    process_topology_t process_topology_2;
    process_topology_t process_topology_3;
} distribution_t;


///
// create 1-, 2- and 3-d cartesian data distributions
//   comm   MPI Communicator
//   d      distribution descriptor
//   n      (global) grid dimensions (3 element array)
//   debug  debugging output
///
void distribution_init(MPI_Comm comm, const int n[], const int n_padded[], distribution_t *d, bool debug);


///
// clean up the data distribution
//   d    distribution descriptor
///
void distribution_fini(distribution_t *d);


///
// assert that the data and processor grids are commensurate
//   d    distribution descriptor
///
void distribution_assert_commensurate(distribution_t *d);


///
// redistribute a 1-d to a 3-d data distribution
//   a    input
//   b    ouput
//   d    distribution descriptor
///
void distribution_1_to_3(const complex_t *a,
                         complex_t *b,
                         distribution_t *d);

///
// redistribute a 3-d to a 1-d data distribution
//   a    input
//   b    ouput
//   d    distribution descriptor
///
void distribution_3_to_1(const complex_t *a,
                         complex_t *b,
                         distribution_t *d);


static inline int distribution_get_nproc_1d(distribution_t *d, int direction)
{
    return d->process_topology_1.nproc[direction];
}

static inline int distribution_get_nproc_2d(distribution_t *d, int direction)
{
    return d->process_topology_2.nproc[direction];
}

static inline int distribution_get_nproc_3d(distribution_t *d, int direction)
{
    return d->process_topology_3.nproc[direction];
}

static inline int distribution_get_self_1d(distribution_t *d, int direction)
{
    return d->process_topology_1.self[direction];
}

static inline int distribution_get_self_2d(distribution_t *d, int direction)
{
    return d->process_topology_2.self[direction];
}

static inline int distribution_get_self_3d(distribution_t *d, int direction)
{
    return d->process_topology_3.nproc[direction];
}

#ifdef __cplusplus
}
#endif

#endif
