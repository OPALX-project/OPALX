#include <assert.h>
#include <complex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <stddef.h>

#include "distribution.h"

#ifndef USE_SLAB_WORKAROUND
#define USE_SLAB_WORKAROUND 1
#endif

enum {
    REDISTRIBUTE_1_TO_3,
    REDISTRIBUTE_3_TO_1
};


///
// return comma or period depending on position in a list
///
static inline char *separator(int i, int n)
{
    return i == (n - 1) ? "." : ", ";
}


///
// create 1-, 2- and 3-d cartesian data distributions
//   comm MPI Communicator
//   d1   1d data distribution
//   d2   2d data distribution
//   d3   3d data distribution
///
void distribution_init(MPI_Comm comm, const int n[], const int n_padded[], distribution_t *d, bool debug)
{
    int nproc;
    int self;
    int ndim = 3;
    int period[3];

    MPI_Comm_rank(comm, &self);
    MPI_Comm_size(comm, &nproc);

    d->debug = debug;
    for (int i = 0; i < 3; ++i) {
        d->n[i] = n[i];
        d->padding[i] = n_padded[i] - n[i];
    }

    // set up process grid with 1d decomposition
    d->process_topology_1.nproc[0] = 0;
    d->process_topology_1.nproc[1] = 1; // don't distribute outer dimensions
    d->process_topology_1.nproc[2] = 1; // don't distribute outer dimensions
    period[0] = period[1] = period[2] = 1;
    MPI_Dims_create(nproc, ndim, d->process_topology_1.nproc);
    if (d->debug && 0 == self) {
        fprintf(stderr, "Process grids:\n");
        fprintf(stderr, "  1d: ");
        for (int i = 0; i < ndim; ++i) {
            fprintf(stderr, "%d%s", d->process_topology_1.nproc[i], separator(i, ndim));
        }
        fprintf(stderr, "\n");
    }
    MPI_Cart_create(comm, ndim, d->process_topology_1.nproc, period, 0, &d->process_topology_1.cart);
    MPI_Cart_get(d->process_topology_1.cart, ndim, d->process_topology_1.nproc, d->process_topology_1.period, d->process_topology_1.self);
    d->process_topology_1.n[0] = n[0] / d->process_topology_1.nproc[0];
    d->process_topology_1.n[1] = n[1] / d->process_topology_1.nproc[1];
    d->process_topology_1.n[2] = n[2] / d->process_topology_1.nproc[2];

    // set up process grid with 2d decomposition
    d->process_topology_2.nproc[0] = 0;
    d->process_topology_2.nproc[1] = 0;
    d->process_topology_2.nproc[2] = 1; // don't distribute outer dimension
    period[0] = period[1] = period[2] = 1;
    MPI_Dims_create(nproc, ndim, d->process_topology_2.nproc);
    if (d->debug && 0 == self) {
        fprintf(stderr, "  2d: ");
        for (int i = 0; i < ndim; ++i) {
            fprintf(stderr, "%d%s", d->process_topology_2.nproc[i], separator(i, ndim));
        }
        fprintf(stderr, "\n");
    }
    MPI_Cart_create(comm, ndim, d->process_topology_2.nproc, period, 0, &d->process_topology_2.cart);
    MPI_Cart_get(d->process_topology_2.cart, ndim, d->process_topology_2.nproc, d->process_topology_2.period, d->process_topology_2.self);
    d->process_topology_2.n[0] = n[0] / d->process_topology_2.nproc[0];
    d->process_topology_2.n[1] = n[1] / d->process_topology_2.nproc[1];
    d->process_topology_2.n[2] = n[2] / d->process_topology_2.nproc[2];

    // set up process grid with 3d decomposition
    d->process_topology_3.nproc[0] = 0;
    d->process_topology_3.nproc[1] = 0;
    d->process_topology_3.nproc[2] = 0;
    period[0] = period[1] = period[2] = 1;
    MPI_Dims_create(nproc, ndim, d->process_topology_3.nproc);
    if (d->debug && 0 == self) {
        fprintf(stderr, "  3d: ");
        for (int i = 0; i < ndim; ++i) {
            fprintf(stderr, "%d%s", d->process_topology_3.nproc[i], separator(i, ndim));
        }
        fprintf(stderr, "\n");
    }
    MPI_Cart_create(comm, ndim, d->process_topology_3.nproc, period, 0, &d->process_topology_3.cart);
    MPI_Cart_get(d->process_topology_3.cart, ndim, d->process_topology_3.nproc, d->process_topology_3.period, d->process_topology_3.self);
    d->process_topology_3.n[0] = n[0] / d->process_topology_3.nproc[0];
    d->process_topology_3.n[1] = n[1] / d->process_topology_3.nproc[1];
    d->process_topology_3.n[2] = n[2] / d->process_topology_3.nproc[2];

    if (d->debug) {
        if (0 == self) {
            fprintf(stderr, "Process map:\n");
        }
        for (int p = 0; p < nproc; ++p) {
            MPI_Barrier(comm);
            if (p == self) {
                fprintf(stderr,
                        "  %d: 1d = (%d, %d, %d), 2d = (%d, %d, %d), 3d = (%d, %d, %d).\n",
                        self,
                        d->process_topology_1.self[0], d->process_topology_1.self[1], d->process_topology_1.self[2],
                        d->process_topology_2.self[0], d->process_topology_2.self[1], d->process_topology_2.self[2],
                        d->process_topology_3.self[0], d->process_topology_3.self[1], d->process_topology_3.self[2]);
            }
        }
    }
}


///
// clean up the data distribution
//   d    distribution descriptor
///
void distribution_fini(distribution_t *d)
{
    MPI_Comm_free(&d->process_topology_1.cart);
    MPI_Comm_free(&d->process_topology_2.cart);
    MPI_Comm_free(&d->process_topology_3.cart);
}


///
// check that the dimensions, n, of an array are commensurate with the
// process grids of this distribution
//   n    (global) grid dimensions
//   d    distribution descriptor
///
void distribution_assert_commensurate(distribution_t *d)
{
    for (int i = 0; i < 3; ++i) {
        assert(0 == (d->n[i] % d->process_topology_1.nproc[i]));
        assert(0 == (d->n[i] % d->process_topology_2.nproc[i]));
        assert(0 == (d->n[i] % d->process_topology_3.nproc[i]));
    }
}


// forward declarations
static void redistribute(const double complex *, double complex *, distribution_t *, int);
static void redistribute_slab(const double complex *, double complex *, distribution_t *, int);


///
// redistribute a 1-d to a 3-d data distribution
//   a    input
//   b    ouput
//   d    distribution descriptor
///
void distribution_1_to_3(const double complex *a,
                         double complex *b,
                         distribution_t *d)
{
    if (USE_SLAB_WORKAROUND) {
        redistribute_slab(a, b, d, REDISTRIBUTE_1_TO_3);
    } else {
        redistribute(a, b, d, REDISTRIBUTE_1_TO_3);
    }
}


///
// redistribute a 3-d to a 1-d data distribution
//   a    input
//   b    ouput
//   d    distribution descriptor
///
void distribution_3_to_1(const double complex *a,
                         double complex *b,
                         distribution_t *d)
{
    if (USE_SLAB_WORKAROUND) {
        redistribute_slab(a, b, d, REDISTRIBUTE_3_TO_1);
    } else {
        redistribute(a, b, d, REDISTRIBUTE_3_TO_1);
    }
}


///
// redistribute between 1- and 3-d distributions.
//   a    input
//   b    ouput
//   d    distribution descriptor
//   dir  direction of redistribution
//
// This actually does the work.
///
static void redistribute(const double complex *a,
                         double complex *b,
                         distribution_t *d,
                         int direction)
{
    int remaining_dim[3];
    MPI_Comm subgrid_cart;
    int subgrid_self;
    int subgrid_nproc;

    int self;

    MPI_Comm_rank(MPI_COMM_WORLD, &self);

    // exchange data with processes in a 2-d slab of 3-d subdomains

    remaining_dim[0] = 0;
    remaining_dim[1] = 1;
    remaining_dim[2] = 1;
    MPI_Cart_sub(d->process_topology_3.cart, remaining_dim, &subgrid_cart);
    MPI_Comm_rank(subgrid_cart, &subgrid_self);
    MPI_Comm_size(subgrid_cart, &subgrid_nproc);

    for (int p = 0; p < subgrid_nproc; ++p) {
        int d1_peer = (subgrid_self + p) % subgrid_nproc;
        int d3_peer = (subgrid_self - p + subgrid_nproc) % subgrid_nproc;
        int coord[2];
        int sizes[3];
        int subsizes[3];
        int starts[3];
        MPI_Datatype d1_type;
        MPI_Datatype d3_type;

        MPI_Cart_coords(subgrid_cart, d1_peer, 2, coord);
        if (0) {
            fprintf(stderr, "%d: d1_peer, d1_coord, d3_peer = %d, (%d, %d), %d\n",
                    self, d1_peer, coord[0], coord[1], d3_peer);
        }

        // create dataypes representing a subarray in the 1- and 3-d distributions

        sizes[0] = d->process_topology_1.n[0];
        sizes[1] = d->process_topology_1.n[1];
        sizes[2] = d->process_topology_1.n[2];
        subsizes[0] = d->process_topology_1.n[0];
        subsizes[1] = d->process_topology_3.n[1];
        subsizes[2] = d->process_topology_3.n[2];
        starts[0] = 0;
        starts[1] = coord[0] * d->process_topology_3.n[1];
        starts[2] = coord[1] * d->process_topology_3.n[2];
        MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE_COMPLEX, &d1_type);
        MPI_Type_commit(&d1_type);

        sizes[0] = d->process_topology_3.n[0];
        sizes[1] = d->process_topology_3.n[1];
        sizes[2] = d->process_topology_3.n[2];
        subsizes[0] = d->process_topology_1.n[0];
        subsizes[1] = d->process_topology_3.n[1];
        subsizes[2] = d->process_topology_3.n[2];
        starts[0] = d3_peer * d->process_topology_1.n[0];
        starts[1] = 0;
        starts[2] = 0;
        MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE_COMPLEX, &d3_type);
        MPI_Type_commit(&d3_type);

        // exchange data

        if (direction == REDISTRIBUTE_3_TO_1) {
            MPI_Sendrecv((void *) a, 1, d3_type, d3_peer, 0,
                         (void *) b, 1, d1_type, d1_peer, 0,
                         subgrid_cart, MPI_STATUS_IGNORE);
        } else if (direction == REDISTRIBUTE_1_TO_3) {
            MPI_Sendrecv((void *) a, 1, d1_type, d1_peer, 0,
                         (void *) b, 1, d3_type, d3_peer, 0,
                         subgrid_cart, MPI_STATUS_IGNORE);
        } else {
            abort();
        }

        // free datatypes

        MPI_Type_free(&d1_type);
        MPI_Type_free(&d3_type);
    }

    MPI_Comm_free(&subgrid_cart);
}


///
// redistribute between 1- and 3-d distributions.
//   a    input
//   b    ouput
//   d    distribution descriptor
//   dir  direction of redistribution
//
// This actually does the work, using slabs of subarrays to work
// around an issue in Open MPI with large non-contiguous datatypes.
///
static void redistribute_slab(const double complex *a,
                              double complex *b,
                              distribution_t *d,
                              int direction)
{
    int remaining_dim[3];
    MPI_Comm subgrid_cart;
    int subgrid_self;
    int subgrid_nproc;
    int self;
    ptrdiff_t d1_slice = d->process_topology_1.n[1] * d->process_topology_1.n[2] * sizeof(double complex);
    ptrdiff_t d3_slice = d->process_topology_3.n[1] * d->process_topology_3.n[2] * sizeof(double complex);

    MPI_Comm_rank(MPI_COMM_WORLD, &self);

    // exchange data with processes in a 2-d slab of 3-d subdomains

    remaining_dim[0] = 0;
    remaining_dim[1] = 1;
    remaining_dim[2] = 1;
    MPI_Cart_sub(d->process_topology_3.cart, remaining_dim, &subgrid_cart);
    MPI_Comm_rank(subgrid_cart, &subgrid_self);
    MPI_Comm_size(subgrid_cart, &subgrid_nproc);

    for (int p = 0; p < subgrid_nproc; ++p) {
        int coord[2];
        int d1_peer = (subgrid_self + p) % subgrid_nproc;
        int d3_peer = (subgrid_self - p + subgrid_nproc) % subgrid_nproc;

        MPI_Cart_coords(subgrid_cart, d1_peer, 2, coord);
        if (0) {
            fprintf(stderr, "%d: d1_peer, d1_coord, d3_peer = %d, (%d, %d), %d\n",
                    self, d1_peer, coord[0], coord[1], d3_peer);
        }

        for (int slice = 0; slice < d->process_topology_1.n[0]; ++slice) {
            int sizes[2];
            int subsizes[2];
            int starts[2];
            MPI_Datatype d1_type;
            MPI_Datatype d3_type;
            ptrdiff_t d1_offset = slice * d1_slice;
            ptrdiff_t d3_offset = (slice + d3_peer * d->process_topology_1.n[0]) * d3_slice;

            // create subarray dataypes representing the slice subarray in the 1- and 3-d distributions

            sizes[0] = d->process_topology_1.n[1];
            sizes[1] = d->process_topology_1.n[2];
            subsizes[0] = d->process_topology_3.n[1];
            subsizes[1] = d->process_topology_3.n[2];
            starts[0] = coord[0] * d->process_topology_3.n[1];
            starts[1] = coord[1] * d->process_topology_3.n[2];
            MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE_COMPLEX, &d1_type);
            MPI_Type_commit(&d1_type);

            MPI_Type_contiguous(d->process_topology_3.n[1] * d->process_topology_3.n[2],
                                MPI_DOUBLE_COMPLEX,
                                &d3_type);
            MPI_Type_commit(&d3_type);

            // exchange data

            if (direction == REDISTRIBUTE_3_TO_1) {
                MPI_Sendrecv((void *) a + d3_offset, 1, d3_type, d3_peer, 0,
                             (void *) b + d1_offset, 1, d1_type, d1_peer, 0,
                             subgrid_cart, MPI_STATUS_IGNORE);
            } else if (direction == REDISTRIBUTE_1_TO_3) {
                MPI_Sendrecv((void *) a + d1_offset, 1, d1_type, d1_peer, 0,
                             (void *) b + d3_offset, 1, d3_type, d3_peer, 0,
                             subgrid_cart, MPI_STATUS_IGNORE);
            } else {
                abort();
            }

            // free datatypes

            MPI_Type_free(&d1_type);
            MPI_Type_free(&d3_type);
        }
    }

    MPI_Comm_free(&subgrid_cart);
}

