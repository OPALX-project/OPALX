#ifndef __HWLOC_COMM_TOPOLOGY__
#define __HWLOC_COMM_TOPOLOGY__

#include "mpi.h"

#include "Comm/Topology/CommTopology.h"

#include <hwloc.h>


class HwLocCommTopology : protected CommTopology {

public:

    HwLocCommTopology(MPI_Comm comm = MPI_COMM_WORLD) : CommTopology(comm)
    {}

    virtual ~HwLocCommTopology()
    {}

    void discover();

private:

    int num_nodes_;
    int num_cores_;
    int num_pus_;

    void printTopologyInformation(hwloc_topology_t topology);
    void printChildren(hwloc_topology_t topology, hwloc_obj_t obj, int depth);
};

#endif
